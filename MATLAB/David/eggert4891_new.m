%For people using octave - need to install mapping, signal, and communications
%% Constants
%Sampling frequency
j=1i;
close all;
Fs = 5000000;

%Number of samples
N=50000;

%Time for samples
Tmax = N/Fs;

%% Reading Scope Data
[x,y] = importAgilentBin("scope_4_1.bin");
x = transpose(x);
y = transpose(y);

% Make time start at 0
x = x+Tmax/2;
figure;
plot (x,y);
legend("Scope Signal");
figure;

%Rough frequency detection
numSamples = 256;
trunced_peaks = abs(fft(y(1:numSamples).^2 ));

plot(trunced_peaks);
legend("FFT");
figure;

%We know our peak will be around 400
%zero out values not close to 400
numdel = round(numSamples/3);
trunced_peaks(1:numdel) = zeros(1,numdel);
trunced_peaks(numSamples-numdel+1:numSamples) = zeros(1,numdel);

plot(trunced_peaks);
legend("FFT");
figure;
[m,fc_actual] = max(trunced_peaks);
%parabolic interpolation for FFTs
d= 0.5*(trunced_peaks(fc_actual-1)-trunced_peaks(fc_actual+1))/(trunced_peaks(fc_actual-1)-2*trunced_peaks(fc_actual)+trunced_peaks(fc_actual+1));
fc_actual = fc_actual + d-1;
fc_actual = fc_actual* Fs/(2*numSamples);
%Frequency c

fc_received = fc_actual;
%Get I and Q
i = y .* cos(2*pi*fc_received * x);
q = y.* sin(2*pi*fc_received * x );
l = i + j*q;

% Downsample by 20

l = decimate(l, 20);
N = length(l);
Fs = Fs / 20;
fc_received = fc_received/20;
%% Filtering
%SRRC filter
sps = ceil(Fs/fc_received);
span = 5;
srrc = rcosdesign(0.5, span, sps, 'sqrt');
l = upfirdn(l, srrc,1,1);
l = l(length(srrc)/2+0.5:end-length(srrc)/2+0.5);

%Low pass filter
%lp_filter = fir1(34, 0.5*fc_received/Fs);
%l = upfirdn(l, lp_filter,1,1);
%l = l(18:end-17);

%Recover I and Q from L = lowpass signal
i = real(l);
q = imag(l);
%i = max(min(0.0001, i), -0.0001)* 3000;
%q = max(min(0.0001, q), -0.0001)* 3000;
hold on;
plot(i);
plot(q);
legend("I - not phase corrected", "Q - not phase corrected");
figure;
%% Phase correction
%Get a good estimate of initial phase
freq_log = zeros(1, N);
phase_log = zeros(1, N);
err_log = zeros(1, N);
freq=0;
%err(k) does not need to be a list - this is just
%here for tracking purposes.
err =zeros(1,N);
A = 0.575*100;
B = 0.255*100;
phase = atan2(imag(l(1)), real(l(1)));
for k = 1:N
  %Adjust phase
  l(k) = l(k) * exp(-j*phase);
  i(k) = real(l(k));
  q(k) = imag(l(k));
  err(k) = i(k)*q(k);
  err_log(k) = err(k);
  freq = freq + B*err(k);
  freq_log(k) = freq*Fs/(2*pi);
  phase = phase + freq + A*err(k);
  phase_log(k) = phase;
end

plot(freq_log);
legend("Frequency offset");
figure;
hold on;
plot(i);
plot(q);
legend("I - phase corrected", "Q - phase corrected");
figure;
plot(movmean(err, 2));
legend("error")
figure;

%% Clock recovery
%create zero crossing signal
zero_cross = sign(i(1:end-1).*i(2:end));
zero_cross = zero_cross - mean(zero_cross);

%Get data rate
zero_fft = abs(fft(zero_cross));
[m, max_ind] = max(zero_fft(1:N/sps*1.5));
sps = N / max_ind;
plot(zero_fft);
legend("FFT - zero crossing");
figure;

%PLL clock recovery
Kp_PLL = 0.1;
Ki_PLL = 0.01;
freq = sps;
phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
real_samp = [];
quad_samp=[];
err = [];
integrator = freq;
for k = 1:N-1
  %Check if we have crossed zero on the phase accumulator
  %If the previous value was around -pi and the current value is around pi - we crossed it.
  if (wrapToPi(phase_acc )< -pi * 0.33 && wrapToPi(prev_acc) > pi*0.33)
    %we have crossed zero - sample
    real_samp(end+1) = i(k);
    quad_samp(end+1) = q(k);
  end
  %Detect zero crossing and use it to tune PLL
  if (zero_cross(k) < -1)
    %we have a zero crossing - adjust freq
    err(end+1) = wrapToPi(phase_acc);
    integrator = integrator + err(end)*Ki_PLL;
    freq = err(end) * Kp_PLL + integrator;
  end
  prev_acc = phase_acc;
  phase_acc = phase_acc + 2*pi/freq;
end
hold on;
plot (real_samp);
plot(quad_samp);
selected = real_samp +j*quad_samp;
legend("Real sampled", "Imaginary sampled");
figure;
scatterplot(real_samp+quad_samp*j);
figure;
plot(err);
legend("PLL error");
figure;

%% Bit conversion
%Conversion to bits
selected_sign = abs(real(selected))./real(selected);
bits_received = double(selected_sign > 0);
key2 = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1,1];
[bit_corr, bit_pos] = xcorr( (bits_received -0.5)*2, key2);
plot(bit_pos, bit_corr);
legend("Bit correlation to key");
shift = 0;

for i=1:size(bit_corr,2)
    if abs(bit_corr(i) + bit_corr(i+15) +bit_corr(i+30) ) > 39

        shift = bit_pos(i)+45;
        if (bit_corr(i) < 0)
            bits_received = bits_received * -1 + 1;
        end
        break;
    end
end
received = bits_received(shift+1:shift+256);
strarr = int2str(received);
disp(strarr);



