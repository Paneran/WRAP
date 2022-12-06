#For people using octave - need to install mapping, signal, and communications
%% Constants
#Sampling frequency
close all;
Fs = 5000000;

#Number of samples
N=50000;

#Time for samples
Tmax = N/Fs;

%% Reading Scope Data
[x,y] = importAgilentBin("scope_4.bin", '1');
x = transpose(x);
y = transpose(y);
# Make time start at 0
x += Tmax/2;
figure;
plot (x,y);
legend("Scope Signal");
figure;

%% Frequency correction

#Detect frequency with FFT
trunced_peaks = abs(fft(y.^2));

#Zero out first few FFT values
trunced_peaks(1:5) = zeros(1,5);
plot(trunced_peaks);
legend("FFT");
figure;
[_,fc_actual] = max(trunced_peaks);
fc_actual *= Fs/(2*N);
#Frequency correction

fc_received = 985850;
#Get I and Q
i = y .* cos(2*pi*fc_received * x);
q = y.* sin(2*pi*fc_received * x );
l = i + j*q;

%% Filtering
#SRRC filter
sps =Fs/fc_received*20;
span = 5;
srrc = rcosdesign(0.5, span, sps, 'sqrt');
l = upfirdn(l, srrc,1,1);
l = l(length(srrc)/2+0.5:end-length(srrc)/2+0.5);

#Low pass filter
%lp_filter = fir1(34, 0.5*fc_received/Fs);
%l = upfirdn(l, lp_filter,1,1);
%l = l(18:end-17);

#Recover I and Q from L = lowpass signal
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
#Get a good estimate of initial phase
phase_acc = zeros(1, N+1);
phase_acc(1) = atan2(q(1), i(1));
#err(k) does not need to be a list - this is just
#here for tracking purposes.
err =zeros(1,N);
#our loop filter has transfer function:
#(0.000439x^2 + 0.000878z + 0.000439)
# ------------------------------------------------------
#(z^2 - 1.951z + 0.9512)
px = 0;
px2 = 0;
py =phase_acc(1) ;
py2 = phase_acc(1) ;
for k = 1:N
  #Adjust phase
  i(k) = real(l(k) * exp(-j*phase_acc(k)));
  q(k) = imag(l(k) * exp(-j*phase_acc(k)));
  err(k) = i(k)*q(k);
  #Implement loop filter
  phase_acc(k+1) = 10*(0.00086024*err(k) + 0.0017*px +0.00086024*px2) +1.9785*py -0.9785*py2 ;

  px2 = px;
  px = err(k);
  py2 = py;
  py = phase_acc(k+1);


end
plot(phase_acc);
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
#create zero crossing signal
zero_cross = sign(i(1:end-1).*i(2:end));
zero_cross = zero_cross - mean(zero_cross);

#Get data rate
zero_fft = abs(fft(zero_cross));
[max, max_ind] = max(zero_fft(1:N/sps*1.5));
sps = N / max_ind;
plot(zero_fft);
legend("FFT - zero crossing");
figure;

#PLL clock recovery
Kp_PLL = 5;
Ki_PLL = 1;
freq = sps;
phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
real_samp = [];
quad_samp=[];
err = [];
integrator = freq;
for k = 1:N-1
  #Check if we have crossed zero on the phase accumulator
  #If the previous value was around -pi and the current value is around pi - we crossed it.
  if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
    #we have crossed zero - sample
    real_samp(end+1) = i(k);
    quad_samp(end+1) = q(k);
  endif
  #Detect zero crossing and use it to tune PLL
  if (zero_cross(k) < -1)
    #we have a zero crossing - adjust freq
    err(end+1) = wrapToPi(phase_acc);
    integrator += err(end)*Ki_PLL;
    freq = err(end) * Kp_PLL + integrator;
  endif
  prev_acc = phase_acc;
  phase_acc += 2*pi/freq;
end
hold on;
plot (real_samp);
plot(quad_samp);
selected = real_samp +j*quad_samp;
legend("Real sampled", "Imaginary sampled");

scatterplot(real_samp+quad_samp*j);
figure;
plot(err);
legend("PLL error");
figure;

%% Bit conversion
#Conversion to bits
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



