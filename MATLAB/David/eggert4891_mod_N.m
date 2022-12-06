% Method:
% -Find the frequency of the carrier signal. 
% -Convert to the in_ph and quad components. 
% -Filter with an SRRC. 
% -Phase correct. Assumes phase shift on each point is different. Fixes
% phase point by point, with a loop filter that stresses smoothness. Some
% digital filter, wonder how he came up with those constants
% -Clock Recovery: Find the correct sampling clock
% watch: https://www.youtube.com/watch?v=HcYFFlsSLrg
%   -Find the symbols per second: Find rising edge of each signal. Take FFT
%   ==> gives the symbol rate. Not sure what david is doing to make this
%   happen/ how it makes sense. Left commented code on how I would do it.
%   Though it gets different results. 
% -PLL:
%   -If there was about a 180 change in phase, thats a bit of interest (not
%   sure if this is what we should be looking for)
%   -Find error and implement loop filter to ensure stability. Find the
%   frequncy. 
% -Find packet header, and extract bits. 
% 
% The code is working. if you take the output and put it into a binary to
% ascii converter, it says "I don't know you tell me eggert!"

% I think this is working. fixed everything. Things that need review/
% thought:
% 1. fc_actual is the peak frequency ~1MHz. I assume this is the carrier
% frequency, not the symbol rate. so sps = Fs/fc_actual*20, where 20 is
% transmitted sps. 
% 2. Where did david find loop filter for phase correction
% 4. If you tune Kp_PLL and Ki_PLL, it affects the convergence of PLL and
% can result in vectors of different sizes. 
% 5. PLL also is not converging. I feel something deeper is wrong. See
% Figure 9

%% Constants
% Sampling frequency
j = 1i;
close all;
clear;
Fs = 5000000;

% Number of samples
N=50000;

% Time for samples
Tmax = N/Fs;

%% Reading Scope Data
[x,y] = importAgilentBin('scope_4_1.bin');
x = transpose(x);
y = transpose(y);
% Make time start at 0
x = x+Tmax/2;
figure;
plot (x,y);
legend("Scope Signal");

%% Frequency correction
% Detect frequency with FFT
trunced_peaks = abs(fft(y.^2));

% Zero out first few FFT values
trunced_peaks(1:5) = zeros(1,5);

L = length(trunced_peaks);
% Find the double sided spectrum
P2 = abs(trunced_peaks/L);
% Plot the single_sided
P1 = P2(1:ceil(L/2));
P1(2:end-1) = 2*P1(2:end-1);
% make frequency scaling
fn = Fs*(0:(L+1)/2-1)/(L+1);
figure;
plot(fn,P1)
legend("FFT");
[M,I] = max(trunced_peaks);
fc_actual = fn(I)/2;


% Frequency correction
% Get I and Q
in_ph = y .* cos(2*pi*fc_actual * x);
q = -y.* sin(2*pi*fc_actual * x );
l = in_ph + j*q;

%% Filtering
% SRRC filter
sps = floor(Fs/fc_actual)*20;
disp(sps)
span = 5;
srrc = rcosdesign(0.5, span, sps, 'sqrt');
l = upfirdn(l, srrc,1,1);
l = l(length(srrc)/2+0.5:end-length(srrc)/2+0.5);

% Low pass filter
%lp_filter = fir1(34, 0.5*fc_received/Fs);
%l = upfirdn(l, lp_filter,1,1);
%l = l(18:end-17);

% Recover I and Q from L = lowpass signal
in_ph = real(l);
q = imag(l);
% i = max(min(0.0001, i), -0.0001)* 3000;
% q = max(min(0.0001, q), -0.0001)* 3000;
hold on;
plot(in_ph);
plot(q);
legend("I - not phase corrected", "Q - not phase corrected");
figure;

%% Phase correction
% Get a good estimate of initial phase
phase_acc = zeros(1, N+1);
phase_acc(1) = atan2(q(1), in_ph(1));
% err(k) does not need to be a list - this is just
% here for tracking purposes.
err =zeros(1,N);
% our loop filter has transfer function:
% (0.000439x^2 + 0.000878z + 0.000439)
% ------------------------------------------------------
% (z^2 - 1.951z + 0.9512)
px = 0;
px2 = 0;
py =phase_acc(1) ;
py2 = phase_acc(1) ;
for k = 1:N
  % Adjust phase
  l(k) = l(k) * exp(-j*phase_acc(k));
  in_ph(k) = real(l(k));
  q(k) = imag(l(k));
  err(k) = in_ph(k)*q(k);
  % Implement loop filter
  phase_acc(k+1) = 10*(0.00086024*err(k) + 0.0017*px +0.00086024*px2) +1.9785*py -0.9785*py2 ;

  px2 = px;
  px = err(k);
  py2 = py;
  py = phase_acc(k+1);
end
plot(phase_acc);
figure;
hold on;
plot(in_ph);
plot(q);
legend("I - phase corrected", "Q - phase corrected");
figure;
plot(movmean(err, 2));
legend("error")
figure;


%% Clock recovery
% create zero crossing signal
zero_cross = sign(in_ph(1:end-1).*in_ph(2:end));

% Get bit rate
zero_fft = abs(fft(zero_cross));
zero_fft(1) = 0;
L = length(zero_cross);
% make frequency scaling
fn = Fs*(0:(L-1)/2)/L;
% given non info points, I = 1 will be highest. exclude it. 
[peaks, loc] = findpeaks(zero_fft(1:end), 'MinPeakHeight', 150);
Rsr = fn(loc(1));
spsr = Fs/Rsr;
plot(zero_fft)
legend("FFT - zero crossing");
figure;

zero_cross = zero_cross - mean(zero_cross);

% Get data rate
zero_fft = abs(fft(zero_cross));
[peaks, loc] = findpeaks(zero_fft(1:end), 'MinPeakHeight', 150);
sps = N / loc(1);
plot(zero_fft);
legend("FFT - zero crossing_2");
figure;

disp(sps)
disp(spsr)
sps = spsr;

% PLL clock recovery
Kp_PLL = 1;
Ki_PLL = 5;
freq = sps;
phase_acc = pi+0.00001;
prev_acc = pi-0.00001;
real_samp = [];
quad_samp=[];
err = [];
integrator = freq;
for k = 1:N-1
    % Check if we have crossed zero on the phase accumulator
    % If the previous value was around -pi and the current value is around pi - we crossed it.
    if (wrapToPi(phase_acc )< -pi * 0.75 && wrapToPi(prev_acc) > pi*0.75)
        % we have crossed zero - sample
        real_samp(end+1) = in_ph(k);
        quad_samp(end+1) = q(k);
    end
    % Detect zero crossing and use it to tune PLL
    if (zero_cross(k) < -1)
        % we have a zero crossing - adjust freq
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

scatterplot(real_samp+quad_samp*j);
figure;
plot(err);
legend("PLL error");
figure;

%% Bit conversion
% Conversion to bits
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
len = shift + 256;
if len > length(bits_received)
    len = length(bits_received);
end
received = bits_received(shift+1:len);
strarr = int2str(received);
disp(strarr);
