clear;
%% Setup
Fs = 1000; % Sampling rate
Tmax = 10; % Max time (the last sample is at t = Tmax-1/Fs)
N = Tmax * Fs; % Number of samples

% Use these guys for plotting
t = linspace(0, Tmax-1/Fs, N); % Time variable
f = linspace(-Fs/2,Fs/2, N); % Frequency variable

Rs = 50; % Symbol rate in symbols/sec (baud)
Ns = Rs * Tmax; % Number of symbols to send
sps = Fs / Rs; % Number of samples per symbol

fc = 100; % Carrier frequency in Hz

% RRC setup
span = 100;
len = sps * span;
rolloff = 0.5;
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

%% Symbol Generation
% Generate random set of bits (0s and 1s)
bits = randi([0, 1], 1, Ns/2);

% Map bits to BPSK symbols (1s and -1s)
second_half = 2 * (bits - 0.5);
key = [1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1]; % length = 15
header = [key,key,key]; % length 45
header_bits = header == 1;
Tx_bits = [header_bits,bits];
quiet = zeros(1,Ns/2-length(header)); %length 250 - 45 = 205
first_half = [quiet,header];
symbs = [first_half,second_half];

% Plot BPSK constellation digram
%%scatterplot(symbs);

% The above should help you get started. Continue doing the rest on your
% own
% The first half of the bit sequence is the packet header
%% Shifted Deltas
deltas = upsample(symbs,sps);

%% Pulse Shaping (RRC)
symbol_wave = upfirdn(deltas, RRC, 1, 1);
%symbol_wave = conv(deltas, RRC); This has the same functionality as the 
%                                 upfirdn with up-sampling and
%                                 down-sampling parameters of 1.
symbol_wave = symbol_wave(1001:length(symbol_wave)-1000);
% symbol_wave deltas has 10000 samples and RRC has 2000. The result of 
% their convolution will have 3000 samples but we want to start at the 
% point where the peak of RRC aligns with the first delta. So the first
% 1000 samples are not useful. Similarly, we want to stop when the peak of
% of RRC aligns with the last delta. The beauty of RRC is that when these 
% other  allignments happen for each delta, the value of RRC is zero at all 
% the deltas. We do not care about the cases when the peak is not aligned 
% with any deltas because when we sample at the receiver, these values will
% be discarded. 


figure('Name','Transmitter Baseband Representation','NumberTitle','off');
stem(t,deltas);
xlabel('Time(s)');
ylabel('Symbols');
hold on
plot(t,symbol_wave);
xlabel('Time(s)');
ylabel('Symbol Waveform');
hold off
%}

%% Carrier Modulation
comp_carr = exp(1i*2*pi*fc*t);
analytic_sig = symbol_wave .* comp_carr;
Tx = real(analytic_sig);

%figure('Name','Transmitter Spectrum','NumberTitle','off');
%plot(f,fftshift(fft(Tx))); 

% The fft function considers the zero frequency to be at the left instead 
% of center. To set the zero frequency at the center, we have to wrap the
% fft output around itself i.e. consider the second half to be the negative 
% frequencies. We do this wrapping by using the fftshift function
%% Noise
std = 0.1;
noise = std * randn(1,N);
Rx = Tx + noise;
%Rx = Tx;

%% Carrier Demodulation
%figure('Name','Receiver Spectrum','NumberTitle','off');
%plot(f,fftshift(fft(Rx.^2)));
%plot(f,(fft(Rx.^2)));

% Frequency Correction
%We suqre Rx so that its frequency domain gets convolved with itself, As a
%result, instead of having two peaks at fc and -fc, we will have peaks at
%the frequencies -2fc, 0, 2fc. We have generated a peak at DC which will
%help us calculate the frequency offset without knowing fc!!!
%[peaks,loc] = findpeaks(real(fft(Rx.^2)),'MinPeakHeight',50);
%fcr = (loc(1)-1)/2;
phase = unifrnd(-pi/2,pi/2);

adj_fft = fftshift(fft(Rx.^2));
[peaks,loc] = findpeaks(abs(adj_fft),'MinPeakHeight',20);
fcr = floor(f(loc(end))/2);

yQ = Rx .* (-sin(2 * pi * fcr * t + phase));
yI = Rx .* cos(2 * pi * fcr * t + phase);
% For simplicity we assume the receiver inherently generates the phase
% offset while processing the data. This is equivalent of saying the
% channel changes the phases and the receiver is ideal. 

% -----------------------------------------------------------------------
% Applying Channel Non-idealities:
% We are applying the non-idealities to the baseband representation since
% applying the non-idealities at the bandpass is more complex and still
% gives us the same result. Remember, bad stuff could have happened around
% fc and the resulting destorted data could have been brought to baseband.
% Or, we could bring the data to baseband and then apply non-idealities to
% generate our destorted data. 

p1 = (2*pi*1+1j*2*pi*40);
s = tf('s');
Hc = 1/(s/p1+1);
Hd = c2d(Hc, 1/Fs, 'tustin');
[a, b] = tfdata(Hd);

yQ = real(filter(a{:}, b{:}, yQ));
yI = real(filter(a{:}, b{:}, yI));
%-------------------------------------------------------------------------
% Cut-off frequency of the lowpass filter must be much less than twice the 
% carrier frequency as high enough to capture the bandwidth of the received
% signal. We choose 110 Hz

Q = 2 .* lowpass(yQ, 110, Fs);
I = 2 .* lowpass(yI, 110, Fs);

%Reconstructing complex baseband signal and applying the second RRC
samplesRaw = I + 1i * Q;
samplesRec = upfirdn(samplesRaw, RRC, 1, 1);
samplesRec = samplesRec(length(RRC)/2+0.5:end-length(RRC)/2+0.5);

%% Phase Correction (Samples)
in_ph = real(samplesRec);
quad = imag(samplesRec);

figure('Name','Receiver Samples Scatter Plot','NumberTitle','off');
scatter(in_ph,quad,'filled');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-2,2,-2,2]);
hold on

% estimating the best fitting line
samplesRec_c = phase_corr_ls(samplesRec);
in_ph_c = real(samplesRec_c);
quad_c = imag(samplesRec_c);

scatter(in_ph_c,quad_c,'filled');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-2,2,-2,2]);
legend('W/O Phase Correction','With Phase Correction')
hold off

%% Symbol Detection
samp_packet_header = upsample(header, sps);
[r, lags] = xcorr(samplesRec_c, samp_packet_header);
r = r(length(samplesRec_c):length(r));
lags = lags(length(samplesRec_c):length(lags));
figure('Name','Correlation Status','NumberTitle','off');
stem(lags,r);
[M, I] =  max(abs(r));
start = lags(I); 

sign = 1;
if M < 0
    sign = -1;
end

% Extraction
endix = start + Ns * sps;
%endix = start + Ns;
if endix > N
    endix = N;
end

data = samplesRec_c(start:sps:endix);
%data = symbsRec_c(start:Ns);
deltasRec = zeros(1, length(data));  % Variable for comparison purposes only
bitsRec = zeros(1,length(data));

in_ph_data = real(data);
quad_data = imag(data);

figure('Name','Data Scatter Plot','NumberTitle','off');
scatter(in_ph_data,quad_data,'filled');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-2,2,-2,2]);
hold on

for i = 1:length(data)
    if in_ph_data(i) > 0
        bitsRec(i) = 1;
        deltasRec(i) = 1;
    else
        bitsRec(i) = 0;
        deltasRec(i) = -1;
    end
end
deltasRec = deltasRec * sign;
bitsRec = bitsRec * sign;

%plotting the received baseband
%deltasRec = upsample(deltasRec,sps);
%figure('Name','Receiver Baseband Representation','NumberTitle','off');
%plot(t,symbsRec);
%xlabel('Time(s)');
%ylabel('Symbols');
%hold on
%stem(t,deltasRec);

%% Error Calculation
errorBits = 0;
for i = 1:295
    if bitsRec(i) ~= Tx_bits(i)
        errorBits = errorBits + 1;
    end
end
display(errorBits)