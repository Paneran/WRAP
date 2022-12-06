clear;
% This version of the code modulates and demodulates the signal.
% The packet header is found. No phase correction or finding of 
% the frequency is done.

%%
Fs = 1000;              % Sampling rate
Tmax = 5;               % Max time
packet_time = 10;       % length of packet (s)
N = packet_time * Fs;   % Number of samples

% Use these guys for plotting
t = linspace(0, packet_time - 1/Fs, N); % Time variable
f = linspace(-Fs/2,Fs/2, N);        % Frequency variable

Rs = 50;                        % Symbol rate in symbols/sec (baud)
Ns = Rs * Tmax;                 % Number of symbols to send
sps = Fs / Rs;                  % Number of samples per symbol
sps = N / (Rs * packet_time);
len_packet = packet_time * Rs;  % Number of possible symbols in packet

fc = 100; % Carrier frequency in Hz

% RRC setup
span = 1000;        % number of symbols to truncate RRC
len = sps * span;   % length of RRC
rolloff = 0.5;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');
rect = sinc(-sps*span/2:sps*span/2);


%% Symbol Generation
% Generate random set of bits
bits = randi([0, 1], 1, Ns);

% Map bits to BPSK symbols
symbs = 2 * (bits - 0.5);

% create packet key
key = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
packet_header = key;
for i=1:2
    packet_header = cat(2, key, packet_header);
end
h = length(packet_header);

% add 3 keys to the front of the data. the packet header. 
symbs = cat(2, packet_header, symbs);

% upsample
packet = cat(2, zeros(1, len_packet - Ns - h), symbs);
deltas = upsample(packet,sps);

%% Pulse Shaping (RRC)
symbol_wave = upfirdn(deltas, RRC, 1, 1);
symbol_wave = symbol_wave(length(RRC)/2+0.5:end-length(RRC)/2+0.5);
%% Carrier Modulation
comp_carr = exp(1i*2*pi*fc*t);
analytic_sig = symbol_wave .* comp_carr;
Tx = real(analytic_sig);

plot(t, Tx);

%% Noise
std = 0.1;
noise = std * randn(1,N);
Rx = Tx + noise;

%% Carrier Demodulation
% threshold varies per data. Taking an FFT in general is not a good idea in
% the future. Consider using a Phase Locked-Loop to find frequency. 
adj_fft = fftshift(fft(Rx.^2));
[peaks,loc] = findpeaks(abs(adj_fft),'MinPeakHeight',20);
fcr = floor(f(loc(end))/2);

phase = unifrnd(-pi,pi);

yQ = Rx .* -sin(2 * pi * fcr * t + phase);
yI = Rx .* cos(2 * pi * fcr * t + phase);

p1 = (2*pi*1+1j*2*pi*40);
s = tf('s');
Hc = 1/(s/p1+1);
Hd = c2d(Hc, 1/Fs, 'tustin');
[a, b] = tfdata(Hd);

yQ = real(filter(a{:}, b{:}, yQ));
yI = real(filter(a{:}, b{:}, yI));

% Cut-off frequency of the lowpass filter must be much less than twice the 
% carrier frequency as high enough to capture the bandwidth of the received
% signal. We choose 110 Hz2

Q = 2 * lowpass(yQ, 110, Fs);
I = 2 * lowpass(yI, 110, Fs);

%Reconstructing complex baseband signal and applying the second RRC
symbsRaw = I + 1i * Q;
symbsRec = upfirdn(symbsRaw, RRC, 1, 1);
symbsRec = symbsRec(length(RRC)/2+0.5:end-length(RRC)/2+0.5);

%% Phase correction

in_ph = real(symbsRec);
quad = imag(symbsRec);

figure('Name','Receiver Samples Scatter Plot','NumberTitle','off');
scatter(in_ph,quad,'filled');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-2,2,-2,2]);
hold on

c = [0.1, 0.00086024, 0.0017, 0.00086024, 1.9785, -0.9785];
symbsRec = phase_corr_ls(symbsRec);
in_ph = real(symbsRec);
quad = imag(symbsRec);

scatter(in_ph,quad,'filled');
xlabel('In-Phase');
ylabel('Quadrature');
axis([-2,2,-2,2]);
legend('W/O Phase Correction','With Phase Correction')
hold off
%% Clock Recovery
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

%% PLL

%% Symbol Detection
% Find beginning of packet. N^2 computations. 
samp_packet_header = upsample(packet_header, sps);
[r, lags] = xcorr(symbsRec, samp_packet_header);
r = r(length(RRC)/2-0.5:length(RRC)-2);
lags = lags(length(RRC)/2-0.5:length(RRC)-2);
[M, I] =  max(abs(r));
start = lags(I);
figure("Name", "Crosscorrelation with Packet Header")
stem(lags, r);
% Extract symbols
l = start+(Ns+h)*sps;
if l > N
    l = N;
end 
data = symbsRec(start:sps:l); % do +1 from start or -1 from end. 
sign = 1;
deltasRec = zeros(1, length(data));  % Variable for comparison purposes only

in_ph = real(data);
quad = imag(data);

for i = 1:length(data)
    if in_ph(i) > 0
        deltasRec(i) = 1;
    else
        deltasRec(i) = -1;
    end
end
header = deltasRec(1:h);
packet_d = deltasRec(h+1:end);
if header * packet_header' < 0
    sign = -1;
end

packet_d = packet_d * sign;
bitsRec = (packet_d+1)/2;

%% Error Calculation
errorBits = 0;
for i = 1:Ns
    if bitsRec(i) ~= bits(i)
        errorBits = errorBits + 1;
    end
end
display(errorBits)
