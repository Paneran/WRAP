clear;
close all;
% This version of the code modulates and demodulates the signal.
% The packet header is found. No phase correction or finding of 
% the frequency is done.

%%
Fs = 5000000;              % Sampling rate
N = 50000;
Tmax = N/Fs;               % Max time
packet_time = Tmax*2;       % length of packet (s)

% Use these guys for plotting
t = linspace(0, packet_time - 1/Fs, 2*N); % Time variable
f = linspace(-Fs/2,Fs/2, N);        % Frequency variable

Rs = 1000000;                        % Symbol rate in symbols/sec (baud)
Ns = Rs * Tmax;                 % Number of symbols to send
sps = Fs / Rs;                  % Number of samples per symbol
len_packet = packet_time * Rs;  % Number of possible symbols in packet

fc = 1000000; % Carrier frequency in Hz

% RRC setup
span = 6;        % number of symbols to truncate RRC
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
phase = unifrnd(-pi,pi);
comp_carr = exp(1i*2*pi*fc*t+phase);
analytic_sig = symbol_wave .* comp_carr;
Tx = real(analytic_sig);

plot(t, Tx);

%% Noise
std = 0.1;
noise = std * randn(1,2*N);
y = Tx + noise;


% increment in time
t = 0:1/Fs:Tmax-1/Fs;
% use costas loop to demodulate
f0 = 1e6; %9.7706e5; % estimated frequency
ph = zeros(1, N+1);
inph = zeros(1, N);
quad = zeros(1, N);
err = zeros(1, N);
% nyquist rate: Fs/2=2500000
% create lowpass as percent of this
f = [0, 0.4, 0.41, 1];
a = [1, 1, 0, 0];
order = 5;
lp = firpm(order, f, a);

for i = order+1:N
    c = 2*cos(2*pi*f0*t(i)+ph(i));
    s = -2*sin(2*pi*f0*t(i)+ph(i));
    % demodulate the signal by mixing with NCO signals
    inph(i) = y(i) * c;
    quad(i)  = y(i) * s;
    % low pass filter to get rid of high frequency components
    inph(i) = mean(filter(lp, 1, inph(i-order:i)));
    quad(i) = mean(filter(lp, 1, quad(i-order:i)));
    
    % mix in_ph and quad to get error signal
    err(i) = inph(i)*quad(i);

    % NCO
    u = 1;
    ph(i+1) = ph(i)-u*err(i);
end
inph(end-9:end)=zeros(1, 10);
figure();
plot(1:N, err)
title("Error Plot");

figure();
scatter(inph, quad)
title("Constellation Diagram")

figure();
plot(t, inph)
title("Recovered Carrier")

figure();
plot(1:N+1, ph)
title("Phase Tracker")

span = 5;
srrc = rcosdesign(0.5, span, 6, 'sqrt');
l = inph+j*quad;
l = upfirdn(l, srrc,1,1);
l = l(length(srrc)/2+0.5:end-length(srrc)/2+0.5);

%% Clock recovery
i = real(l);
q = imag(l);
%create zero crossing signal
zero_cross = sign(i(1:end-1).*i(2:end));
zero_cross = zero_cross - mean(zero_cross);

%Get data rate
sps = 5;
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



