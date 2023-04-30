close all
clear all

%% transmit wave
key = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
bits = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,1];
%bits = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,0,1];

FS = 4000000;
FC = 1000000;
RS = 50000;
sps = FS/RS;

% RRC setup
span = 5;           % number of symbols to truncate RRC
rolloff = 1;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

symbs = (bits-0.5)*2;
packet = symbs;
N = length(packet);

packet_length = N*sps;

deltas = upsample(packet,sps);
symbol_wave = conv(deltas, RRC);
symbol_wave = symbol_wave(length(RRC)/2+0.5:end-length(RRC)/2+0.5);
t = linspace(0, packet_length-1, packet_length);
wave = symbol_wave .* cos(2*pi*FC/FS*t);

wave_out = floor((wave+1)*2^11);
transmit=readmatrix('MyFile.csv');

figure();
subplot(2, 1, 1);
plot(wave_out);
title("Transmit Wave: MATLAB");
subplot(2, 1, 2);
plot(transmit(1:end-1));
title("Transmit Wave: C");

%% Normalize Wave

norm_wave_M = (wave_out-mean(wave_out))/std(wave_out)^2 *3;
norm_wave_C=readmatrix('norm_wave.csv');

figure();
subplot(2, 1, 1);
plot(norm_wave_M);
title("Norm Wave: MATLAB");
subplot(2, 1, 2);
plot(norm_wave_C);
title("Norm Wave: C");


%% Costas Loop
df = 0;
f0 = FC+df; %1.05e6; % estimated frequency
order = 5;

phase = 0;
inph = zeros(1, order+1);
quad = zeros(1, order+1);
inph_ = zeros(1, order+1);
quad_ = zeros(1, order+1);
error = 0; 
% create lowpass as percent of this
lp = [0.238053579626575	0.0644231488551473	0.0684583456734438	0.0684583456734438	0.0644231488551473	0.238053579626575];
integrator = 0;

kp = 8.5;
ki = 0.1;

for i = order+1:packet_length+order
    k = i - order;
    c = 2*cos(2*pi*f0/FS*(k-1)+phase);
    s =-2*sin(2*pi*f0/FS*(k-1)+phase);
    % demodulate the signal by mixing with NCO signals
    inph_(order+1) = norm_wave_M(k) * c;
    quad_(order+1) = norm_wave_M(k) * s;
    % low pass filter to get rid of high frequency components
    I = filter(lp, 1, inph_);
    Q = filter(lp, 1, quad_);
    inph(k) = I(order+1);
    quad(k) = Q(order+1);
    
    inph_(1:end-1) = inph_(2:end);
    quad_(1:end-1) = quad_(2:end);

    % mix in_ph and quad to get error signal
    error = I(order+1) * Q(order+1);

    % NCO
    integrator = integrator + ki * error;
    phase = phase+kp*error+integrator;
end

figure();
plot(quad);
hold on;
plot(inph);
hold off;

costa_loop_C=readmatrix('costa_samples.csv');
figure();
subplot(2, 1, 1);
plot(inph);
title("Costa Sample: MATLAB");
subplot(2, 1, 2);
plot(costa_loop_C);
title("Costa Sample: C");

%% SRRC filtering
filtered_samps_M = conv(inph, RRC, 'same');
srrc_samples=readmatrix('SRRC.csv');

figure();
subplot(2, 1, 1);
plot(filtered_samps_M);
title("filtered RRC: MATLAB");
subplot(2, 1, 2);
plot(srrc_samples);
title("Costa Sample: C");

%% Timing Recovery
i = filtered_samps_M;
q = quad;
%create zero crossing signal
zero_cross = ~(sign(i(1:end-1).*i(2:end))+1);

%PLL clock recovery pll
Kp_PLL = 0.01;
Ki_PLL = 0.05;
prev_acc = -0.00001;
real_samp = [];
quad_samp=[];
err = [];
integrator = sps;
for k = 1:packet_length-1
    phase_acc = prev_acc + 2*pi/sps;
    %Check if we have crossed zero on the phase accumulator
    %If the previous value was around -pi and the current value is around pi - we crossed it.
    p1 = wrapToPi(phase_acc);
    p2 = wrapToPi(prev_acc);
    if (p1< -pi * 0.75 && p2 > pi*0.75)
        %we have crossed zero - sample
        real_samp(end+1) = i(k);
        quad_samp(end+1) = q(k);
    end
    %Detect zero crossing and use it to tune PLL
    if (zero_cross(k))
        %we have a zero crossing - adjust freq
        err(end+1) = wrapToPi(phase_acc);
        integrator = integrator + err(end)*Ki_PLL;
        sps = err(end) * Kp_PLL + integrator;
    end
    prev_acc = phase_acc;
end

sampled=readmatrix('sampled.csv');
sampled_M = sign(real_samp);
figure();
subplot(2, 1, 1);
scatter(sampled_M, zeros(length(real_samp)));
title("Sampled: MATLAB");
subplot(2, 1, 2);
scatter(sampled, zeros(length(sampled)));
title("Sampled: C");

%{
figure();
plot (real_samp);
hold on;
plot(quad_samp);
hold off;
selected = real_samp +1j*quad_samp;
legend("Real sampled", "Imaginary sampled");
scatterplot(real_samp+quad_samp*1j);
figure;
plot(err);
legend("PLL error");
%}
