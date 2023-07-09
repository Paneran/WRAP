close all
clear all

%% transmit wave
key = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
bits = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,1];
%bits = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,0,1];

FS = 5000000;
FC = 1000000;
RS = 50000;
sps = FS/RS;

% RRC setup
span = 5;           % number of symbols to truncate RRC
rolloff = 0.5;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

transmit=readmatrix('data/stm_buffer_7.csv');
dim = size(transmit);
packet_length = dim(1)*dim(2);
transmit=reshape(transmit, packet_length, 1);
figure();
plot(transmit);
plot(transmit(1:end-1));
title("Transmit Wave: C");

%% Normalize Wave

norm_wave_M = (transmit-mean(transmit))/std(transmit) /25;

figure();
plot(norm_wave_M);
title("Norm Wave: MATLAB");


%% Costas Loop
df = 0;
f0 = FC+df; %1.05e6; % estimated frequency
f = [0, 0.1, 0.15, 1];
a = [1, 1, 0, 0];
order = 11;
lp = firpm(order, f, a);

phase = 0;
inph = zeros(1, order+1);
quad = zeros(1, order+1);
inph_ = zeros(1, order+1);
quad_ = zeros(1, order+1);
error = 0; 
% create lowpass as percent of this
%lp = [0.238053579626575	0.0644231488551473	0.0684583456734438	0.0684583456734438	0.0644231488551473	0.238053579626575];
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

figure();
plot(inph);
title("Costa Sample: MATLAB");

%% SRRC filtering
filtered_samps_M = conv(inph, RRC, 'same');
figure();
plot(filtered_samps_M);
title("filtered RRC: MATLAB");

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

figure();
sampled_M = sign(real_samp);
scatter(sampled_M, zeros(length(real_samp)));
title("Sampled: MATLAB");

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
