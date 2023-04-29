close all
clear all

Array=csvread('costa_samples.csv');
figure();

plot(Array(1000:1100))

bits = [0,1,0,0,1,0,0,1,0,0,...
        1,0,0,0,0,0,0,1,1,0,...
        0,1,0,0,0,1,1,0,1,1,...
        1,1,0,1,1,0,1,1,1,0,...
        0,0,1,0,0,1,1,1,0,1];

sps = 20;

% RRC setup

%{
span = 5;           % number of symbols to truncate RRC
rolloff = 0;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

bits = (bits-0.5)*2;
packet = bits;
deltas = upsample(packet,sps);
symbol_wave = conv(deltas, RRC);
symbol_wave = symbol_wave(length(RRC)/2+0.5:end-length(RRC)/2+0.5);
figure();
plot(symbol_wave);

%}