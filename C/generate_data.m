Fs = 5000000; % Sampling rate
Tmax = 0.1; % Max time
% Use these guys for plotting

Rs = 50000; % Symbol rate in symbols/sec (baud)
%Ns = Rs * Tmax; % Number of symbols to send
Ns = 15 + (8*1);
sps = Fs/Rs; % Number of samples per symbol
N = Ns*sps; % Number of samples


fc = 1000000; % Carrier frequency in Hz

%packet =[0 1 0 1 0 1 1 1 0 1 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 1 0 1 0 0 0 0];
packet = [0 1 1 0 1 0 0 0];
packet_symbs = 2 * (packet - 0.5);
packet_header = [1 1 1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1 1];
packet_w_header = [packet_header packet_symbs];
symbs = packet_w_header;
bits = (symbs + 1)/2;

deltas = upsample(symbs,sps);

window_rrc = rcosdesign(0.5,4,sps,'sqrt');
windowed_deltas = conv(deltas, window_rrc, 'same');
matlab_transmit = windowed_deltas .* cos(2 * pi * fc/Fs * 0:(N-1));
matlab_transmit = round(2048 + 2048 * matlab_transmit / max(abs(matlab_transmit)));