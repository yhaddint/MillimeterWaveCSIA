clear;clc;
%
Ts = 1e-13; % time resolution in simulation is 0.1ps
N_rx = 16; % 16 receiver antenna
tone_num = 128; % simulation uses 128 tones with different phases
tone_range = linspace(59, 61, tone_num)*1e9; % Freq of tones in [Hz]
t = (0:2e6).'*Ts; % simulating 200us

% Generate 128 tones in the Tx
for tone_idx = 1:tone_num
    tone_freq = tone_range(tone_idx);
    tone_phase = exp(1j*2*pi*rand); % some random starting phase of each tone
    tone(:,tone_idx) = 1/sqrt(tone_num) * exp(1j*2*pi*tone_freq*t) * tone_phase;
end

% Tx signal with 1e-13 sample duration (RF simulating, i.e. carrier)
sig_tx = sum(tone,2); % transmit signal

%% This simulation supports 167 different angle candidates within [-90,90]
angle_cand = zeros(167,1);
for kk=1:167
    angle_cand(kk) = asin(2*3e8/0.005*(kk-84)*1e-13)/pi*180;
end

% figure
% plot(angle_cand,(1:167)-84,'o')
% xlabel('True AoA')
% ylabel('Over-the-Air Delay b/w Adjacent Ant [1e-13]')
% grid on

% This script simulating 60GHz and 
%   -83.0000  -84.8736
%   -82.0000  -79.7369
%   -81.0000  -76.4095
%    ...        ...
% It means we simulate AoA -84.8736 deg, the wave arrives at each antenna with
% time different 83*1e-13 second for a 60GHz critically spaced ULA
delay_air_LUT = [(1:167).'-84, angle_cand];

%% Simulating true time delay over air due to propagation

% Randomly create an AoA that is supported by this script
idx = randi([1,167]);
AoA_true = delay_air_LUT(idx,2);
delay_air = delay_air_LUT(idx,1);

% Initialize a zero matrix, each column for an antenna
% It should have plenty space to apply true time delay operation (shifting)
rx_sig = zeros(length(t) + 100 * N_rx, N_rx); 

% Apply true time delay in the received antennas
if delay_air>0 % right antenna sees more delay of signal prop. from the air
    for nn=1:N_rx
        delay_to_n_ant = delay_air*(nn-1);
        rx_sig_seg = (1:length(t)) + delay_to_n_ant;
        rx_sig(rx_sig_seg, nn) = sig_tx;
    end
elseif delay_air<0 % right antenna sees less delay of signal prop. from the air
    for nn=1:N_rx
        delay_to_n_ant = -delay_air*(nn-1);
        rx_sig_seg = (1:length(t)) + delay_to_n_ant;
        rx_sig(rx_sig_seg,(N_rx-nn+1)) = sig_tx;
    end
else
    for nn=1:N_rx
        rx_sig_seg = (1:length(t));
        rx_sig(rx_sig_seg, nn) = sig_tx;
    end
end

%% Simulating true time delay introduced by TTD circuitry 

% 0.5ns extra delay by TTD front-end
TTD_delay = round(0.5e-9/Ts);

% Initialize a zero matrix, each column for an antenna
% It should have plenty space to apply true time delay operation (shifting)
rx_TTD_sig = zeros(size(rx_sig,1) + TTD_delay * N_rx, N_rx); 

% Apply TTD in each antenna
for nn=1:N_rx
    TTD_to_n_ant = TTD_delay * (nn-1);
    rx_sig_TTD_seg = (1:size(rx_sig,1)) + TTD_to_n_ant;
    rx_TTD_sig(rx_sig_TTD_seg, nn) = rx_sig(:,nn);
end

% Combine signal from N_rx antenna
rx_TTD_sig_comb = sum(rx_TTD_sig,2);

%% Post processing (spectrogram)
for tone_idx = 1:tone_num
    pow(tone_idx) = abs(rx_TTD_sig_comb(1:length(t))'*tone(:,tone_idx));
end
%% Plot 'PSD'
figure
plot(tone_range/1e9, 20*log10(pow))
grid on
xlabel('Freq [GHz]')
ylabel('PSD [dB]')

%% Freq to angle map
[~, maxidx] = max(pow);
peak_freq = tone_range(maxidx);
AoA_est = -asin(2*(mod(peak_freq/2e9+0.5,1)-0.5))/pi*180;
fprintf('True AoA = %.2f deg, Est. AoA = %.2f deg\n', AoA_true, AoA_est)
