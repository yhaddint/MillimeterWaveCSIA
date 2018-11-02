% This script evaluates beam pattern of different beamformer used in
% initial access

clear;clc;
font_size = 14;

%% plot LS codebook
N_element = 256;
M_burst = 16; % T/Rx partial, so it's smaller than 64
W_mat = get_IA_BF(N_element, M_burst, 'sector_LS'); % Tx beamformer in IA stage

angle_range = 90;
FF = zeros(N_element, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:N_element-1).'*sin((kk - angle_range -1 )/180*pi));
end
for mm=1:M_burst
    pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
    % ------ beam pattern test -------
    angle_test = 0;
    % xdata = linspace(-pi/2,pi/2,181);
    xdata = linspace(-pi/2,pi/2,181);
    %
    figure(100)
    subplot(131)
    polarplot(xdata+pi,pattern,'linewidth',2);hold on
    set(gca,'FontSize', font_size)
    title('LS-Sec. BF [16]')
    rlim([-20,20])
    % thetalim([-90 90])
    thetalim([90,270])
    ax = gca;
    % ax.ThetaTick = -90:22.5:90;
    ax.ThetaTick = 90:22.5:270;
    ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
    grid on
end

%% plot FSM-KW codebook
N_element = 256;
M_burst = 16; % T/Rx partial, so it's smaller than 64
W_mat = get_IA_BF(N_element, M_burst, 'sector_FSM_KW'); % Tx beamformer in IA stage

angle_range = 90;
FF = zeros(N_element, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:N_element-1).'*sin((kk - angle_range -1 )/180*pi));
end
for mm=1:M_burst
    pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
    % ------ beam pattern test -------
    angle_test = 0;
    % xdata = linspace(-pi/2,pi/2,181);
    xdata = linspace(-pi/2,pi/2,181);
    %
    figure(100)
    subplot(132)
    polarplot(xdata+pi,pattern,'linewidth',2);hold on
    set(gca,'FontSize', font_size)
    title('FSM-Sec. BF [33]')
    rlim([-20,20])
    % thetalim([-90 90])
    thetalim([90,270])
    ax = gca;
    % ax.ThetaTick = -90:22.5:90;
    ax.ThetaTick = 90:22.5:270;
    ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
    grid on
end



%% plot PN codebook
N_element = 256;
M_burst = 16; % T/Rx partial, so it's smaller than 64
W_mat = get_IA_BF(N_element, M_burst, 'PN'); % Tx beamformer in IA stage

angle_range = 90;
FF = zeros(N_element, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:N_element-1).'*sin((kk - angle_range -1 )/180*pi));
end
for mm=1:M_burst
    pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
    % ------ beam pattern test -------
    angle_test = 0;
    % xdata = linspace(-pi/2,pi/2,181);
    xdata = linspace(-pi/2,pi/2,181);
    %
    figure(100)
    subplot(133)
    polarplot(xdata+pi,pattern,'linewidth',2);hold on
    set(gca,'FontSize', font_size)
    title('PN BF (prop.)')
    rlim([-20,20])
    % thetalim([-90 90])
    thetalim([90,270])
    ax = gca;
    % ax.ThetaTick = -90:22.5:90;
    ax.ThetaTick = 90:22.5:270;
    ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
    grid on
end



