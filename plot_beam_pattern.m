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

%% plot Hadamard codebook
N_element = 16;
M_burst = N_element; % T/Rx partial, so it's smaller than 64
temp = hadamard(N_element); % Tx beamformer in IA stage
W_mat = temp./norm(temp,'fro')*sqrt(N_element);

angle_range = 90;
FF = zeros(N_element, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:N_element-1).'*sin((kk - angle_range -1 )/180*pi));
end
for mm=1:16
    pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
    % ------ beam pattern test -------
    angle_test = 0;
    % xdata = linspace(-pi/2,pi/2,181);
    xdata = linspace(-pi/2,pi/2,181);
    %
    figure(101)
    subplot(1,4,floor((mm-1)/4)+1)
    polarplot(xdata+pi,pattern,'linewidth',2);hold on
    set(gca,'FontSize', font_size)
    title('Hadamard BF')
    rlim([-20,20])
    % thetalim([-90 90])
    thetalim([90,270])
    ax = gca;
    % ax.ThetaTick = -90:22.5:90;
    ax.ThetaTick = 90:22.5:270;
    ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
    grid on
end

%% plot FSM-KW codebook in 3D
clear;clc;
plot3Dpolar = 1; 
N_element = 256;
N_az = 32;
N_el = 8;

M_burst = 16; % T/Rx partial, so it's smaller than 64
M_burst_az = 8;
M_burst_el = 2;

az_lim = pi/3;
el_lim = pi/6;

W_mat = get_IA_BF_3D(N_az, N_el,...
                     M_burst_az, M_burst_el,...
                     'sector_FSM_KW', az_lim, el_lim); % Tx beamformer in IA stage

angle_range = 90;
angle_num = (angle_range*2+1);
FF_az = zeros(N_az, (angle_range*2+1));
FF_el = zeros(N_el, (angle_range*2+1));

for kk = 1:angle_num
    FF_az(:,kk) = exp(1j*pi*(0:N_az-1).'*sin((kk - angle_range -1 )/180*pi));
    FF_el(:,kk) = exp(1j*pi*(0:N_el-1).'*sin((kk - angle_range -1 )/180*pi));
end

for kk=1:angle_num^2
    kk_el = floor((kk-1)/angle_num)+1;
    kk_az = kk - (kk_el-1)*angle_num;
    FF(:,kk) = kron( FF_el(:,kk_el), FF_az(:,kk_az) );
end

for mm=1:16
    pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
    % ------ beam pattern test -------
    pattern_3D = reshape(pattern, angle_num, angle_num);
    [XDATA, YDATA] = meshgrid(-angle_range:angle_range,-angle_range:angle_range);
    
    figure(99)
    subplot(M_burst_el, M_burst_az, mm)
    mesh(XDATA, YDATA, pattern_3D);hold on
    xlabel('El. [deg]')
    ylabel('Az. [deg]')
    view(180, 90);
    
    if plot3Dpolar
        figure((100+mm))
        subplot(M_burst_az, M_burst_el,mm)
        
        % need some weird shift to look correct
        theta_range  = 90+(-angle_range:angle_range);
        phi_range = -angle_range:angle_range;
        patternCustom(pattern_3D, theta_range, phi_range)
    end
end

