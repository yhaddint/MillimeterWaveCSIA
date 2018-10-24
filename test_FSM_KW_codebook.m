% This script evaluates beam pattern of different beamformer used in
% tracking

clear;clc;
% --- codebook construction -----
Nt = 16;
A_stopband = 10; % attenuation outside mainlobe (dB)
beam_width = 5/180*pi;
steer_dir = 0/180*pi;
vec1 = get_FSM_KW_codebook( steer_dir, beam_width, Nt, A_stopband);

%% ------- test hybrid approximation-------
% phase_bits = 10; % phase shifter quantization; set to 10 if don't care
% vec = get_hybrid_approx( vec1, phase_bits );
% vec_norm = vec./norm(vec);
% 
% % ------ spatial angle vector of ULA -------
% angle_range = 90;
% for kk = 1:(angle_range*2+1)
%     FF(:,kk) = exp(1j*pi*(0:Nt-1).'*sin((kk - angle_range -1 )/180*pi));
% end
% response = 20*log10(abs(FF'*vec1)+1e-1);
% 
% % ------ beam pattern test -------
% angle_test = 0;
% % xdata = linspace(-pi/2,pi/2,181);
% xdata = linspace(-pi/2,pi/2,181);
% figure
% plot(xdata/pi*180,response)
% grid on
% % ----- beam peattern plot
% figure
% polarplot(xdata+pi,response,'linewidth',2);hold on
% rlim([-10,10])
% % thetalim([-90 90])
% thetalim([90,270])
% ax = gca;
% % ax.ThetaTick = -90:22.5:90;
% ax.ThetaTick = 90:22.5:270;
% ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
% grid on
%% test entire codebook
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
    figure(99)
    plot(xdata/pi*180,pattern);hold on
    grid on
    %
    figure(100)
    polarplot(xdata+pi,pattern,'linewidth',2);hold on
    rlim([-20,20])
    % thetalim([-90 90])
    thetalim([90,270])
    ax = gca;
    % ax.ThetaTick = -90:22.5:90;
    ax.ThetaTick = 90:22.5:270;
    ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
    grid on
end



