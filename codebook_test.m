% This script evaluates beam pattern of different beamformer used in IA

clear;clc;
% --- codebook construction -----
Nt = 32;
M = 8;
type = 'sector';
[ BF0 ] = get_IA_BF( Nt, M, type );

% ------ spatial angle vector of ULA -------
angle_range = 90;
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:Nt-1).'*sin((kk - angle_range -1 )/180*pi));
end

% ------ beam pattern test -------
angle_test = 0;
xdata = linspace(-pi/2,pi/2,181);
figure
subplot(121)
polarplot(xdata,abs(FF'*BF0(:,1)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,3)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,5)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,7)),'linewidth',2);hold on
thetalim([-90 90])
ax = gca;
ax.ThetaTick = -90:22.5:90;
ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
grid on

subplot(122)
polarplot(xdata,abs(FF'*BF0(:,2)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,4)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,6)),'linewidth',2);hold on
polarplot(xdata,abs(FF'*BF0(:,8)),'linewidth',2);hold on
thetalim([-90 90])
ax = gca;
ax.ThetaTick = -90:22.5:90;
ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
grid on
grid on


