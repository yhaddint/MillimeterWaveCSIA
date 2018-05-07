% This script evaluates beam pattern of different beamformer used in IA

clear;clc;
% --- codebook construction -----
Nt = 64;
A_stopband = 30; % attenuation outside mainlobe (dB)
vec = get_FSM_KW_codebook( 45/180*pi, 0/180*pi, Nt, A_stopband);
% vec = exp(1j*pi*(0:Nt-1).'*sin(45/180*pi));
vec_norm = vec./norm(vec);

% ------- test hybrid approximation-------
phase_bits = 3;
[ v_approx ] = get_hybrid_approx( vec_norm, phase_bits );

% ------ spatial angle vector of ULA -------
angle_range = 90;
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:Nt-1).'*sin((kk - angle_range -1 )/180*pi));
end
response = 20*log10(abs(FF'*v_approx)+1e-1);

% ------ beam pattern test -------
angle_test = 0;
% xdata = linspace(-pi/2,pi/2,181);
xdata = linspace(-pi/2,pi/2,181);
figure
plot(xdata/pi*180,response)
grid on
%%
figure
polarplot(xdata,response,'linewidth',2);hold on
rlim([-20,20])
thetalim([-90 90])
grid on



