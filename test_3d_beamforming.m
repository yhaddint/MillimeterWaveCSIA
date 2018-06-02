% This script evaluates beam pattern of different beamformer used in IA

clear;clc;
% --- codebook construction -----
Nt_x = 16;
Nt_y = 16;
A_stopband = 30; % attenuation outside mainlobe (dB)
% vec_x = get_FSM_KW_codebook( 15/180*pi, 0/180*pi, Nt_x, A_stopband);
% vec_y = get_FSM_KW_codebook( 15/180*pi, 0/180*pi, Nt_y, A_stopband);
vec_x = exp(1j*pi*(0:Nt_x-1).'*sin(45/180*pi));
vec_y = exp(1j*pi*(0:Nt_y-1).'*sin(45/180*pi));
vec = kron(vec_x,vec_y);
vec_norm = vec./norm(vec);

% ------- test hybrid approximation-------
phase_bits = 10;
[ v_approx ] = get_hybrid_approx( vec_norm, phase_bits );

% ------ spatial angle vector of ULA -------
angle_range = 90;
for k1 = 1:(angle_range*2+1)
    for k2 = 1:(angle_range*2+1)
        ang_sig_x = exp(1j*pi*(0:Nt_x-1).'*sin((k1 - angle_range -1 )/180*pi));
        ang_sig_y = exp(1j*pi*(0:Nt_y-1).'*sin((k2 - angle_range -1 )/180*pi));
        response(k1,k2) = 20*log10(abs(kron(ang_sig_x,ang_sig_y)'*v_approx+1e-1));
    end
end
% response = 20*log10(abs(FF'*v_approx)+1e-1);

%% ------ beam pattern test -------
angle_test = 0;
% xdata = linspace(-pi/2,pi/2,181);
xdata = linspace(-pi/2,pi/2,181).';
ydata = linspace(-pi/2,pi/2,181).';
[X,Y] = meshgrid(xdata,ydata);
figure
s = surf(X/pi*180,Y/pi*180,response);
s.EdgeColor  = 'none'
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
zlabel('Gain (dB)')
colorbar
colormap('jet')
caxis([-20,15])
zlim([-20,20])
%%
patternCustom(response,xdata,ydata)