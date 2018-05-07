% This script is used to verify the main theory in [Sohabi et al 2016]
% It claims, any complex vector v can be exactly expressed by alpha(v1+v2)
% where v1 and v2 have all elements with unit magnitude

clear;clc;
N = 64;
v = randn(N,1) + 1j*randn(N,1);
v_abs = abs(v);
v_angle = phase(v);
vmax  = max(v_abs);
angle1 = v_angle - acos(v_abs./(2*vmax));
angle2 = v_angle + acos(v_abs./(2*vmax));
error = v - vmax*(exp(1j*angle1)+exp(1j*angle2));
figure;
plot(abs(error))