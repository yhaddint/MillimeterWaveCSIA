clear;clc;
f0 = 58e9;
c = 3e8;
Nrx = 36;
lambda0 = c/f0;
d0 = lambda0/2;
freq_range = linspace(58,64,10)*1e9;
% 
theta_range = linspace(-pi/10,pi/3,1e3);
AoA = 40/180*pi;
w_vec = exp(1j*pi*(0:Nrx-1)'*sin(AoA));
gain = zeros(length(freq_range), length(theta_range));

for freq_idx = 1:length(freq_range)
    freq = freq_range(freq_idx);
    lambda = c/freq;
    for theta_idx = 1:length(theta_range)
        arx = exp(1j*2*pi*d0/lambda*(0:Nrx-1)'*sin(theta_range(theta_idx)));
        gain(freq_idx,theta_idx) = arx'*w_vec;
    end
end
%
figure
for freq_idx = 1:length(freq_range)
    plot(theta_range/pi*180, (abs(gain(freq_idx,:))));hold on
end
xlabel('Angle [deg]')
ylabel('Gain [dB]')
grid on
%% theoretical analysis
figure
semilogy(theta_range/pi*180, (asin(sin(theta_range)*(1+0.002/2))-theta_range)/pi*180,'linewidth',2);
hold on
semilogy(theta_range/pi*180, (asin(sin(theta_range)*(1+0.015/2))-theta_range)/pi*180,'linewidth',2);
hold on
semilogy(theta_range/pi*180, (asin(sin(theta_range)*(1+0.05/2))-theta_range)/pi*180,'linewidth',2);
hold on 
semilogy(theta_range/pi*180, (asin(sin(theta_range)*(1+0.1/2))-theta_range)/pi*180,'linewidth',2);
hold on
grid on
set(gca,'FontSize',14)
xlabel('True Angle \theta [deg]')
ylabel('Squinted Angle asin((1+\epsilon)sin(\theta))-\theta [deg]')
legend('\epsilon = 0.002','\epsilon = 0.015','\epsilon = 0.05','\epsilon = 0.1')
xlim([0,60])
ylim([0.01,10])
yticks([0.01,0.1,1,10])

%%
for freq_idx = 1:length(freq_range)
    freq = freq_range(freq_idx);
    lambda = c/freq;
    
    arx_true = exp(1j*2*pi*d0/lambda*(0:Nrx-1)'*sin(AoA));
    
    % no squint assumption
    arx_dict_0th = exp(1j*pi*(0:Nrx-1)'*sin(AoA));
    gain_0th(freq_idx) = (arx_dict_0th)'*arx_true;
    
    
    % 1st order approximation
    arx_dict_1st = exp(1j*pi*(0:Nrx-1)'*sin(AoA)).*(1+...
                    1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA));
%     gain_1st(freq_idx) = (arx_dict_1st./norm(arx_dict_1st))'*arx_true;
    gain_1st(freq_idx) = (arx_dict_1st)'*arx_true;

    
    % 2nd order approximation
    arx_dict_2nd = exp(1j*pi*(0:Nrx-1)'*sin(AoA)).*(1+...
                    (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA))+...
                    (1/2) * (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^2);
    gain_2nd(freq_idx) = (arx_dict_2nd)'*arx_true;
    
    % 3rd order approximation
    arx_dict_3rd = exp(1j*pi*(0:Nrx-1)'*sin(AoA)).*(1+...
                    (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA))+...
                    (1/2)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^2+...
                    (1/2/3)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^3+...
                    (1/2/3/4)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^4+...
                    (1/2/3/4/5)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^5+...
                    (1/2/3/4/5/6)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^6+...
                    (1/2/3/4/5/6/7)* (1j*pi*(freq-f0)/f0*(0:Nrx-1).'*sin(AoA)).^7);
    gain_3rd(freq_idx) = (arx_dict_3rd)'*arx_true;
    
    % infinate terms
    arx_dict_inf = exp(1j*pi*(0:Nrx-1)'*sin(AoA)).*exp(1j*pi*(freq-f0)/f0*(0:Nrx-1)'*sin(AoA));
    gain_inf(freq_idx) = (arx_dict_inf)'*arx_true;
    
end

figure
plot(freq_range/1e9,20*log10(abs(gain_0th)),'linewidth',2);hold on
plot(freq_range/1e9,20*log10(abs(gain_1st)),'linewidth',2);hold on
plot(freq_range/1e9,20*log10(abs(gain_2nd)),'linewidth',2);hold on
plot(freq_range/1e9,20*log10(abs(gain_3rd)),'linewidth',2);hold on
plot(freq_range/1e9,20*log10(abs(gain_inf)),'linewidth',2);hold on

grid on
xlabel('Freq [GHz]')
ylabel('Matching Score')
legend('Squint nonaware','1st Approx.','2st Approx.','3rd Approx.','Inf Terms')
ymax = max(20*log10(abs(gain_1st)));
% ylim([ymax-60,ymax+10])
