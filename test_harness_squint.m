clear;clc;

fc_num = 16;
fc_range = linspace(60, 62, fc_num+1)*1e9;
f0 = 60e9;
c_speed = 3e8;
lambda0 = c_speed/f0;
dtau_az = 1/(2e9);
dtau_el = 0.01765*16/c_speed;

N_az = 16;
N_el = 16;
theta_az_num = 100;
theta_az_range = linspace(-pi/2,pi/2,theta_az_num);
theta_el_num = 100;
theta_el_range = linspace(-pi/2,pi/2,theta_az_num);

for ff=1:fc_num
    fc = fc_range(ff);
    lambda = c_speed/fc;
    dphi_az = 2*pi*fc*dtau_az;
    dphi_el = 2*pi*fc*dtau_el;
    
    array_win = kaiser(N_az,3)./norm(kaiser(N_az,3))*sqrt(N_az);
    v_az_vec(:,ff) = exp(1j*(0:N_az-1).'*dphi_az);%.*array_win;
    v_el_vec(:,ff) = exp(1j*(0:N_el-1).'*dphi_el);

    
    for tt=1:theta_az_num
        h_az(:,tt) = exp(1j*2*pi*lambda0/2/lambda*(0:N_az-1).'*...
                        sin(theta_az_range(tt)))/sqrt(N_az)^2;
        h_el(:,tt) = exp(1j*2*pi*lambda0/2/lambda*(0:N_el-1).'...
                        *sin(theta_el_range(tt)))/sqrt(N_el)^2;

        gain_az(ff,tt) = abs(h_az(:,tt)'*v_az_vec(:,ff))^2;
        gain_el(ff,tt) = abs(h_el(:,tt)'*v_el_vec(:,ff))^2;
        
    end
end
%%
fig = figure('units','inch','position',[0,0,16,4]);
% subplot(211)
% f_idx_temp = reshape((1:64).',16,4);
% f_idx = reshape(f_idx_temp(1:4,:),16,1);
f_idx = 1:16;
plot(theta_az_range/pi*180, 10*log10(gain_az(f_idx,:)),'linewidth',2);
set(gca,'FontSize',14)
xlabel('Theta Azimuth [deg]')
ylabel('Gain [dB]')
ylim([-30,0])
xlim([-90,90])
xticks(-90:30:90)
grid on
legendtext = {};
for ff=1:length(f_idx)
    legendtext{ff} = [num2str(fc_range(f_idx(ff))/1e9,'%.2f') 'GHz'];
end

% legendtext = {'a','b'};
legend(legendtext,'Location','eastoutside','FontSize',11)
% subplot(212)
% f_idx = 1:4;
% plot(theta_el_range/pi*180, 10*log10(gain_el(f_idx,:)),'linewidth',2);
% set(gca,'FontSize',14)
% xlabel('Theta Elevation [deg]')
% ylabel('Gain [dB]')
% ylim([-30,0])
% xlim([-90,90])
% xticks(-90:30:90)
% grid on
% legendtext = {};
% for ff=1:fc_num
%     legendtext{ff} = [num2str(fc_range(f_idx(ff))/1e9,'%.1f') 'GHz'];
% end
% 
% % legendtext = {'a','b'};
% legend(legendtext)
