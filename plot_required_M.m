clear;clc
%
blue = [0    0.4470    0.7410];
red = [   0.8500    0.3250    0.0980];
yellow = [0.9290    0.6940    0.1250];
purple = [ 0.4940    0.1840    0.5560];

open('Nt32Nr16_alignment_vs_SNR_and_bursts.fig');
h = gcf
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); 

for ii=1:12
SNR(:,ii) = (dataObjs{2}(ii).XData).';
Alignment(:,ii) = (dataObjs{2}(ii).YData).';
end
%%
DoubleAlign = (Alignment(:,1:6).*Alignment(:,7:12))';
%
% figure
% plot(SNR(1:9,1)',DoubleAlign(:,1:9))
%%
M_and_align0 = [[16,20,24,34,48,68]',DoubleAlign(:,1:9)];

load Alignment_vs_M_L1
M_and_align1 = [M_range',SNR_vs_Align_L1'];
load Alignment_vs_M_L1_new
M_and_align2 = temp;
new_data = sort([M_and_align0;M_and_align1;M_and_align2],1);

%%
% figure
% plot(SNR(1:9,1)',new_data(:,2:end),'linewidth',2)
% legendtext = [];
% for M_index = 1:length(new_data(:,1))
%     legendtext = [legendtext;'M=',num2str(new_data(M_index,1))];
% end
% legend(legendtext)
%%
for ss=1:length(SNR_range)
    minindex = min(find(new_data(:,ss+1)>0.95));
    if minindex > 0
        Required_M(ss) = new_data(minindex,1)
    else
        Required_M(ss) = 0;
    end
end
figure(99)
p1 = semilogy([SNR_range,12.5,15,17.5,20], [Required_M,ones(1,4)*Required_M(end)],...
            'o',...
            'linewidth',2,...
            'markersize',8,...
            'color',blue);
hold on
grid on
% ylim([10,200])
% xlim([-10,20])
% xlabel('SNR [dB]')
% ylabel('Required M')
%%
load Alignment_vs_M_L2.mat
for ss=1:length(SNR_range)
    minindex = min(find(Align_vs_SNR_L2(:,ss+1)>0.95));
    if minindex > 0
        Required_M_L2(ss) = Align_vs_SNR_L2(minindex,1);
    else
        Required_M_L2(ss) = 0;
    end
end
figure(99)
p2 = semilogy([SNR_range,12.5,15,17.5,20], [Required_M_L2,ones(1,4)*Required_M_L2(end)],...
            's',...
            'linewidth',2,...
            'markersize',8,...
            'color',red);
hold on
grid on
% ylim([10,200])
% xlim([-10,20])
% xlabel('SNR [dB]')
% ylabel('Required M')
%%
load Alignment_vs_M_L3.mat
for ss=1:length(SNR_range)
    minindex = min(find(Align_vs_SNR_L3(:,ss+1)>0.95));
    if minindex > 0
        Required_M_L3(ss) = Align_vs_SNR_L3(minindex,1);
    else
        Required_M_L3(ss) = 0;
    end
end
figure(99)
p3 = semilogy([SNR_range,12.5,15,17.5,20], [Required_M_L3,ones(1,4)*Required_M_L3(end)],...
            'x',...
            'linewidth',2,...
            'markersize',8,...
            'color',yellow);
hold on
grid on
% ylim([10,256])
% xlim([-10,20])
% xlabel('SNR [dB]')
% ylabel('Required Measurement')
%%
load Alignment_vs_M_L4.mat
for ss=1:length(SNR_range)
    minindex = min(find(Align_vs_SNR_L4(:,ss+1)>0.95));
    if minindex > 0
        Required_M_L4(ss) = Align_vs_SNR_L4(minindex,1);
    else
        Required_M_L4(ss) = 0;
    end
end
figure(99)
p4 = semilogy([SNR_range,12.5,15,17.5,20], [Required_M_L4,ones(1,4)*Required_M_L4(end)],...
            '*',...
            'linewidth',2,...
            'markersize',8,...
            'color',purple);
hold on
grid on
ylim([10,256])
xlim([-10,20])
xlabel('SNR [dB]')
ylabel('Required Measurement')
%% Curve fitting
figure(99)
xdata = -10:0.1:2.5;
ydata = 2.^(-xdata/6+5.2);
semilogy(xdata,ydata,'--','linewidth',1,'color',blue);hold on
ydata = 2.^(-xdata/6+5.5);
semilogy(xdata,ydata,'--','linewidth',1,'color',red);hold on
ydata = 2.^(-xdata/6+5.8);
semilogy(xdata,ydata,'--','linewidth',1,'color',yellow);hold on
ydata = 2.^(-xdata/6+6);
semilogy(xdata,ydata,'--','linewidth',1,'color',purple);hold on

xdata = 7.5:0.1:20;
ydata = ones(length(xdata),1)*Required_M(end);
semilogy(xdata,ydata,'--','linewidth',1,'color',blue);hold on
ydata = ones(length(xdata),1)*Required_M_L2(end);
semilogy(xdata,ydata,'--','linewidth',1,'color',red);hold on
ydata = ones(length(xdata),1)*Required_M_L3(end);
semilogy(xdata,ydata,'--','linewidth',1,'color',yellow);hold on
ydata = ones(length(xdata),1)*Required_M_L4(end);
semilogy(xdata,ydata,'--','linewidth',1,'color',purple);hold on

legend([p1(1),p2(1),p3(1),p4(1)],'L=1 Paths','L=2 Paths','L=3 Paths','L=4 Paths');


y_wanted = 2.^(4:8);
yticks(y_wanted)
yticklabels({'16','32','64','128','256'})








