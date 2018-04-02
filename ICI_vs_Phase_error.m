clear;clc;
rng(3)
% system parameters
Nt = 16;
fft_size = 1024;
DFT = dftmtx(fft_size);
Ts = 1/(240e6);
fc = 28e9;
MCtimes = 5e1;
connected_LO = 1;
OFDM_blk_num = 5e1;

% % Connected LO CFO evaluation
% %-------------------------------------
% % CFO
% %-------------------------------------
% CFO_rms = 8e3;
% for MCindex = 1:MCtimes
%     if connected_LO
%         CFO = CFO_rms;
%     else
%         CFO = (rand(Nt,1)*2-1)*sqrt(3) * CFO_rms;
%     end
% 
%     PN_seq_save = ones(Nt,1);
%     for tt = 1:OFDM_blk_num
%         PE_seq_DA(:,1) = PN_seq_save;
%         for ll=1:fft_size-1
%             PE_seq_DA(:,ll+1) = PE_seq_DA(:,ll).*exp(1j*CFO*Ts);
%         end
%         PN_seq_save = PE_seq_DA(:,end);
%         for nn=1:Nt
%             gain(nn) = DFT(:,1)'*(PE_seq_DA(nn,:).'.*DFT(:,1))/fft_size;
%             ICI(nn) = DFT(:,1)'*(PE_seq_DA(nn,:).'.*DFT(:,2))/fft_size;
%         end
%         gain_time_evo(tt,MCindex) = abs(sum(gain));
%         ICI_time_evo(tt,MCindex) = abs(sum(ICI));
%     end
% end
% gain_time_evo_mean = 10*log10(mean(gain_time_evo,2));
% ICI_time_evo_mean = 10*log10(mean(ICI_time_evo,2));
% 
% figure
% plot(1:OFDM_blk_num,gain_time_evo_mean);hold on
% plot(1:OFDM_blk_num,ICI_time_evo_mean);hold on
% grid on
% xlabel('OFDM Block Number')
% ylabel('Gain of Signal and ICI (dB)')

%% Connected LO PN evaluation
%-------------------------------------
% Phase Noise Specification
% FOM = L(f0,df) + 20log10(df/f0)+10log10(P_VCO)
%-------------------------------------

VCO_FOM = -114 + 20*log10(1e6/28e9) + 10*log10(27)+20;
P_VCO = 27; % mW as unit
VCO_c = 10^(VCO_FOM/10)/P_VCO; % parameter c in PN specs
PN_sigma2 = VCO_c*4*pi^2*fc^2*Ts; % Time domain variance of PN Wiener process

% Connected LO
PN_sigma = sqrt(PN_sigma2); % Pick PN power and take square root
for MCindex = 1:MCtimes
    PN_seq_save = ones(Nt,1);
    for tt = 1:OFDM_blk_num
        PE_seq_DA(:,1) = PN_seq_save;
        for ll=1:fft_size-1
            if connected_LO
                PE_seq_DA(:,ll+1) = PE_seq_DA(:,ll).*exp(1j*ones(Nt,1)*randn*PN_sigma);
            else
                PE_seq_DA(:,ll+1) = PE_seq_DA(:,ll).*exp(1j*randn(Nt,1)*PN_sigma);
            end
        end
        PN_seq_save = PE_seq_DA(:,end);
        for nn=1:Nt
            gain(nn) = DFT(:,1)'*(PE_seq_DA(nn,:).'.*DFT(:,1))/fft_size;
            ICI(nn) = DFT(:,1)'*(PE_seq_DA(nn,:).'.*DFT(:,2))/fft_size;
        end
        gain_time_evo(tt,MCindex) = abs(sum(gain));
        ICI_time_evo(tt,MCindex) = abs(sum(ICI));
    end
end
gain_time_evo_mean = 10*log10(mean(gain_time_evo,2));
ICI_time_evo_mean = 10*log10(mean(ICI_time_evo,2));
%%
figure
plot(1:OFDM_blk_num,gain_time_evo_mean);hold on
plot(1:OFDM_blk_num,ICI_time_evo_mean);hold on
grid on
xlabel('OFDM Block Number')
ylabel('Gain(dB)')