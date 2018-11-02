% Evaluation of kappa, parameter that indicate SNR drop due to STO and CFO

clear;clc;
P = 128; % total number of subcarrier
K_range = [0,8,16,32,64]; % overlap of beamformer window between BS and UE
TS_OFDM = 1/240e3; % OFDM duration in [sec]
Ts = TS_OFDM/P; % sample duration in [sec]
fc = 28; % carrier in [GHz]
CFO_ppm_range = 10.^(linspace(-2,1,100)); % CFO range in [ppm]
CFO_range = CFO_ppm_range*fc*1e3*Ts*2*pi; % CFO range in [rad/sample]

figure
for kk=1:length(K_range)
    K = K_range(kk);
    for ss = 1:length(CFO_range)

        CFO = CFO_range(ss);
        result(ss) = (2-real(exp(1j*CFO*K))-real(exp(1j*CFO*(P-K))))/(P^2-P^2*real(exp(1j*CFO)));
    end
    semilogx(CFO_ppm_range,10*log10(result),'linewidth',2)
    hold on
end
grid on
xlabel('Absolute CFO [ppm]')
ylabel('Eff. SNR Drop \zeta(\epsilon_F,\epsilon_T) [dB]')
legend('K=0','K=5','K=10','K=30','K=60')
ylim([-10,0])
xticks([0.01,0.1,1,10])
xticklabels({'0.01','0.1','1','10'})

%%
SNR_dB = -40:0.01:0;
Nc = 1;
M = 64;
P = 128;
L = 1;
kappa = 1;
NB = 128*4;
PFA = 0.05;
SNR = 10.^(SNR_dB/10);
PMD_PT = qfunc((SNR - sqrt(Nc/M)/P*qfuncinv(PFA))./sqrt(11*SNR.^2/M/L+Nc/P^2/M));
% PMD_NT = qfunc((kappa*SNR + sqrt(Nc/M)/P*qfuncinv(1/NB)-0.78*sqrt(Nc/M)/P/qfuncinv(1/NB)*log(-log(PFA)))...
%                 ./sqrt(11*SNR.^2*kappa^2/M/L+Nc/P^2/M));
PMD_NT = qfunc((kappa*SNR - sqrt(Nc/M)/P*qfuncinv(1/NB)+0.78*sqrt(Nc/M)/P/qfuncinv(1/NB)*log(-log(PFA)))...
                ./sqrt(11*SNR.^2*kappa^2/M/L+Nc/P^2/M));

figure
semilogy(SNR_dB,PMD_PT,'-','linewidth',2);hold on
% semilogy(SNR_dB,PMD_NT,'--','linewidth',2)
grid on
xlabel('SNR [dB]')
ylabel('PMD')
legend('PT','NT')
%%
Pfa = 0.05;
STO_max_range = 20:1000;
for ss=1:length(STO_max_range)
    STO_max = STO_max_range(ss);
    TH_factor(ss) = qfuncinv(1/STO_max)-0.78*log(-log(1-Pfa))/qfuncinv(1/STO_max);
end
% figure
plot(STO_max_range,TH_factor,'linewidth',2);hold on
plot(STO_max_range,qfuncinv(Pfa)*ones(length(STO_max_range),1),'linewidth',2);hold on
grid on
xlabel('STO MAX')
ylabel('TH factor')
ylim([0,max(TH_factor)*1.5]);
