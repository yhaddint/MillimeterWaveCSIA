clear;clc;

Tb = 74.6e-6; % Constant
BWSS = 57.4e6; % Constant
BWtot = 250e6; % Constant
TRS = 1.25e-3;
NUE_range = 1:50;

SNR_range = -30:2.5:-10
% PMD_range = [0.99,0.9,0.65,0.25,0.05,0.01,0.01,0.01,0.01];
% PMA_range = [0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7];

PMD_range = [0.99,0.97,0.9,0.5,0.1,0.01,0.01,0.01,0.01];
PMA_range = [0.7,0.7,0.7,0.3,0.1,0.05,0.03,0.025,0.01,0.01];


% TSS_range = [10,15,20,40,80,160]*1e-3;
TSS = 20e-3
M_burst = 64; % Varying
latency = zeros(length(NUE_range),1);
for ii=1:length(SNR_range)

%     TSS = TSS_range(ii); % Varying 
    NUE = 50; % Varying
    PMD = PMD_range(ii);
    PMA = PMA_range(ii);
    
    NRS = max(floor((TSS-Tb*M_burst)/TRS),1); % adaptively adjusted
    
    % compute overhead
    OH(ii) = (Tb*BWSS*M_burst)/(TSS*BWtot);
    K=0:100;
    latency_IA = sum(PMD.^K*(1-PMD).*K*TSS);
    Kcyc = floor((NUE-1)/NRS);
    Kres = NUE-Kcyc*NRS;
    if Kcyc==0
        latency(ii) = 0;
        for qq=1:NUE
            latency(ii) = latency(ii)+(qq*TRS)/NUE;
        end
    else
        for kk=1:Kcyc
            for qq=1:NRS
                latency(ii) = latency(ii)+((kk-1)*TSS+qq*TRS)/NUE;
            end
        end
        for qq=1:Kres
            latency(ii) = latency(ii)+(Kcyc*TSS+qq*TRS)/NUE;
        end
    end
%     latency(ii) = latency_SS(ii)+latency_RS(ii);
    latency_tot(ii) = latency_IA+PMA*latency(ii);
end
figure;
% plot(OH*100,latency_tot/1e-3,'-o','linewidth',2);hold on
plot(SNR_range,latency_tot/1e-3,'-o','linewidth',2);hold on

grid on
xlabel('SNR (dB)')
ylabel('Access Latency (ms)')
