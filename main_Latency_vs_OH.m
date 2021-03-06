clear;clc;

Tb = 74.6e-6; % Constant
BWSS = 57.4e6; % Constant
BWtot = 250e6; % Constant
TRS = 1.25e-3;
NUE_range = 1:50;

TRS_range = 2.^(-1*(linspace(1,6,20)))*10*1e-3;
% TRS_range = 2.^(-1*(1:6))*10*1e-3;

% TSS_range = [10,15,20,40,80,160]*1e-3;
TSS = 20e-3
BWRS = BWtot
DRS = 50e-6;
M_burst = 64; % Varying
latency = zeros(length(NUE_range),1);

% PMD = 0.05;
% PMA = 0.05;

PMD = 0.05;
PMA = 0.7; 

for ii=1:length(TRS_range)

%     TSS = TSS_range(ii); % Varying 
    NUE = 20; % Varying
    
    TRS = TRS_range(ii);
    NRS = max(floor((TSS-Tb*M_burst)/TRS),1); % adaptively adjusted
    
    % compute overhead
    OH(ii) = (Tb*BWSS*M_burst + NRS*DRS*BWRS)/(TSS*BWtot);
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
%%
figure;
% plot(OH*100,latency_tot/1e-3,'-o','linewidth',2);hold on
plot(OH*100,latency_tot/1e-3,'linewidth',2);hold on

grid on
xlabel('Overhead Ratio (%)')
ylabel('Access Latency (ms)')
