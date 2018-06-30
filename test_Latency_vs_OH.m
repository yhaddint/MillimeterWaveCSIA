clear;clc;

Tb = 74.6e-6; % Constant
BWSS = 57.4e6; % Constant
BWtot = 250e6; % Constant
TRS = 1.25e-3;
NUE_range = 1:50;


PMD = 0.15;
PMA = 0.15;

% TSS_range = [10,15,20,40,80,160]*1e-3;
TSS = 20e-3
M_burst = 64; % Varying
latency = zeros(length(NUE_range),1);
for ii=1:length(NUE_range)

%     TSS = TSS_range(ii); % Varying 
    NUE = NUE_range(ii); % Varying
    
    NRS = max(floor((TSS-Tb*M_burst)/TRS),1); % adaptively adjusted
    
    % compute overhead
    OH(ii) = (Tb*BWSS*M_burst)/(TSS*BWtot);
    K=0:5;
    latency_IA = sum(PMD*((1-PMD).^K).*K*TSS);
    Kcyc(ii) = floor((NUE-1)/NRS);
    Kres = NUE-Kcyc(ii)*NRS;
    if Kcyc(ii)==0
        latency(ii) = 0;
        for qq=1:NUE
            latency(ii) = latency(ii)+(qq*TRS)/NUE;
        end
    else
        for kk=1:Kcyc(ii)
            for qq=1:NRS
                latency(ii) = latency(ii)+((kk-1)*TSS+qq*TRS)/NUE;
            end
        end
        for qq=1:Kres
            latency(ii) = latency(ii)+(Kcyc(ii)*TSS+qq*TRS)/NUE;
        end
    end
%     latency(ii) = latency_SS(ii)+latency_RS(ii);
    latency_tot(ii) = latency_IA+PMA*latency(ii);
end
figure;
% plot(OH*100,latency_tot/1e-3,'-o','linewidth',2);hold on
plot(NUE_range,latency_tot/1e-3,'-o','linewidth',2);hold on

grid on
xlabel('Connected Number of UEs')
ylabel('Access Latency (ms)')
