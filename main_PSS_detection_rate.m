clear;clc;

% ------ script control parameters -------
rng(2)
plot_cdf = 0;
MCtimes = 1e3; % Monte Carlo simulation
Pfa = 0.01; % target pfa in setting threshold
SNR_data_range = -30:1.5:-12;
SNR_range = SNR_data_range + 0; % NYU paper has SNR_data/SNR_IA mapping
SNR_num = length(SNR_range);

% parameter for evaluation of theoretical equations
Nt = 256;
Nr = 32;
BW = 2*28.8e6; % sample rate in [Hz]
Ts = 1/BW; % sample duration in [sec]
fc = 28; % carrier freq in [GHz]
P = 127; % subcarrier number in PPS
SC = 256; % subcarrier number in SS signal
OFDM_sym_num = 4; % num of OFDM in each SS burst
M = 64; % number of bursts
Nc = 4; % maximum multipath delay in [sample]
CP = 8; % cyclic prefix in [sample]
STO_max = 1000; % knowledge of maximum STO in [sample]

%% Monte Carlo evaluations
STO = 170; %[170 (3us), 568(10us), 960 (critical BF mismatch)]; Type of sample timing offset: value or 'zero'/'random'
if STO > SC*(OFDM_sym_num-1)
    K = SC*OFDM_sym_num-STO; 
else
    K = 0;
end % BF window mismatch, used in theo. eq.
BFtype = 'PN'; % Type of beamformer in IA: 'sector_LS','sector_FSM_KW' 'PN', or 'directional'
STOinfo = 0; % Assuming perfect knowledge of ZC-corr peak location or not
M_burst = [16, 4]; % Number of bursts in IA; For directional use [M_Tx_BF,M_Rx_BF] for beams in BS and UE, otherwise M equals to M_Tx_BF*M_Rx_BF
CFO = 5; % determinstic CFO in ppm
CFO_samp = CFO*fc*1e3*Ts*2*pi; % phase rototion per sample in [rad]; 1e3 due to ppm is 1e6 and GHz is 1e9;

% pre-computer sector beam if used
switch BFtype
    case 'PN'
        Tx_sec_codebook = [];
        Rx_sec_codebook = [];
    case 'sector_LS'
        Tx_sec_codebook = get_IA_BF(Nt, M_burst(1), BFtype); % Tx beamformer in IA stage
        Rx_sec_codebook = get_IA_BF(Nr, M_burst(2), BFtype); % Tx beamformer in IA stage
    case 'sector_FSM_KW'
        Tx_sec_codebook = get_IA_BF(Nt, M_burst(1), BFtype); % Tx beamformer in IA stage
        Rx_sec_codebook = get_IA_BF(Nr, M_burst(2), BFtype); % Tx beamformer in IA stage
end

% ------------ MC iterations (each has all SNRs)--------------
for MCindex = 1:MCtimes
    clc
    fprintf('Iteration %d:\n',MCindex);
%     if MCindex==2
%         debug_flag = 1;
%     end
    [ peak_pow_H1(:,MCindex),...
      peak_pow_H0(:,MCindex),...
      peakindex_H1(:,MCindex)] = run_PSS_detection_CFO( SNR_range,...
                                                        STO,...
                                                        STOinfo,...
                                                        BFtype,...
                                                        M_burst,...
                                                        CFO,...
                                                        Tx_sec_codebook,...
                                                        Rx_sec_codebook);
end

%% % --------- Detection based on emprical threshold --------
for ss = 1:SNR_num
    temp_H0 = sort(peak_pow_H0(ss,:),'descend');
    TH(ss) = temp_H0(MCtimes*Pfa);
    Pm_sim(ss) = sum(peak_pow_H1(ss,:)<TH(ss))/MCtimes;
    switch STOinfo
        case 0
            PTOm_sim(ss) = sum((peakindex_H1(ss,:)<(STO+1-CP))|(peakindex_H1(ss,:)>(STO+1))|(peak_pow_H1(ss,:)<TH(ss)))/MCtimes;
    end
    Pfa_sim(ss) = sum(peak_pow_H0(ss,:)>TH(ss))/MCtimes;
end

%% ---------- Detection statistic in H0, sim. vs theo. (for debug) ------------------
for ss = 1:SNR_num
    noise_pow = 10^(-SNR_range(ss)/10);
    
    % plot cdf of detection statistics (emprical)
    if plot_cdf
        figure
        [a,b] = ecdf(peak_pow_H0(ss,:));
        plot(b,a);hold on
        grid on
    end
    
    mu(ss) = Nc*noise_pow/127; % As derived in appendix
    sigma(ss) = noise_pow/127*sqrt(Nc/M_burst(1)/M_burst(2)); % As derived in appendix
    
    switch STOinfo
        case 1
            x = linspace(mu(ss)-4*sigma(ss),mu(ss)+4*sigma(ss),1e3);
%             y = gampdf(x,M_burst(1)/2,2*mu(ss)/M_burst(1)); % for small M_burst(1), gamma dist. should be a better fit but does not work out
            y = normpdf(x,mu(ss),sigma(ss));
        case 0
            mu_max(ss) = mu(ss) + sigma(ss)*qfuncinv(1/STO_max); % it seems a factor of 0.9 gives more fit
            sigma_max(ss) = sigma(ss)/qfuncinv(1/STO_max);
            x = linspace(mu_max(ss)-4*sigma_max(ss),mu_max(ss)+4*sigma_max(ss),1e3);
            y = evpdf(-x,-mu_max(ss),sigma_max(ss));
    end
    for xx = 1:length(x)
        theo_cdf(xx) = sum(y(1:xx))*(x(2)-x(1));
    end
    
    [~,TH_theo_index] = min(abs(theo_cdf-(1-Pfa)));
%     TH_theo0(ss) = x(TH_theo_index);

%     TH_theo(ss) = noise_pow*(Nc/P+sqrt(Nc/M)/P*qfuncinv(Pfa));
    % from derived equation (perfect timing)
    
    xi = qfuncinv(1/STO_max)-0.78*log(-log(1-Pfa))/qfuncinv(1/STO_max);
    TH_theo(ss) = noise_pow*(Nc/P+sqrt(Nc/M/P^2)*xi);
    % from derived equation (no timing)
    
    % plot emp. cdf of detection statistics and theo cdf in H0
    if plot_cdf
        plot(x,theo_cdf);hold on
        title(num2str(SNR_range(ss)))
        legend('sim.','theo.')
    end
end

if plot_cdf
    figure
    plot(TH);hold on
    plot(TH_theo);hold on
    legend('sim.','theo.')
    title('optimal threshold for target pfa')
end

%% --------- Detection based on theoretical threshold --------
% for ss = 1:SNR_num
%     Pm_sim(ss) = sum(peak_pow_H1(ss,:)<TH_theo(ss))/MCtimes;
%     PTOm_sim(ss) = sum((peakindex_H1(ss,:)<(STO+1-CP))|(peakindex_H1(ss,:)>(STO+1))|(peak_pow_H1(ss,:)<TH_theo(ss)))/MCtimes;
%     Pfa_sim(ss) = sum(peak_pow_H0(ss,:)>TH_theo(ss))/MCtimes;
% end

%% ------- Detection statistics in H1, sim. vs theo. (for debug) --------------
clearvars theo_cdf_H1
for ss = 1:SNR_num
    noise_pow = 10^(-SNR_range(ss)/10);
    [a,b] = ecdf(peak_pow_H1(ss,:));

    % plot cdf of detection statistics (emprical)
    if plot_cdf
        figure
        plot(b,a);hold on
        grid on
    end

    
    switch STOinfo
        case 1
            switch STO
                case 'zero'
                    mu_H1 = 1 + mu(ss); % it seems a factor of 0.95 gives more fit
                    sigma_H1 = sqrt(1/M + sigma(ss)^2); 
                    x = linspace(mu_H1-4*sigma_H1,mu_H1+4*sigma_H1,1e3);
                    y = normpdf(x,mu_H1,sigma_H1);
                    for xx = 1:length(x)
                        theo_cdf_H1(xx) = sum(y(1:xx))*(x(2)-x(1));
                    end
                    [~, H1_theo_index(ss)] = min(abs((x-TH_theo(ss))));
%                     Pm_theo(ss) = theo_cdf_H1(H1_theo_index(ss));
                    Pm_theo(ss) = qfunc((mu_H1-TH_theo(ss))/sigma_H1);
                case 'random'
                    mu_H1 = 0.6681 + mu_max(ss); % mean value of t^2-2t+1 when t \in [0,1]
                    sigma_H1 = sqrt((sqrt(11)/sqrt(127))^2 + (sigma_max(ss))^2); 
                    x = linspace(mu_H1-4*sigma_H1,mu_H1+4*sigma_H1,1e3);
                    y = normpdf(x,mu_H1,sigma_H1);
                    for xx = 1:length(x)
                        theo_cdf_H1(xx) = sum(y(1:xx))*(x(2)-x(1));
                    end
                    [~, H1_theo_index(ss)] = min(abs((x-TH_theo(ss))));
                    Pm_theo(ss) = theo_cdf_H1(H1_theo_index(ss));
            end
        case 0 % it's conservative to consider true correlation peak (detected peak is always higher than it!)
            
            kappa = (2-real(exp(1j*CFO_samp*K))-real(exp(1j*CFO_samp*(P-K))))/(P^2-P^2*real(exp(1j*CFO_samp)));
            mu_H1 = kappa*1 + mu(ss); % mean value of t^2-2t+1 when t \in [0,1]
            sigma_H1 = sqrt(3*kappa^2/M + (sigma(ss))^2); 
            x = linspace(mu_H1-4*sigma_H1,mu_H1+4*sigma_H1,1e3);
            y = normpdf(x,mu_H1,sigma_H1);
            for xx = 1:length(x)
                theo_cdf_H1(xx) = sum(y(1:xx))*(x(2)-x(1));
            end
            [theo_cdf_H1,x] = get_emp_maxrv_cdf(mu_H1,sigma_H1,mu_max(ss),sigma_max(ss));
            [~, H1_theo_index(ss)] = min(abs((x-TH_theo(ss))));
            Pm_theo(ss) = theo_cdf_H1(H1_theo_index(ss));
    end
    
    % plot cdf of detection statistics (theo)
    if plot_cdf
        plot(x,theo_cdf_H1)
        legend('sim.','theo.')
        title(num2str(SNR_range(ss)))
    end
end
%%
SNR_dB_data_fine = -30:0.1:0;
Pm_theo = zeros(length(SNR_dB_data_fine),1);
SNR_data_fine = 10.^(SNR_dB_data_fine/10);
switch STOinfo
    case 1
        switch STO
            case 'zero'
                for ss=1:length(SNR_data_fine)
                    SNR = SNR_data_fine(ss);
                    Pm_theo(ss) = qfunc((SNR - sqrt(Nc/M)/127*qfuncinv(Pfa))./sqrt(2*SNR.^2/M/1+Nc/127^2/M));
                end
        end
    case 0
        switch STO
            case 'random'
            otherwise
                kappa = (2-real(exp(1j*CFO_samp*K))-real(exp(1j*CFO_samp*(P-K))))/(P^2-P^2*real(exp(1j*CFO_samp)));
                xi = qfuncinv(1/STO_max)-0.78*log(-log(1-Pfa))/qfuncinv(1/STO_max);
                for ss=1:length(SNR_data_fine)
                    SNR = SNR_data_fine(ss);
                    Pm_theo(ss) = qfunc((kappa*SNR - sqrt(Nc/M/P^2)*xi)./sqrt(2*kappa^2*SNR.^2/M+Nc/P^2/M));
                end
        end
end

%% Figure
mycolor = [         0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250;
                0.4940    0.1840    0.5560];
figure
subplot(211)
    
switch STOinfo
    case 0
    semilogy(SNR_data_range,PTOm_sim,'-x','linewidth',2,'markersize',8);hold on
end

semilogy(SNR_data_range,Pm_sim,'o',...
            'linewidth',2,'markersize',8,'color',mycolor(2,:));hold on
% semilogy(SNR_data_range,Pm_theo,...
%               'linewidth',2,'color',mycolor(2,:));hold on
semilogy(SNR_dB_data_fine,Pm_theo,...
            'linewidth',2,'color',mycolor(2,:));hold on
grid on
ylim([0.001,1])
xlabel('SNR (dB)')
ylabel('Miss Detection of PSS')
legend('PTOM Sim.','PM Sim.','PM Theo.')
subplot(212)
plot(SNR_data_range,Pfa_sim,'-o');hold on
grid on
xlabel('SNR (dB)')
ylabel('False Alarm Rate')
ylim([0,0.1])




