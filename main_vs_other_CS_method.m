%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
path_num = 1; % Num of rays in a cluster
Nr = 32; % Number of antenna in Rx
Nt = 128;
M = 64; % Length of training
MCtimes = 20; % Num of Monte Carlo Sim.
AOAspread2 = 0;
AOAspread = 0;
AODspread2 = 0;
AODspread = 0;
SNR_num = 1;
SNR_range = linspace(80,80,SNR_num);
fc = 28e9;
BW = 2048e6; % IA bandiwdth
Ts = 1/BW; % Sample duration
P = 512; % number of subcarrier
df = BW/P;
DFT = dftmtx(P);

%-------- dictionary generation -------------
cand_num_r = 65;
cand_num_t = 255;
dict_num = cand_num_r*cand_num_t;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Nt-1)'*pi*sin(cand_angle_t));


% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCindex = 1:MCtimes
    
    clc; fprintf('iteration %d:\n',MCindex);
       
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);
    
    Measure_mat = kron(transpose(F)*conj(cand_ARV_t),W'*cand_ARV_r);
    select_row = zeros(1,M);
    for ii=1:M
        select_row(ii) = (ii-1)*M+ii;
    end
    Measure_mat_new = Measure_mat(select_row,:);
    for cc=1:cand_num_r*cand_num_t
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    end
    
    
%     probe_Tx_BF = ones(Nt,M);
%     F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);   

    % AoA of rays with disired seperation
    phi = zeros(path_num,1);
    phi0(MCindex) = 0/180*pi;%(rand*90-45)/180*pi;
    phi = phi0(MCindex) + randn(path_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta = zeros(path_num,1);
    theta0(MCindex) = cand_angle_t(35);%-43.9370/180*pi;%(rand*90-45)/180*pi; 
    % Index 25 in the dictionary
    theta = theta0(MCindex) + randn(path_num,1) * AODspread;
    
    % Find True AoA/AoD in grid (for debug)
    [~,row_true] = min(abs(cand_angle_r - phi0(MCindex)));
    [~,col_true] = min(abs(cand_angle_t - theta0(MCindex)));
    index_true = (col_true-1)*cand_num_r + row_true;

%     % Gain
%     g_cmplx = exp(1j*rand(ray_num,1)*2*pi)/sqrt(ray_num);
%     g = g_cmplx;
    % Rotate of ray
    tau = rand*(90e-9);
    tau_samp(MCindex) = tau/Ts*2*pi;
%     g_ray = (randn+1j*randn)/sqrt(2);
    g_ray = 1;%exp(1j*rand*2*pi);

    % Pre-compute some vectors/matrices in FIM
    for pathindex = 1:path_num

        % Spatial response in T/Rx
        arx(:,pathindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(pathindex)))/sqrt(Nr);       
        atx(:,pathindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(pathindex)))/sqrt(Nt);
        
        % Delay response and its derivative over tau
        fvec(:,pathindex) = exp(-1j * (0:P-1)' * tau_samp(MCindex) / P);
        
    end
    
    % Received signals (Using unit freq training symbol)
    tau_num = 500;
    delay_cand = linspace(0,100,tau_num)*1e-9/Ts*2*pi;
    for tt=1:tau_num
        delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P));
    end
    
    sig_rx = zeros(P*M,1);
    for ll=1:path_num
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            sig_freq = zeros(P,1);
            for pp = 1:P
                atx_new(:,pp) = exp(1j * pi *(1+(pp-1)*df/fc) * (0:Nt-1)' * sin(theta(ll)))/sqrt(Nt);
                sig_freq(pp) = g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx_new(:,pp))...
                *fvec(pp);
%                 sig_freq(pp) = conj(F(:,mm)'*atx_new(:,pp)) * fvec(pp);
            end
            
            sig_rx(index) = DFT'* sig_freq;
        end
    end
    
    awgn = (randn(P*M,1)+1j*randn(P*M,1))/sqrt(2);
    
    % ------------- For loop for Various SNR ----------------------
    for ss = 1:SNR_num
        
        % SNR and Adding AWGN to Rx Signal
        sigman2 = 10^(-SNR_range(ss)/10);
        sig_noisy = sig_rx + awgn * sqrt(sigman2);
        
        % ----- Benchmark Method #1 (OMP per Subcarrier) ---------
        sig_PtimeM = DFT * reshape(sig_noisy,P,M);
        for pp=1:round(P/16)
            SC_idx = (pp-1)*16+1;
            sig_for_OMP = sig_PtimeM(SC_idx,:).';
            for dd=1:dict_num
                score_OMP(dd) = abs(sig_for_OMP'*(Measure_mat_new(:,dd)))...
                    /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
            end
            [~,bestindex_temp] = max(abs(score_OMP));
            bestrow = floor((bestindex_temp-1)/cand_num_r)+1;
            bestcol = bestindex_temp-(bestrow-1)*cand_num_r;
            bestAOA_temp(pp) = (bestcol-1)*AOAstep-60*pi/180;
            bestAOD_temp(pp) = (bestrow-1)*AODstep-60*pi/180;
        end
%         bestAOA_BM = mean(bestAOA_temp);
%         bestAOD_BM = mean(bestAOD_temp);
        AOA_error_BM(MCindex,ss) = sqrt(mean(abs(bestAOA_temp - phi0(MCindex)).^2));
        AOD_error_BM(MCindex,ss) = sqrt(mean(abs(bestAOD_temp - theta0(MCindex)).^2));
        
        % --------Benchmark Method #2 (???) ---------------
        
        % -------- Genie Verification on Dictionary Adaptation ------
        sig_SC_sum = sum(sig_PtimeM, 1);
        Gamma2 = (1-(alpha2)^(P))/(1-alpha2);
        
        alpha_new(1) = alpha2;
        Gamma(1) = conj(Gamma2);
        for nn=1:Nt-1
            alpha_new(nn+1) = exp(1j*pi*(df/fc)*(nn)*sin(theta0(MCindex)))...
                    *exp(1j*tau_samp(MCindex)/P);
            
            Gamma(nn+1) = (1-(alpha_new(nn+1))^(P))/(1-alpha_new(nn+1));
        end

        for mm=1:M
            sig_model_verify(mm) = P * g_ray(1) * (W(:,mm)'*arx(:,1)) * conj(F(:,mm)'*(diag(Gamma)*atx(:,1)));
            sig_model_verify2(mm) = P*Gamma2*(g_ray(1) * (W(:,mm)'*arx(:,1)) * conj(F(:,mm)'*(atx(:,1))));

        end
        best_fit = pinv(sig_model_verify')*(sig_SC_sum');
%         figure
%         subplot(211)
%         plot(abs(sig_model_verify - sig_SC_sum))
%         hold on
%         subplot(212)
%         plot(abs(sig_SC_sum))
%         grid on
        
        % -------- Proposed Method w/ First Delay Matching Pursuit ------------
        sig_ave = mean(reshape(sig_noisy,P,M),2);
        for tt=1:tau_num
            score(tt) = abs(sum(sig_ave.*conj(delay_mtx(:,tt)))/norm(delay_mtx(:,tt))^2);
        end
        
        % Plug-in Matching Pursuit index for CFO and Delay est.
        [maxscore,maxindex] = max(score);
        maxindex_est(ss,MCindex) = maxindex;
        delay_est(ss,MCindex) = delay_cand(maxindex);
        
        
        % ---------------   AoA/AoD Estimation   ----------------
        sig_desymb = zeros(M,1);
        
        % Evaluate average subcarrier
        sig_ave_SC = (sum(DFT *reshape(sig_noisy,P,M),1)).';
        
        % Dictionary adaptation
%         alpha2 = exp(-1j*delay_est(ss,MCindex)/P);
        alpha2 = exp(-1j*tau_samp(MCindex)/P);

        Gamma2 = (1-(alpha2)^(P))/(1-alpha2);
        alpha_new(1) = alpha2;
        Gamma(1) = conj(Gamma2);
        for dd=8078
            rowdd = floor((dd-1)/cand_num_r)+1;
            coldd = dd-(rowdd-1)*cand_num_r;
            for nn=1:Nt-1
%                 alpha_new(nn+1) = exp(1j*pi*(df/fc)*(nn)*sin(cand_angle_t(rowdd)))...
%                         *exp(1j*delay_est(ss,MCindex)/P);
                alpha_new(nn+1) = exp(1j*pi*(df/fc)*(nn)*sin(cand_angle_t(rowdd)))...
                        *exp(1j*tau_samp(MCindex)/P);
                    

                Gamma(nn+1) = (1-(alpha_new(nn+1))^(P))/(1-alpha_new(nn+1));
            end
            
            dict_adapt = zeros(M,1);
            arx_dd = exp(1j * pi * (0:Nr-1)' * sin(cand_angle_r(coldd)))/sqrt(Nr);       
            atx_dd = exp(1j * pi * (0:Nt-1)' * sin(cand_angle_t(rowdd)))/sqrt(Nt);
            for mm=1:M
                dict_adapt(mm) = P*(W(:,mm)'*arx_dd) * conj(F(:,mm)'*(diag(Gamma)*atx_dd));
            end
            best_fit = pinv(dict_adapt)*(sig_ave_SC)
            
            score_final(dd) = abs(sig_ave_SC'* dict_adapt)/(dict_adapt'*dict_adapt);
        end
        
%         % Old Matching pursuit for all AoA/AoD pair in dictionary
%         for dd=1:dict_num
%             
%             score_final(dd) = abs(sig_desymb'*(Measure_mat_new(:,dd)))...
%                 /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
%         end
        [~,bestindex_comp(MCindex)] = max(abs(score_final));
        
%         bestindex_comp(MCindex) = index_true;% debug. comment in main script
        bestrow = floor((bestindex_comp(MCindex)-1)/cand_num_r)+1;
        bestcol = bestindex_comp(MCindex)-(bestrow-1)*cand_num_r;
        bestAOA(MCindex,ss) = (bestcol-1)*AOAstep-60*pi/180;
        bestAOD(MCindex,ss) = (bestrow-1)*AODstep-60*pi/180;
        
        phi_last = bestAOA(MCindex,ss);
        theta_last= bestAOD(MCindex,ss);
        
        AOA_error_nocomp(MCindex,ss) = abs(phi_last - phi0(MCindex));
        AOD_error_nocomp(MCindex,ss) = abs(theta_last - theta0(MCindex));
        
    end
end
%%
for ss=1:SNR_num
AOAalign_comp_mean(ss) = sum((AOA_error_nocomp(:,ss)/pi*180)<(105/Nr),1)/MCtimes;
AODalign_comp_mean(ss) = sum((AOD_error_nocomp(:,ss)/pi*180)<(105/Nt),1)/MCtimes;
Align_comp_mean(ss) = sum(((AOD_error_nocomp(:,ss)/pi*180)<(105/Nt)&...
                           (AOA_error_nocomp(:,ss)/pi*180)<(105/Nr)),1)/MCtimes;

end

figure
% plot(SNR_range,AOAalign_comp_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_comp_mean,'--','linewidth',2);hold on

semilogy(SNR_range,1-AOAalign_comp_mean,'-','linewidth',2);hold on
semilogy(SNR_range,1-AODalign_comp_mean,'-','linewidth',2);hold on
semilogy(SNR_range,1-Align_comp_mean,'-','linewidth',2);hold on

grid on
xlabel('SNR (dB)')
ylabel('Misalignment Rate')
legend('AoA','AoD','Alignment')
%% RMSE evaluation
portion = 0.9; % for debug 
for ss=1:SNR_num
    delay_error = abs(delay_est(ss,:)-tau_samp).^2;
    RMSE_delay(ss) = sqrt(mean(get_portion(delay_error,portion)))*(Ts/(2*pi)/1e-9);
%     RMSE_delay(ss) = sqrt(mean((delay_error)))*(Ts/(2*pi)/1e-9);

    
    input_sig = abs(AOD_error_nocomp(:,ss)).^2;
    RMSE_theta(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
%     RMSE_theta(ss) = sqrt(mean((input_sig)))*(1/pi*180);

    
    input_sig = abs(AOA_error_nocomp(:,ss)).^2;
%     RMSE_phi(ss) = sqrt(mean((input_sig)))*(1/pi*180);
    RMSE_phi(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
    
    input_sig = abs(AOD_error_BM(:,ss)).^2;
    RMSE_theta_BM(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
    
    input_sig = abs(AOA_error_BM(:,ss)).^2;
    RMSE_phi_BM(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);


end
figure
subplot(211)
semilogy(SNR_range,RMSE_delay);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of normalized delay [ns]')
title('delay est')


subplot(212)
semilogy(SNR_range,RMSE_theta,'linewidth',2);hold on
semilogy(SNR_range,RMSE_phi,'linewidth',2);hold on
semilogy(SNR_range,RMSE_theta_BM,'linewidth',2);hold on
semilogy(SNR_range,RMSE_phi_BM,'linewidth',2);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA/AoD [deg]')
title('angle est')
legend('AoD','AoA','AoD (BM)','AoA (BM)')

%%
for ss=1:length(SNR_range)
    for MCindex = 1:MCtimes
        AOA_error_nofine(MCindex,ss) = abs(bestAOA(MCindex,ss) - phi0(MCindex));
        AOD_error_nofine(MCindex,ss) = abs(bestAOD(MCindex,ss) - theta0(MCindex));
    end
        
        input_sig = abs(AOD_error_nofine(:,ss)).^2;
        RMSE_theta_nofine(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
%     RMSE_theta(ss) = sqrt(mean((input_sig)))*(1/pi*180);

    
        input_sig = abs(AOA_error_nofine(:,ss)).^2;
%     RMSE_phi(ss) = sqrt(mean((input_sig)))*(1/pi*180);
        RMSE_phi_nofine(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
end

figure
semilogy(SNR_range,RMSE_theta_nofine);hold on
semilogy(SNR_range,RMSE_phi_nofine);hold on
grid on
xlabel('SNR [dB]')
ylabel('RMSE of AoA/AoD [deg]')
title('angle est w/o refinement')

