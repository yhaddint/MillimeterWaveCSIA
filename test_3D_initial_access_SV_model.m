%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3);                                     %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
path_num = 1;                               % Num of rays in a cluster

fc = 28e9;                                  % carrier freq [Hz]

Nr = 8;                                     % Num of antenna in UE/Rx (total)
Nr_az = 8;                                  % Num of antenna in UE/Rx (azimuth)
Nr_el = 1;                                  % Num of antenna in UE/Rx (elevation)

Nt = 128;                                   % Num of antenna in BS/Tx (total)                           
Nt_az = 16;                                 % Num of antenna in BS/Tx (azimuth) 
Nt_el = 8;                                  % Num of antenna in BS/Tx (elevation) 

M = 64;                                     % Length of SS bursts (IA OFDM symbols)
MCtimes = 5e1;                             % Num of Monte Carlo Sim.

AOAspread2 = 0;                             % Intra-cluster AoA spread square 
AOAspread = 0;                              % Intra-cluster AoA spread RMS 
AODspread2 = 0;                             % Intra-cluster AoD spread RMS 
AODspread = 0;                              % Intra-cluster AoD spread RMS 
SNR_num = 1;                               % Num of Monte Carlo Sim.
SNR_range = linspace(30,30,SNR_num);
BW = 57.6e6;                                % IA bandiwdth [Hz]
Ts = 1/BW;                                  % Sample duration
Nb = 512;                                   % Sample per SS burst
CFO_ppm = 0;                                % CFO in ppm
CFO = (fc/1e6*CFO_ppm);                     % CFO with unit Hz
eF = CFO*Ts*2*pi;                           % CFO normalized with Ts
P = 128;                                    % Number of subcarrier for PSS
DFT = dftmtx(P);
to_est_CFO = 1;                             % Assuming perfect knowledge of CFO or not
max_ite_num = 1e3;                          % Number iteration in refinement steps
refine_CFO = 1;                             % Turn on refine CFO when coarse estimation of AoA/AoD (long long time!!!)

az_lim = pi/3;
el_lim = pi/6;

%-------------------------------------
% Phase Noise Specification
% FOM = L(f0,df) + 20log10(df/f0)+10log10(P_VCO)
%-------------------------------------
VCO_FOM = -114 + 20*log10(1e6/fc)...
            + 10*log10(27);                 % phase noise FOM -114dBc@1MHz w. 27mW
P_VCO = 500;                                % Scaling PN var by changing VCO power [mW]
VCO_c = 10^(VCO_FOM/10)/P_VCO;              % parameter c in PN specs
PN_sigma2 = VCO_c*4*pi^2*fc^2*Ts;           % Time domain variance of PN Wiener process
PN_sigma = sqrt(PN_sigma2);                 % Weiner process RMS imcrement


%-------- CS dictionary generation -------------
cand_num_r_az = 17;                         % Grid size for AoA Azimuth (2 times Nr_az)
cand_num_r_el = 1;                          % Grid size for AoA Elevation (1 by default)

cand_num_t_az = 33;                         % Grid size for AoD Azimuth (2 times Nt_az)
cand_num_t_el = 17;                         % Grid size for AoD elevation (2 times Nt_el)

dict_num = cand_num_r_az * cand_num_r_el * cand_num_t_az * cand_num_t_el;

cand_y = zeros(M, dict_num);

cand_angle_r_az = linspace(-az_lim, az_lim, cand_num_r_az);
AOAstep_az = cand_angle_r_az(2)-cand_angle_r_az(1);

cand_angle_t_el = linspace(-el_lim, el_lim, cand_num_t_el);
AODstep_el = cand_angle_t_el(2)-cand_angle_t_el(1);

cand_angle_t_az = linspace(-az_lim, az_lim, cand_num_t_az);
AODstep_az = cand_angle_t_az(2)-cand_angle_t_az(1);

cand_ARV_r_az = exp(1j*(0:Nr_az-1)'*pi*sin(cand_angle_r_az));
cand_ARV_r = cand_ARV_r_az;

cand_ARV_t_az = exp(1j*(0:Nt_az-1)'*pi*sin(cand_angle_t_az));
cand_ARV_t_el = exp(1j*(0:Nt_el-1)'*pi*sin(cand_angle_t_el));

% Planar geometry is [ANT_row1.'; ANT_row2.' cdots, ANT_row_Nr_el.']
cand_ARV_t = kron(cand_ARV_t_el, cand_ARV_t_az);

% test scenario when assuming CFO is known
phase_error_mat = kron(exp(1j*eF*Nb*(0:M-1)),exp(1j*eF*(0:P-1)'));
phase_error = reshape(phase_error_mat,M*P,1);


% ------- Precompute codebook for Sector Search -------------
% number of sector beams in steering
M_BS_burst_az = 8;
M_BS_burst_el = 2;
M_UE_burst_az = 4;

% angle grid for sector beams; sector approach gives estimator from grid
BS_az_grid = (az_lim)*linspace(-1+(1/M_BS_burst_az),1-(1/M_BS_burst_az),M_BS_burst_az);
BS_el_grid = (el_lim)*linspace(-1+(1/M_BS_burst_el),1-(1/M_BS_burst_el),M_BS_burst_el);
UE_az_grid = (az_lim)*linspace(-1+(1/M_UE_burst_az),1-(1/M_UE_burst_az),M_UE_burst_az);

% codebook for sector IA approach
W_sec_mat = get_IA_BF(Nr_az, M_UE_burst_az, 'sector_FSM_KW'); % Rx BF in IA stage
F_sec_mat = get_IA_BF_3D(Nt_az, Nt_el,...
                     M_BS_burst_az, M_BS_burst_el,...
                     'sector_FSM_KW', az_lim, el_lim); % Tx beamformer in IA stage


% ---- signal length parameters --------
ZC_root = 29; % ZC root, a coprime number with ZC length
ZC_N = 127; % ZC sequence length
seq = lteZadoffChuSeq(ZC_root,ZC_N);
                
debug_flag = 0;           

% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCidx = 1:MCtimes
    
    clc
    fprintf('iteration %d:\n',MCidx);

    % ------- Compressive Approach -------------
%     debug_flag = 0;
%     if MCindex==1
%         debug_flag = 1;
%     end
    
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
    for cc=1:cand_num_r_az*cand_num_t_az
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    end
    
    % weighted average
    for dd=1:dict_num
        index_new = [abs(Measure_mat_new(1:M-1,dd)),abs(Measure_mat_new(2:M,dd))];
        CFO_weight = min(index_new,[],2);
        CFO_weight_norm = CFO_weight/sum(CFO_weight);
        CFO_select(:,dd) = CFO_weight_norm>1/M;
    end

    
%     probe_Tx_BF = ones(Nt,M);
%     F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);   

    % AoA of rays with disired seperation
    phi_az = zeros(path_num,1);
    phi0_az(MCidx) = 0.2618;%(rand*90-45)/180*pi;
    phi_az = phi0_az(MCidx) + randn(path_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta_az = zeros(path_num,1);
    theta_el = zeros(path_num,1);

    theta0_az(MCidx) = 0.1309;%(rand*90-45)/180*pi;
    theta0_el(MCidx) = 0.2618;%(rand*30-15)/180*pi;

    theta_az = theta0_az(MCidx) + randn(path_num,1) * AODspread;
    theta_el = theta0_el(MCidx) + randn(path_num,1) * AOAspread;

    
    % Find True AoA/AoD in grid (for debug)
    [~,row_true] = min(abs(cand_angle_r_az - phi0_az(MCidx)));
    [~,col_az_true] = min(abs(cand_angle_t_az - theta0_az(MCidx)));
    [~,col_el_true] = min(abs(cand_angle_t_el - theta0_el(MCidx)));
    
    index_true = ((col_el_true-1) * cand_num_t_az + col_az_true - 1) * cand_num_r_az + row_true;

%     % Gain
%     g_cmplx = exp(1j*rand(ray_num,1)*2*pi)/sqrt(ray_num);
%     g = g_cmplx;
    % Rotate of ray
    tau = rand*(90e-9);
    tau_samp(MCidx) = tau/Ts*2*pi;
%     g_ray = (randn+1j*randn)/sqrt(2);
    g_ray = exp(1j*rand*2*pi);

    % Pre-compute some vectors/matrices in FIM
    for pathindex = 1:path_num

        % Spatial response and its derivative over phi
        arx_az(:,pathindex) = exp(1j * pi * (0:Nr_az-1)' * sin(phi_az(pathindex)))/sqrt(Nr_az);
        arx(:,pathindex) = arx_az(:,pathindex);

        % Spatial response and its derivative over theta
        atx_az(:,pathindex) = exp(1j * pi * (0:Nt_az-1)' * sin(theta_az(pathindex)))/sqrt(Nt_az);
        atx_el(:,pathindex) = exp(1j * pi * (0:Nt_el-1)' * sin(theta_el(pathindex)))/sqrt(Nt_el);
        atx(:,pathindex) = kron( atx_el(:,pathindex), atx_az(:,pathindex) );
        
        % Delay response and its derivative over tau
        fvec(:,pathindex) = exp(-1j * (0:P-1)' * tau_samp(MCidx) / P);
%         Dfvec(:,pathindex) = -1j * (0:P-1)'/P;
%         dfvec(:,pathindex) = Dfvec(:,pathindex).*fvec(:,pathindex);
    end
    
    % ------ Phase noise process -----------
    PN_seq = zeros(Nb * M, 1);
    PN_seq(1) = 1;
    for ll=1:length(PN_seq)
        PN_seq(ll+1) = PN_seq(ll).*exp(1j * randn * PN_sigma);
    end
    
    
    % About CFO and its derivative
    qvec = exp(1j * (0:P-1)' * eF);
    for mm=1:M
        timeindex = ((mm-1)*Nb+0):((mm-1)*Nb+P-1);
        Dqvec(:,mm) = 1j * (timeindex).';
        dqvec(:,mm) = exp(1j*Nb*(mm-1)*eF)*qvec.*Dqvec(:,mm);
    end
    
    
    % Pre-compute some equations, vecx, vecalpha, vecbeta
    % vecx is derivative of phi_0 in f(r,eta)
    vaa = zeros(M,1);
    vad = zeros(M,1);
    vda = zeros(M,1);
    
    % Zero initialization of vectors for FIM
    vdcfo = zeros(P*M,1);
    vdtheta = zeros(P*M,1);
    vdphi = zeros(P*M,1);
    vdtau = zeros(P*M,1);
    vdalpha = zeros(P*M,1);
    vdbeta = zeros(P*M,1);
    
    % Received signals (Using random symbol or ZC sequence)
    symb = [seq;1]; %exp(1j*rand(P,1)*2*pi);
    tau_num = 500;
    delay_cand = linspace(0,100,tau_num)*1e-9/Ts*2*pi;
    for tt=1:tau_num
        delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P).*symb);
    end
    

    sig_rx = zeros(P*M, 1);
    for ll=1:path_num
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            indexPN = (mm-1)*Nb+1:(mm-1)*Nb+P;
            sig_rx(index) = (g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*eF)).*PN_seq(indexPN);
        end
    end
    
    awgn = (randn(P*M,1)+1j*randn(P*M,1))/sqrt(2);
    
    % ------- Precompute for Sector Search -------------
    sig_rx_sec = zeros(P*M, 1);
    for ll=1:path_num
        for mm=1:M
            
            mm_BS = floor((mm-1)/M_UE_burst_az)+1;
            mm_UE = mm - (mm_BS-1)*M_UE_burst_az;
            
            index = (mm-1)*P+1:mm*P;
            indexPN = (mm-1)*Nb+1:(mm-1)*Nb+P;
            sig_rx_sec(index) = (g_ray(ll) * (W_sec_mat(:,mm_UE)'*arx(:,ll)) * conj(F_sec_mat(:,mm_BS)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*eF)).*PN_seq(indexPN);
        end
    end

    
    % ------------- For loop for Various SNR ----------------------
    for ss = 1:SNR_num
        
        % SNR and Adding AWGN to Rx Signal
        sigman2 = 10^(-SNR_range(ss)/10);
        sig_noisy = sig_rx + awgn * sqrt(sigman2);
        
        % -------- Sector Approach -------
        sig_sec_ave = mean(abs(reshape(sig_rx_sec + awgn * sqrt(sigman2), P, M)).^2,1);
        [~,best_sec_idx] = max(sig_sec_ave);
        
        % Use SS burst index to get sector index in UE, BS_az, and BS_el
        best_sec_BS_idx = floor((best_sec_idx-1)/M_UE_burst_az)+1;
        best_sec_UE_idx = best_sec_idx - (best_sec_BS_idx-1)*M_UE_burst_az;
        best_sec_BS_el_idx = floor((best_sec_BS_idx-1)/M_BS_burst_az)+1;
        best_sec_BS_az_idx = best_sec_BS_idx - (best_sec_BS_el_idx-1)*M_BS_burst_az;
        
        % use sector index to get propagation angle etimation
        sec_UE_az_est = UE_az_grid(best_sec_UE_idx);
        sec_BS_az_est = BS_az_grid(best_sec_BS_az_idx);
        sec_BS_el_est = BS_el_grid(best_sec_BS_el_idx);
        
        % Error Evaluation
        AOA_az_error_sec(MCidx,ss) = abs(sec_UE_az_est - phi0_az(MCidx));
        AOD_az_error_sec(MCidx,ss) = abs(sec_BS_az_est - theta0_az(MCidx));
        AOD_el_error_sec(MCidx,ss) = abs(sec_BS_el_est - theta0_el(MCidx));
        
        
        % -------- CFO Estimation and Delay Matching Pursuit ------------
        sig_ave = mean(reshape(sig_noisy,P,M),2);
        for tt=1:tau_num
%             if MCindex==1 & tt==6
%                 apple=1;
%             end
            
            % If selected delay candidate is true, sig.*conj(p_kx) is like
            % tone signal
            sig_tone = sig_ave.*conj(delay_mtx(:,tt))./abs(delay_mtx(:,tt)).^2;
            
%             % Method1 ML estimation of CFO; Kay's paper on freq. est. of tone
            for pp = 1:P-1
                CFO_hat(pp) = angle(sig_tone(pp+1).*conj(sig_tone(pp)));
            end
            CFO_est = sum(CFO_hat)/(P-1);
            
            % Method2 ML estimation of CFO; Kay's paper on freq. est. of tone
            for pp = 1:P-1
                CFO_hat(pp) = sig_tone(pp+1).*conj(sig_tone(pp));
            end
            CFO_est = angle(sum(CFO_hat)/(P-1));
            CFO_est = 0;
            
            % Use the ML CFO to evaluate efficacy of delay candidate
            CFO_est1 = angle(mean(exp(1j*CFO_hat)));
%             CFO_est = eF; % For debug; plug-in true CFO
            sig_deCFO = sig_ave.*exp(-1j*CFO_est*(0:P-1).');
            score(tt) = abs(sum(sig_deCFO.*conj(delay_mtx(:,tt)))/norm(delay_mtx(:,tt))^2);
        end
        
        % Plug-in Matching Pursuit index for CFO and Delay est.
        [maxscore,maxindex] = max(score);
        maxindex_est(ss,MCidx) = maxindex;
        delay_est(ss,MCidx) = delay_cand(maxindex);
        
        % watch score to debug
%         figure;plot(delay_cand/2/pi*Ts/1e-9,score)
        
        % ---------------   AoA/AoD Estimation   ----------------
        sig_desymb = zeros(M,1);
        
        % Evaluate average sample over each SS burst
        for mm=1:M
            indextemp = (mm-1)*P+1:mm*P;
            sig_desymb(mm) = mean(sig_noisy(indextemp).*(conj(delay_mtx(:,maxindex))...
                ./abs(delay_mtx(:,maxindex)))).';
        end
        
        % Use true CFO (for debug)
        if to_est_CFO==0
            phase_error_mat = exp(1j*eF*Nb*(0:M-1).');
            phase_error = phase_error_mat;
        end
        
        % Matching pursuit for all AoA/AoD pair in dictionary
        for dd=1:dict_num
            
            % Debug flag, used by seting AoA=AoD = 0
%             if MCindex==26 & dd==3691
%                 apple=1;
%             end
            
%             sig_cand = zeros(P*M,1);
%             for mm=1:M
%                 indextemp = (mm-1)*P+1:mm*P;
%                 sig_cand(indextemp) = delay_mtx(:,maxindex)*Measure_mat_new(mm,dd);
%             end
%             sig_cand = kron(Measure_mat_new(:,dd),delay_mtx(:,maxindex));
            
            % use estimated CFO or perfect CFO to comp phase error

            if to_est_CFO
                
                sig_burst = sig_desymb.*conj(Measure_mat_new(:,dd))./abs(Measure_mat_new(:,dd));
                
                % adjust N/A numbers
                sig_burst(isnan(sig_burst))=0;
                
%                 %  ---------  Quasi-ML Method 1  ---------
%                 CFO_hat_new1 = zeros(M-1,1);
%                 for mm=1:M-1
%                     if CFO_select(mm,dd)
%                     CFO_hat_new1(mm) = angle(sig_burst(mm+1).*conj(sig_burst(mm)));
%                     end
%                 end
%                 CFO_est = sum(CFO_hat_new1)/sum(CFO_select(:,dd))/Nb;

%                 
                %  ---------  Quasi ML Method 2  ---------
                CFO_hat_new2 = zeros(M-1,1);
                for mm=1:M-1
                    if CFO_select(mm,dd)
                        CFO_hat_new2(mm) = sig_burst(mm+1).*conj(sig_burst(mm));
                    end
                end
                CFO_est = angle(sum(CFO_hat_new2(CFO_select(:,dd))./abs(CFO_hat_new2(CFO_select(:,dd))))/sum(CFO_select(:,dd)))/Nb;
                
                % It seems refine CFO is necessary at beginning
                if refine_CFO
                    CFO_range = linspace(CFO_est*0.7,CFO_est*1.4,10);
                    score_CFO = zeros(10,1);
                for zz=1:10
%                     CFO_estold = CFO_est;
%                     phase_error_mat = exp(1j*CFO_estold*Nb*(0:M-1).');
%                     x_vec_CFO = Measure_mat_new(:,dd).*phase_error_mat;
%                     y_vec_CFO = sig_desymb;
%                     opt_alpha = 1;
%                     error_CFO_est =  y_vec_CFO - opt_alpha*x_vec_CFO;
%                     sig_deF = Measure_mat_new(:,dd).*(1j*(Nb*(0:M-1).').*exp(1j*CFO_estold*Nb*(0:M-1).'));           
%                     Q_cal = [real(sig_deF);imag(sig_deF)];
%                     tr = [real(error_CFO_est);imag(error_CFO_est)];
%                     deF = pinv(Q_cal) * tr;
%                     CFO_est = CFO_estold + deF;
                    phase_est_test = exp(1j*CFO_range(zz)*Nb*(0:M-1).');
                    score_CFO(zz) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_est_test));
%                     phase_error_mat = exp(1j*CFO_estold*Nb*(0:M-1).');
                end
                [~,CFO_index] = max(score_CFO);
                CFO_est = CFO_range(CFO_index);
                end

                % ----- true CFO (for debug) ---------
%                 CFO_est = eF;
                
                %  --------- Simpler approach of estimating CFO  (no good)---------
%                 CFO_est = angle(sum(sig_burst(find(CFO_select(:,dd)+1)))...
%                     *conj(sum(sig_burst(find(CFO_select(:,dd))))))/Nb;

                %  --------- simpler approach 2 (no good) ---------
%                 pad = 511;
%                 x_temp = sig_burst.*conj(Measure_mat_new(:,dd));
%                 CFO_est = max(abs(fft([x_temp;zeros(pad*M,1)])))*(2*pi)/((pad+1)*Nb);
                
                phase_error_mat = exp(1j*CFO_est*Nb*(0:M-1).');
                phase_error = phase_error_mat;
                CFO_final(dd) = CFO_est;
            end
            score_final(dd) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_error))...
                /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
        end
        [~,bestindex_comp(MCidx)] = max(abs(score_final));
        
%         bestindex_comp(MCindex) = index_true;% debug. comment in main script
        bestrow = floor((bestindex_comp(MCidx)-1)/cand_num_r_az)+1;
        bestcol = bestindex_comp(MCidx)-(bestrow-1)*cand_num_r_az;
        
        bestrow_el = floor((bestrow-1)/cand_num_t_az)+1;
        bestrow_az = bestrow - (bestrow_el-1)*cand_num_t_az;
        
        bestAOA_az(MCidx,ss) = (bestcol-1)*AOAstep_az-60*pi/180;
        
        bestAOD_az(MCidx,ss) = (bestrow_az-1)*AODstep_az-60*pi/180;
        bestAOD_el(MCidx,ss) = (bestrow_el-1)*AODstep_el-30*pi/180;
        
        % Plot score for debug
        if 0%debug_flag
            figure;plot(abs(score_final));hold on;plot(index_true,abs(score_final(index_true)),'o','linewidth',2,'markersize',10);
            grid on
            legend('scores','true AoA/AoD pair')
            
            CFO_est_range = linspace(0.0010,0.0020,200);
            for xx=1:200
                CFO_est = CFO_est_range(xx) ;
            phase_error_mat = exp(1j*CFO_est*Nb*(0:M-1).');
%             phase_error_mat = exp(1j*eF*Nb*(0:M-1).');
            phase_error = phase_error_mat;
            score_test_CFO(xx) = abs(sig_desymb'*(Measure_mat_new(:,dd).*phase_error));
            end
            
            figure
            plot(CFO_est_range,score_test_CFO)

        end
        
%         if debug_flag
%             figure;
%             subplot(211)
%             plot(CFO_hat_new1(CFO_select(:,dd))/Nb/eF)
%             title('single sample pair angle est')
%             grid on
%             subplot(212);
%             plot(abs(CFO_hat_new2(CFO_select(:,dd))))
%             title('sample mag')
%             grid on
%         end
        
        % -------- Refinement to improve accuracy ----------------
%         
%         ii = 0;
%         phi_hat = zeros(max_ite_num+1,1);
%         theta_hat = zeros(max_ite_num+1,1);
%         eF_hat = zeros(max_ite_num+1,1);
%         alpha_hat = zeros(max_ite_num+1,1);
%         
        phi_az_hat(1) = bestAOA_az(MCidx,ss);
        theta_az_hat(1)= bestAOD_az(MCidx,ss);
        theta_el_hat(1)= bestAOD_el(MCidx,ss);
        eF_hat(1) = (CFO_final(bestindex_comp(MCidx))*Nb+2*pi)/Nb;
%         stop_sign = 1;
%         while stop_sign==0
%             ii = ii + 1;
%             %---------------------------------------------
%             % Alpha estimation using previous coeff.
%             %---------------------------------------------
%             arx_hat = exp(1j*(0:Nr-1)'*pi*sin(phi_hat(ii)))/sqrt(Nr);
%             atx_hat = exp(1j*(0:Nt-1)'*pi*sin(theta_hat(ii)))/sqrt(Nt);
%             phase_error_refine = kron( exp(1j*eF_hat(ii)*Nb*(0:M-1)).',...
%                                        exp(1j*eF_hat(ii)*(0:P-1).'));
% 
%             H_cal = diag(W' * arx_hat * atx_hat' * F);
%             sig_alpha = kron(H_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             alpha_hat(ii) = pinv(sig_alpha) * sig_noisy;
%             error(ii) = norm(alpha_hat(ii)*sig_alpha - sig_noisy);
%             
%             % determine when to stop
%             if ii>1
%                 if error(ii) > error(ii-1)
%                     stop_sign = 1;
%                 end
%                 if abs(error(ii)-error(ii-1))/abs(error(ii-1))<1e-4
%                     stop_sign = 1;
%                 end
%                 if ii > max_ite_num
%                     stop_sign = 1;
%                 end
%             end
%             
%             
%             %---------------------------------------------
%             % CFO estimation using previous coeff.
%             %---------------------------------------------
%             deF = 1j*kron(Nb*(0:M-1).',(0:P-1).')...
%                 .*kron(exp(1j*eF_hat(ii)*Nb*(0:M-1).'),exp(1j*eF_hat(ii)*(0:P-1).'));
% 
%             H_cal = alpha_hat(ii) * diag(W' * arx_hat * atx_hat' * F);
%             sig_deF = kron(H_cal,delay_mtx(:,maxindex)).*deF;
%             sig_eF = kron(H_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             
%             yr = sig_noisy - sig_eF;
%             Q_cal = [real(sig_deF);imag(sig_deF)];
%             tr = [real(yr);imag(yr)];
% 
%             deF = pinv(Q_cal) * tr;
%             eF_hat(ii+1) = eF_hat(ii) + deF;
%             
%             phase_error_refine = kron( exp(1j*eF_hat(ii+1)*Nb*(0:M-1)).',...
%                                        exp(1j*eF_hat(ii+1)*(0:P-1).'));
%             
%             %---------------------------------------------
%             % Phi estimation using previous coeff.
%             %---------------------------------------------      
%             dPhi = exp(1j*pi*(0:Nr-1).'*sin(phi_hat(ii)))/sqrt(Nr).*(1j*pi*(0:Nr-1).'*cos(phi_hat(ii)));
%             
%             D_cal = alpha_hat(ii) * diag(W' * dPhi * atx_hat' * F);
%             H_cal = alpha_hat(ii) * diag(W' * arx_hat * atx_hat' * F);
%             sig_dphi = kron(D_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             sig_phi = kron(H_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             
%             yr = sig_noisy - sig_phi;
%             Q_cal = [real(sig_dphi);imag(sig_dphi)];
%             tr = [real(yr);imag(yr)];
%             dphi = pinv(Q_cal) * tr;
%             phi_hat(ii+1) = phi_hat(ii) + dphi;
%             
%             arx_hat = exp(1j*(0:Nr-1)'*pi*sin(phi_hat(ii+1)))/sqrt(Nr);
%             
%             %---------------------------------------------
%             % Theta estimation using previous coeff.
%             %---------------------------------------------      
%             dTheta = exp(1j*pi*(0:Nt-1).'*sin(theta_hat(ii)))/sqrt(Nt).*(1j*pi*(0:Nt-1).'*cos(theta_hat(ii)));
% 
%             D_cal = alpha_hat(ii) * diag(W' * arx_hat * dTheta' * F);
%             sig_dtheta = kron(D_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             H_cal = alpha_hat(ii) * diag(W' * arx_hat * atx_hat' * F);
%             sig_theta = kron(H_cal,delay_mtx(:,maxindex)).*phase_error_refine;
%             
%             yr = sig_noisy - sig_theta;
%             Q_cal = [real(sig_dtheta);imag(sig_dtheta)];
%             tr = [real(yr);imag(yr)];
% 
%             dtheta = pinv(Q_cal) * tr;
%             theta_hat(ii+1) = theta_hat(ii) + dtheta;
% 
%         end
%         theta_last = theta_hat(ii-1); %theta_hat(1); %
%         phi_last = phi_hat(ii-1); %phi_hat(1);
%         eF_last = eF_hat(ii-1); %eF_hat(1);
        
        theta_az_last = theta_az_hat(1); %
        theta_el_last = theta_el_hat(1); %
        phi_az_last = phi_az_hat(1);
        eF_last = eF_hat(1);
        
        AOA_az_error_comp(MCidx,ss) = abs(phi_az_last - phi0_az(MCidx));
        AOD_az_error_comp(MCidx,ss) = abs(theta_az_last - theta0_az(MCidx));
        AOD_el_error_comp(MCidx,ss) = abs(theta_el_last - theta0_el(MCidx));
        
        CFO_error(MCidx,ss) = abs(eF_last- eF);
%         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);

        % RMSE evaluation from CRLB perspective
%         CRLB_CFO(ss,MCindex) = sqrt(temp(1,1))*(1/Ts/2/pi);
%         CRLB_theta(ss,MCindex) = sqrt(temp(2,2))*(1/pi*180);
%         CRLB_phi(ss,MCindex) = sqrt(temp(3,3))*(1/pi*180);
    end
end
%%
for ss=1:SNR_num
AOA_az_align_comp_mean(ss) = sum((AOA_az_error_comp(:,ss)/pi*180)<(105/Nr_az),1)/MCtimes;
AOD_az_align_comp_mean(ss) = sum((AOD_az_error_comp(:,ss)/pi*180)<(105/Nt_az),1)/MCtimes;
AOD_el_align_comp_mean(ss) = sum((AOD_el_error_comp(:,ss)/pi*180)<(105/Nt_el),1)/MCtimes;

Align_comp_mean(ss) = sum(((AOD_az_error_comp(:,ss)/pi*180)<(105/Nt_az)&...
                           (AOD_el_error_comp(:,ss)/pi*180)<(105/Nt_el)&...
                           (AOA_az_error_comp(:,ss)/pi*180)<(105/Nr_az)),1)/MCtimes;

end

figure
% plot(SNR_range,AOAalign_comp_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_comp_mean,'--','linewidth',2);hold on

plot(SNR_range,AOA_az_align_comp_mean,'-','linewidth',2);hold on
plot(SNR_range,AOD_az_align_comp_mean,'-','linewidth',2);hold on
plot(SNR_range,AOD_el_align_comp_mean,'-','linewidth',2);hold on
plot(SNR_range,Align_comp_mean,'-','linewidth',2);hold on

grid on
xlabel('SNR (dB)')
ylabel('Misalignment Rate')
legend('AoA (az)','AoD (az)','AOD (el)','Full Alignment')
%% RMSE evaluation
portion = 1; % for debug 
for ss=1:SNR_num
    delay_error = abs(delay_est(ss,:)-tau_samp).^2;
    RMSE_delay(ss) = sqrt(mean(get_portion(delay_error,portion)))*(Ts/(2*pi)/1e-9);
%     RMSE_delay(ss) = sqrt(mean((delay_error)))*(Ts/(2*pi)/1e-9);

    
    input_sig = abs(CFO_error(:,ss)).^2;
    RMSE_CFO(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/Ts/2/pi);
%     RMSE_CFO(ss) = sqrt(mean((input_sig)))*(1/Ts/2/pi);

    
    input_sig = abs(AOD_az_error_nocomp(:,ss)).^2;
    RMSE_theta(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
%     RMSE_theta(ss) = sqrt(mean((input_sig)))*(1/pi*180);

    
    input_sig = abs(AOA_error_nocomp(:,ss)).^2;
%     RMSE_phi(ss) = sqrt(mean((input_sig)))*(1/pi*180);
    RMSE_phi(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
end
figure
subplot(311)
semilogy(SNR_range,RMSE_delay);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of normalized delay [ns]')
title('delay est')

subplot(312)
semilogy(SNR_range,RMSE_CFO);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of CFO est [Hz]')
title('CFO est')

subplot(313)
semilogy(SNR_range,RMSE_theta);hold on
semilogy(SNR_range,RMSE_phi);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA/AoD [deg]')
title('angle est')

%%
for ss=1:length(SNR_range)
    for MCidx = 1:MCtimes
        AOA_error_nofine(MCidx,ss) = abs(bestAOA(MCidx,ss) - phi0_az(MCidx));
        AOD_error_nofine(MCidx,ss) = abs(bestAOD(MCidx,ss) - theta0_az(MCidx));
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

%% Plot CRLB of angle est./tracking
% figure
% subplot(211)
% semilogy(SNR_range,mean(CRLB_theta,2));hold on
% semilogy(SNR_range,mean(CRLB_phi,2));hold on
% grid on
% legend('Theta','Phi')
% xlabel('Point-to-Point SNR [dB]')
% ylabel('RMSE of AoA/AoD [deg]')
% 
% subplot(212)
% semilogy(SNR_range,mean(CRLB_CFO,2));hold on
% grid on
% title('CRLB of CFO')
% xlabel('Point-to-Point SNR [dB]')
% ylabel('RMSE of CFO [Hz]')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
%% best should be 3691
for mm=1:10
    tempchan(mm) = (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll));
end
[Measure_mat_new(1:10,3691) 22.6175*tempchan.']
