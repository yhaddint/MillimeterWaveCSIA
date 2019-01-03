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
MCtimes = 2e1;                             % Num of Monte Carlo Sim.

AOAspread2 = 0;                             % Intra-cluster AoA spread square 
AOAspread = 0;                              % Intra-cluster AoA spread RMS 
AODspread2 = 0;                             % Intra-cluster AoD spread RMS 
AODspread = 0;                              % Intra-cluster AoD spread RMS 
SNR_num = 4;                               % Num of Monte Carlo Sim.
SNR_range = linspace(-30,0,SNR_num);
BW = 57.6e6;                                % IA bandiwdth [Hz]
Ts = 1/BW;                                  % Sample duration
Nb = 512;                                   % Sample per SS burst
CFO_ppm = 1;                                % CFO in ppm
CFO = (fc/1e6*CFO_ppm);                     % CFO with unit Hz
eF = CFO*Ts*2*pi;                           % CFO normalized with Ts
CFO_samp = eF;                              % same thing
P = 128;                                    % Number of subcarrier for PSS
DFT = dftmtx(P);
to_est_CFO = 1;                             % Assuming perfect knowledge of CFO or not
max_ite_num = 1e3;                          % Number iteration in refinement steps
refine_CFO = 1;                             % Turn on refine CFO when coarse estimation of AoA/AoD (long long time!!!)

Nc = 4;                                     % maximum multipath delay in [samples]
STO_max = 500;                             % range of maximum integer offset

az_lim = pi/3;                              % Az. range limit (by default -60 to 60 deg)
el_lim = pi/6;                              % El. range limit (by default -30 to 30 deg)

CP = 8;                                     % cyclic prefix in [samples]
OFDM_sym_num = 4;                           % Num of OFDM in each burst
burst_N = P * OFDM_sym_num;                 % Num of sample in each SS burst

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


%-------- CS Dictionary generation -------------
cand_num_r_az = 2 * Nr_az+ 1;               % Grid size for AoA Azimuth (2 times Nr_az)
cand_num_r_el = 1;                          % Grid size for AoA Elevation (1 by default)

cand_num_t_az = 2 * Nt_az+ 1;               % Grid size for AoD Azimuth (2 times Nt_az)
cand_num_t_el = 2 * Nt_el+ 1;               % Grid size for AoD elevation (2 times Nt_el)

dict_num = cand_num_r_az * cand_num_r_el * cand_num_t_az * cand_num_t_el;

cand_y = zeros(M, dict_num);

OMP_grid_rx_az = linspace(-az_lim, az_lim, cand_num_r_az);
AOAstep_az = OMP_grid_rx_az(2) - OMP_grid_rx_az(1);

OMP_grid_tx_el = linspace(-el_lim, el_lim, cand_num_t_el);
AODstep_el = OMP_grid_tx_el(2) - OMP_grid_tx_el(1);

OMP_grid_tx_az = linspace(-az_lim, az_lim, cand_num_t_az);
AODstep_az = OMP_grid_tx_az(2) - OMP_grid_tx_az(1);

grid_ARV_r_az = exp(1j*(0:Nr_az-1)'*pi*sin(OMP_grid_rx_az));
grid_ARV_r = grid_ARV_r_az;

grid_ARV_t_az = exp(1j*(0:Nt_az-1)'*pi*sin(OMP_grid_tx_az));
grid_ARV_t_el = exp(1j*(0:Nt_el-1)'*pi*sin(OMP_grid_tx_el));

% Planar geometry is [ANT_row1.'; ANT_row2.' cdots, ANT_row_Nr_el.']
grid_ARV_t = kron(grid_ARV_t_el, grid_ARV_t_az);

% Scenario when assuming CFO is known (debug)
phase_error_mat = kron(exp(1j*eF*Nb*(0:M-1)),exp(1j*eF*(0:P-1)'));
phase_error = reshape(phase_error_mat,M*P,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       IA Sector Sounding BF and related
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of sector beams in steering (total M=64)
M_BS_burst_az = 8;
M_BS_burst_el = 2;
M_UE_burst_az = 4;

M_burst = [M_BS_burst_az*M_BS_burst_el, M_UE_burst_az];

% angle grid for sector beams; sector approach gives estimator from grid
BS_az_grid = (az_lim)*linspace(-1+(1/M_BS_burst_az),1-(1/M_BS_burst_az),M_BS_burst_az);
BS_el_grid = (el_lim)*linspace(-1+(1/M_BS_burst_el),1-(1/M_BS_burst_el),M_BS_burst_el);
UE_az_grid = (az_lim)*linspace(-1+(1/M_UE_burst_az),1-(1/M_UE_burst_az),M_UE_burst_az);

% codebook for sector IA approach
W_sec_mat = get_IA_BF(Nr_az, M_UE_burst_az, 'sector_FSM_KW'); % Rx BF in IA stage
F_sec_mat = get_IA_BF_3D(Nt_az, Nt_el,...
                     M_BS_burst_az, M_BS_burst_el,...
                     'sector_FSM_KW', az_lim, el_lim); % Tx beamformer in IA stage


% ---- PSS and signal length parameters --------
ZC_root = 29;                               % ZC root, a coprime number with ZC length
ZC_N = 127;                                 % ZC sequence length
seq = lteZadoffChuSeq(ZC_root,ZC_N);        % Generate ZC sequence

% Add CP in PSS
seq_1DC = ifft(seq)*sqrt(ZC_N+1);           % T-domain sig. for ZC corr in detction & STO estimation; DC subcarrier is not null
burst_sig = [seq_1DC(end-CP+1:end);...      % CP
             seq_1DC;...                    % ZC seq. in T-domain
             zeros(burst_N-(ZC_N+CP),1)];   % SSS/PBCH but modeled as zero here
Tx_sig_CP = repmat(burst_sig,M,1);          % Repetition for M burst
burst_length = length(burst_sig);           % Num. of samples in M ZC burst
Tx_sig_length = length(Tx_sig_CP);          % Num. of samples in M ZC burst (leave space for STO)
ZC_t_domain = conj(flipud(seq_1DC));        % ZC sequence used for correlation in Rx


% Debug flag used somewhere
debug_flag = 0;           

% For loop of Monte Carlo Sim. (iid. Realization of Channel, Noise, and Sounding BF)
for MCidx = 1:MCtimes
    
    % Print iteration num.
    clc; fprintf('Monte Carlo Run %d out of %d\n',MCidx, MCtimes);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Channel parameter and related (SV model)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % AoA of rays with disired seperation
    phi_az = zeros(path_num,1);
    phi0_az(MCidx) = (rand*90-45)/180*pi;%20/180*pi;%(rand*90-45)/180*pi;0.2618;%
    phi_az = phi0_az(MCidx) + randn(path_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta_az = zeros(path_num,1);
    theta_el = zeros(path_num,1);

    theta0_az(MCidx) = (rand*90-45)/180*pi;%20/180*pi;%; 0.1309;%
    theta0_el(MCidx) = (rand*30-15)/180*pi;%0;%%;0.2618;%

    theta_az = theta0_az(MCidx) + randn(path_num,1) * AODspread;
    theta_el = theta0_el(MCidx) + randn(path_num,1) * AOAspread;
    
    % Find closest AoA/AoD in grid (for debug)
    [~,row_true] = min(abs(OMP_grid_rx_az - phi0_az(MCidx)));
    [~,col_az_true] = min(abs(OMP_grid_tx_az - theta0_az(MCidx)));
    [~,col_el_true] = min(abs(OMP_grid_tx_el - theta0_el(MCidx)));
    
    % Find index in dictionary that correponding to closest AoA/AoD
    index_true = ((col_el_true-1) * cand_num_t_az + col_az_true - 1) * cand_num_r_az + row_true;

    % Rotate of ray
    tau = rand*(90e-9);
    pathdelay = [0];
    tau_samp(MCidx) = tau/Ts*2*pi;
    g_ray = exp(1j*rand*2*pi);
    
    % Pre-compute some vectors/matrices
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
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       IA PN Sounding BF and related
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);
    
    Measure_mat = kron(transpose(F)*conj(grid_ARV_t),W'*grid_ARV_r);
    
    % Only each pair of m-th column of F and W are useful
    select_row = zeros(1,M);
    for ii=1:M
        select_row(ii) = (ii-1) * M + ii;
    end
    Measure_mat_new = Measure_mat(select_row,:);
    
    % Precompute normalization
    for cc=1:cand_num_r_az*cand_num_t_az
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Waveform, Phase Error and related
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % weighted average
    for dd=1:dict_num
        index_new = [abs(Measure_mat_new(1:M-1,dd)),abs(Measure_mat_new(2:M,dd))];
        CFO_weight = min(index_new,[],2);
        CFO_weight_norm = CFO_weight/sum(CFO_weight);
        CFO_select(:,dd) = CFO_weight_norm>1/M;
    end

    % ------ Phase noise process -----------
    PN_seq = zeros(Nb * M, 1);
    PN_seq(1) = 1;
    for ll=1:length(PN_seq)
        PN_seq(ll+1) = PN_seq(ll).*exp(1j * randn * PN_sigma);
    end
    
    % ------ About CFO and its derivative --------
    qvec = exp(1j * (0:P-1)' * eF);
    
    % Received signals (When Ideal CP-Removal is Used)
    symb = [seq;1]; %exp(1j*rand(P,1)*2*pi);
    tau_num = 500;
    delay_cand = linspace(0,100,tau_num)*1e-9/Ts*2*pi;
    for tt=1:tau_num
        delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P).*symb);
    end
    
    % ------- Precompute for CS Search (When Ideal CP-Removal is Used) -------------
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
    
    % ------- Precompute for Sector Search (When Ideal CP-Removal is Used) -------------
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
    
    % ------- Precompute for Sector Search (When Nonideal CP-Removal is Used) -------------
    sig_rx_sec2 = zeros(P*M, 1);
    STO = 100; % timing offset as the sample number
    
    % received sample number 
    Rx_sig_length = burst_N * M + ZC_N - 1 + STO; % signal length after ZC correlation;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     IA Waveform (Sector Approach, when Nonideal CP-Removal is Used)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for path_index = 1:path_num
        % ------- Channel Generation --------
        H_chan = get_H_NB_3D(g_ray(path_index),...
                             phi_az(path_index),...
                             theta_az(path_index),...
                             theta_el(path_index),...
                             1,...                  % cluster number
                             1,...                  % ray number
                             Nt_az, Nt_el, Nr);     % Generate discrete time domain frequency-flat channel
        H_chan0 = H_chan./norm(H_chan,'fro')*sqrt(Nt*Nr/path_num); % H per multipath

        % ----- sector version received signal generation ------
        precoder_index_old = 0;
        combiner_index_old = 0;
        for nn=1:Tx_sig_length

            precoder_index = floor( (nn-1) / (burst_length*M_burst(2)) )+1;
            combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
            combiner_index = mod(combiner_index_raw-1,M_burst(2))+1;

            if (precoder_index ~= precoder_index_old) || (combiner_index ~= combiner_index_old)

                w_vec = W_sec_mat(:,combiner_index);
                v_vec = F_sec_mat(:,precoder_index);
                g_effective = (w_vec'*H_chan0*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
%             index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
%             Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
        Rx_sig0_sec(:,path_index) = g_save_debug.' .* Tx_sig_CP;
    end
    
    % ----- Summation over all delay tap for freq selective sim --------
    Rx_sig_sec = zeros(burst_length*M,1);
    for path_index = 1:path_num
        timewindow0 = (1+pathdelay(path_index)):burst_length*M;
        timewindow1 = 1:(burst_length*M-pathdelay(path_index));
        Rx_sig_sec(timewindow0,1) = Rx_sig_sec(timewindow0,1) + Rx_sig0_sec(timewindow1,path_index);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     IA waveform (PN approach, when Nonideal CP-Removal is Used)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nc_acc_mtx = toeplitz([1,zeros(1,burst_length-Nc)]',[ones(1,Nc),zeros(1,burst_length-Nc)]);
    for path_index = 1:path_num
        % ------- Channel Generation --------
        H_chan = get_H_NB_3D(g_ray(path_index),...
                             phi_az(path_index),...
                             theta_az(path_index),...
                             theta_el(path_index),...
                             1,...                  % cluster number
                             1,...                  % ray number
                             Nt_az, Nt_el, Nr);     % Generate discrete time domain frequency-flat channel
        H_chan0 = H_chan./norm(H_chan,'fro')*sqrt(Nt*Nr/path_num); % H per multipath

        % ----- sector version received signal generation ------
        precoder_index_old = 0;
        combiner_index_old = 0;
        for nn=1:Tx_sig_length

            precoder_index = floor( (nn-1) / burst_length )+1;
            combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
            combiner_index = mod(combiner_index_raw-1,M)+1;

            if (precoder_index ~= precoder_index_old) || (combiner_index ~= combiner_index_old)

                w_vec = W(:,combiner_index);
                v_vec = F(:,precoder_index);
                g_effective = (w_vec'*H_chan0*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
%             index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
%             Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
        Rx_sig0_PN(:,path_index) = g_save_debug.' .* Tx_sig_CP;
    end
    
    % ----- summation over all delay tap for freq selective sim --------
    Rx_sig_PN = zeros(burst_length*M,1);
    for path_index = 1:path_num
        timewindow0 = (1+pathdelay(path_index)):burst_length*M;
        timewindow1 = 1:(burst_length*M-pathdelay(path_index));
        Rx_sig_PN(timewindow0,1) = Rx_sig_PN(timewindow0,1) + Rx_sig0_PN(timewindow1,path_index);
    end
    
    % ------- AWGN -------
    noise_CP = (randn(Tx_sig_length,1)+1j*randn(Tx_sig_length,1))/sqrt(2);
    noise_at_STO = (randn(STO,1)+1j*randn(STO,1))/sqrt(2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     Detection and Timing Est. from Waveform (both approach)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ----- Initializations of vec, mat -----
%     Rx_sig = zeros(Tx_sig_length, 1); % received signal in t domain
    corr_out_H1 = zeros(Rx_sig_length + burst_N - ZC_N + 1, 1); % pad a zero at the end
    corr_out_H0 = zeros(Rx_sig_length + burst_N - ZC_N + 1,1); % pad a zero at the end
    
    % ----- For loop for SNR (Both BF Approach; In Detection) -----
    BFtype = 'PN';
    for ss = 1:SNR_num
        
        noise_pow = 10^(-SNR_range(ss)/10);
        awgn_CP = noise_CP * sqrt(noise_pow);
        
        % case 'sector'
        Rx_sig_H1_sec = Rx_sig_sec.*exp(1j*CFO_samp*(0:length(Rx_sig_sec)-1).') + awgn_CP ;
        Rx_sig_H0_sec = awgn_CP;
        
        % case 'PN'
        Rx_sig_H1 = Rx_sig_PN.*exp(1j*CFO_samp*(0:length(Rx_sig_PN)-1).') + awgn_CP ;
        Rx_sig_H0 = awgn_CP;


        % ------ T Domain ZC Correlation -------

        Rx_sig_H0_wSTO = [noise_at_STO * sqrt(noise_pow); Rx_sig_H0];
        Rx_sig_H1_wSTO = [noise_at_STO * sqrt(noise_pow); Rx_sig_H1];
        corr_out_H1_STO = abs(conv(ZC_t_domain,Rx_sig_H1_wSTO)/ZC_N).^2; % corr rx t-domain sig with ZC
        corr_out_H0_STO = abs(conv(ZC_t_domain,Rx_sig_H0_wSTO)/ZC_N).^2; % corr rx t-domain sig with ZC

        Rx_sig_H0_wSTO_sec = [noise_at_STO * sqrt(noise_pow); Rx_sig_H0_sec];
        Rx_sig_H1_wSTO_sec = [noise_at_STO * sqrt(noise_pow); Rx_sig_H1_sec];
        corr_out_H1_STO_sec = abs(conv(ZC_t_domain,Rx_sig_H1_wSTO_sec)/ZC_N).^2; % corr rx t-domain sig with ZC
        corr_out_H0_STO_sec = abs(conv(ZC_t_domain,Rx_sig_H0_wSTO_sec)/ZC_N).^2; % corr rx t-domain sig with ZC
        
        % ----- Multi-Peak Detection w. PN beam ---------
        % Practical scenario where peak location is unknown
        post_corr_ED_H1 = sum(reshape([corr_out_H1_STO(ZC_N+CP:end);...
            zeros(burst_length*(M+2)-length(corr_out_H1_STO(ZC_N+CP:end)),1)],...
            burst_length,M+2),2)/M;

        post_corr_ED_H0 = sum(reshape([corr_out_H0_STO(ZC_N+CP:end);...
            zeros(burst_length*(M+2)-length(corr_out_H0_STO(ZC_N+CP:end)),1)],...    
            burst_length,M+2),2)/M;

        ave_Nc_H0 = Nc_acc_mtx*post_corr_ED_H0;
        ave_Nc_H1 = Nc_acc_mtx*post_corr_ED_H1;
        [peak_pow_H1(ss) peakindex_H1(ss)] = max(ave_Nc_H1(1:STO_max));
        peak_pow_H0(ss) = max(ave_Nc_H0(1:STO_max));

        % ----- Single-Peak Detection w. directional beams---------
        % Practical scenario where peak location is unknown
        peak_pow_H1_sec(ss) = max(corr_out_H1_STO_sec);
        peak_pow_H0_sec(ss) = max(corr_out_H0_STO_sec);

        
        % SNR and Adding AWGN to Rx Signal (when Ideal CP-removal Sig. for ChEst)
%         sigman2 = 10^(-SNR_range(ss)/10);
%         sig_noisy = sig_rx + awgn * sqrt(sigman2); % Use 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     Sig. Rearrangement from Est STO (CP-Removal, PSS extraction)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sig_rx_from_dec = zeros(P*M,1);
        sig_rx_from_dec_sec = zeros(P*M,1);
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            index_from_rx_sig = peakindex_H1(ss) + CP + ((mm-1)*Nb:(mm-1)*Nb+P-1);
            sig_rx_from_dec(index) = Rx_sig_H1_wSTO(index_from_rx_sig);
            
            % there must be some way for timing est. of sec beam
            % but here I just use ideal STO as "optimistic" results 
            index_from_rx_sig_sec = (STO+1) + CP + ((mm-1)*Nb:(mm-1)*Nb+P-1);
            sig_rx_from_dec_sec(index) = Rx_sig_H1_wSTO_sec(index_from_rx_sig);
        end
        sig_noisy = sig_rx_from_dec;
        sig_noisy_sec = sig_rx_from_dec_sec;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     BF Training (Sector Approach)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         sig_sec_ave = mean(abs(reshape(sig_rx_sec + awgn * sqrt(sigman2), P, M)).^2,1);
        sig_sec_ave = mean(abs(reshape(sig_noisy_sec, P, M)).^2,1);
        [~,best_sec_idx] = max(sig_sec_ave);
        
        % From SS burst index to get sector index in UE, BS_az, and BS_el
        best_sec_BS_idx = floor((best_sec_idx-1)/M_UE_burst_az)+1;
        best_sec_UE_idx = best_sec_idx - (best_sec_BS_idx-1)*M_UE_burst_az;
        best_sec_BS_el_idx = floor((best_sec_BS_idx-1)/M_BS_burst_az)+1;
        best_sec_BS_az_idx = best_sec_BS_idx - (best_sec_BS_el_idx-1)*M_BS_burst_az;
        
        % From sector index to index in the angle estimation
        sec_UE_az_est = UE_az_grid(best_sec_UE_idx);
        sec_BS_az_est = BS_az_grid(best_sec_BS_az_idx);
        sec_BS_el_est = BS_el_grid(best_sec_BS_el_idx);
        
        % Error Evaluation
        AOA_az_error_sec(MCidx,ss) = abs(sec_UE_az_est - phi0_az(MCidx));
        AOD_az_error_sec(MCidx,ss) = abs(sec_BS_az_est - theta0_az(MCidx));
        AOD_el_error_sec(MCidx,ss) = abs(sec_BS_el_est - theta0_el(MCidx));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     BF Training (Compressive Approach)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -------- Delay Matching Pursuit (Ignore CFO) ------------
        sig_ave = mean(reshape(sig_noisy,P,M),2);
        for tt=1:tau_num
            score(tt) = abs(sum(sig_ave.*conj(delay_mtx(:,tt)))/norm(delay_mtx(:,tt))^2);
        end
        
        % Plug-in Matching Pursuit index for CFO and Delay est.
        [maxscore,maxindex] = max(score);
        maxindex_est(ss,MCidx) = maxindex;
        delay_est(ss,MCidx) = delay_cand(maxindex);
        
        % watch score to debug
%         figure;plot(delay_cand/2/pi*Ts/1e-9,score)
        
        % ---------------   Compensate Impact of Delay Taps   ----------------
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
        
        %------  Matching Pursuit LOOP (AoA/AoD pair in dictionary) ------
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
            
            % Estimated CFO or assuming perfect CFO info in debug
            if to_est_CFO
                
                sig_burst = sig_desymb.*conj(Measure_mat_new(:,dd))./abs(Measure_mat_new(:,dd));
                
                % adjust N/A numbers
                sig_burst(isnan(sig_burst))=0;
                
%                 %  ---------  Quasi-ML CFO est. Method 1 (M2 is faster)  ---------
%                 CFO_hat_new1 = zeros(M-1,1);
%                 for mm=1:M-1
%                     if CFO_select(mm,dd)
%                     CFO_hat_new1(mm) = angle(sig_burst(mm+1).*conj(sig_burst(mm)));
%                     end
%                 end
%                 CFO_est = sum(CFO_hat_new1)/sum(CFO_select(:,dd))/Nb;

%                 
                %  ---------  Quasi ML CFO est. Method 2  ---------
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

                    phase_est_test = exp(1j*CFO_range(zz)*Nb*(0:M-1).');
                    score_CFO(zz) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_est_test));

                end
                [~,CFO_index] = max(score_CFO);
                CFO_est = CFO_range(CFO_index);
                end

                % ----- true CFO (for debug) ---------
%                 CFO_est = eF;
                
                phase_error_mat = exp(1j*CFO_est*Nb*(0:M-1).');
                phase_error = phase_error_mat;
                CFO_final(dd) = CFO_est;
            end
            score_final(dd) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_error))...
                /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
        end
        
        % Get index in dictionary w. highest MP peak
        [~,bestindex_comp(MCidx)] = max(abs(score_final));
        
%         bestindex_comp(MCindex) = index_true; % debug. comment in main script

        % From index to on-grid angle estimator
        bestrow = floor((bestindex_comp(MCidx)-1)/cand_num_r_az)+1;
        bestcol = bestindex_comp(MCidx)-(bestrow-1)*cand_num_r_az;
        
        bestrow_el = floor((bestrow-1)/cand_num_t_az)+1;
        bestrow_az = bestrow - (bestrow_el-1)*cand_num_t_az;
        
        bestAOA_az(MCidx,ss) = (bestcol-1)*AOAstep_az-az_lim;
        
        bestAOD_az(MCidx,ss) = (bestrow_az-1)*AODstep_az-az_lim;
        bestAOD_el(MCidx,ss) = (bestrow_el-1)*AODstep_el-el_lim;
        
        % Plot score for debug
        if 0%debug_flag
            figure;
            plot(abs(score_final));hold on;
            plot(index_true,abs(score_final(index_true)),'o','linewidth',2,'markersize',10);
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

        % -------- Error Evaluation in beam training of CSIA ---------       
        AOA_az_error_comp(MCidx,ss) = abs(bestAOA_az(MCidx,ss) - phi0_az(MCidx));
        AOD_az_error_comp(MCidx,ss) = abs(bestAOD_az(MCidx,ss) - theta0_az(MCidx));
        AOD_el_error_comp(MCidx,ss) = abs(bestAOD_el(MCidx,ss) - theta0_el(MCidx));
        CFO_error(MCidx,ss) = abs((CFO_final(bestindex_comp(MCidx))*Nb+2*pi)/Nb - eF);

    end
end
%% Alignment evaluation

critical_width = 60; % 3dB beam width (half) 105/2 [deg]

for ss=1:SNR_num
AOA_az_align_comp_mean(ss) = sum((AOA_az_error_comp(:,ss)/pi*180)<(critical_width/Nr_az),1)/MCtimes;
AOD_az_align_comp_mean(ss) = sum((AOD_az_error_comp(:,ss)/pi*180)<(critical_width/Nt_az),1)/MCtimes;
AOD_el_align_comp_mean(ss) = sum((AOD_el_error_comp(:,ss)/pi*180)<(critical_width/Nt_el),1)/MCtimes;

Align_comp_mean(ss) = sum(((AOD_az_error_comp(:,ss)/pi*180)<(critical_width/Nt_az)&...
                           (AOD_el_error_comp(:,ss)/pi*180)<(critical_width/Nt_el)&...
                           (AOA_az_error_comp(:,ss)/pi*180)<(critical_width/Nr_az)),1)/MCtimes;
                       
Align_sec_mean(ss) = sum((( AOD_az_error_sec(:,ss)/pi*180)<(critical_width/Nt_az)&...
                           (AOD_el_error_sec(:,ss)/pi*180)<(critical_width/Nt_el)&...
                           (AOA_az_error_sec(:,ss)/pi*180)<(critical_width/Nr_az)),1)/MCtimes;

end

figure
% plot(SNR_range,AOAalign_comp_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_comp_mean,'--','linewidth',2);hold on

% plot(SNR_range,AOA_az_align_comp_mean,'-','linewidth',2);hold on
% plot(SNR_range,AOD_az_align_comp_mean,'-','linewidth',2);hold on
% plot(SNR_range,AOD_el_align_comp_mean,'-','linewidth',2);hold on
plot(SNR_range,Align_comp_mean,'-o','linewidth',2);hold on
plot(SNR_range,Align_sec_mean,'--x','linewidth',2);hold on



grid on
xlabel('SNR (dB)')
ylabel('Misalignment Rate')
% legend('AoA (az)','AoD (az)','AOD (el)','Full Alignment')
legend('CSIA, mmMAGIC UMi LOS','DIA, mmMAGIC UMi LOS')


