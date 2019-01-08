%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3);                                     %random seed
%-------------------------------------
% System Parameters
%-------------------------------------
path_num = 20;                              % Num of rays in a cluster

fc = 28e9;                                  % carrier freq [Hz]

Nr    = 16;                                 % Num of antenna in UE/Rx (total)
Nr_az = 4;                                  % Num of antenna in UE/Rx (azimuth)
Nr_el = 4;                                  % Num of antenna in UE/Rx (elevation)

Nt    = 64;                                 % Num of antenna in BS/Tx (total)                           
Nt_az = 4;                                  % Num of antenna in BS/Tx (azimuth) 
Nt_el = 16;                                 % Num of antenna in BS/Tx (elevation) 

M = 64;                                     % Length of SS bursts (IA OFDM symbols)
MCtimes = 2e2;                              % Num of Monte Carlo Sim.

BW_data = 400e6;                            % noise BW in evalution SNR of data stage
NF = 4;                                     % Receiver Noise Figure in [dB]
Tx_pow_dBW = 16;                            % Transmit power in [dBW], i.e., 46dBm - 30dB
SNR_num = 1;                                % Num of Monte Carlo Sim.
SNR_range = linspace(60,60,SNR_num);        % SNR range in evaluation
BW = 57.6e6;                                % IA bandiwdth [Hz]
Ts = 1/BW;                                  % Sample duration
Nb = 512;                                   % Sample per SS burst
CFO_ppm = 1;                                % CFO in ppm
CFO = (fc/1e6*CFO_ppm);                     % CFO with unit [Hz]
eF = CFO*Ts*2*pi;                           % CFO normalized with Ts
CFO_samp = eF;                              % same; phase rotate per sample due to CFO
P = 128;                                    % Number of subcarrier for PSS
DFT = dftmtx(P);
to_est_CFO = 1;                             % Assuming perfect knowledge of CFO or not
max_ite_num = 1e3;                          % Number iteration in refinement steps
refine_CFO = 1;                             % Turn on refine CFO when coarse estimation of AoA/AoD (long long time!!!)

Nc = 31;                                    % maximum multipath delay in [samples]
STO_max = 470;                              % range of maximum integer offset

CP = 32;                                    % cyclic prefix in [samples]
OFDM_sym_num = 4;                           % Num of OFDM in each burst
burst_N = P * OFDM_sym_num;                 % Num of sample in each SS burst
channel_model = 'QuaDRiGa';                 % use 'SV' or 'QuaDRiGa' model


% Genie detection threshold from running H0 scenario
% Here I've used predefined value to save time
% [CS_TH, sec_TH] = get_Genie_Detection_TH();
CS_TH = 1.5397e-13;
sec_TH = 6.9196e-14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SV model w. mmMAGIC UMi NLOS scene parameter 
%  For intra-cluster setting test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AOAspread_az = 22.1/180*pi;                 % Intra-cluster AoA spread square 
AOAspread_el = 10.0/180*pi;                 % Intra-cluster AoA spread RMS 
AODspread_az = 5.40/180*pi;                 % Intra-cluster AoD spread RMS 
AODspread_el = 0.30/180*pi;                 % Intra-cluster AoD spread RMS 
tauspread    = 10e-9;                       % intra-cluster delay spread RMS [second]

az_lim = pi/3;                              % Az. range limit (by default -60 to 60 deg)
el_lim = pi/6;                              % El. range limit (by default -30 to 30 deg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Phase Noise Specification
% FOM = L(f0,df) + 20log10(df/f0)+10log10(P_VCO)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VCO_FOM = -114 + 20*log10(1e6/fc)...
            + 10*log10(27);                 % phase noise FOM -114dBc@1MHz w. 27mW
P_VCO = 500;                                % Scaling PN var by changing VCO power [mW]
VCO_c = 10^(VCO_FOM/10)/P_VCO;              % parameter c in PN specs
PN_sigma2 = VCO_c*4*pi^2*fc^2*Ts;           % Time domain variance of PN Wiener process
PN_sigma = sqrt(PN_sigma2);                 % Weiner process RMS imcrement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           CS Dictionary generation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_num_r_az = 2 * Nr_az + 1;              % Grid size for AoA Azimuth (2 times Nr_az)
grid_num_r_el = 2 * Nr_el + 1;              % Grid size for AoA Elevation (2 times Nr_el)
grid_num_t_az = 2 * Nt_az + 1;              % Grid size for AoD Azimuth (2 times Nt_az)
grid_num_t_el = 2 * Nt_el + 1;              % Grid size for AoD elevation (2 times Nt_el)

dict_num = grid_num_r_az * grid_num_r_el * grid_num_t_az * grid_num_t_el;

cand_y = zeros(M, dict_num);

OMP_grid_rx_az = linspace(-az_lim, az_lim, grid_num_r_az);
AOAstep_az = OMP_grid_rx_az(2) - OMP_grid_rx_az(1);

OMP_grid_rx_el = linspace(-el_lim, el_lim, grid_num_r_el);
AOAstep_el = OMP_grid_rx_el(2) - OMP_grid_rx_el(1);

OMP_grid_tx_el = linspace(-el_lim, el_lim, grid_num_t_el);
AODstep_el = OMP_grid_tx_el(2) - OMP_grid_tx_el(1);

OMP_grid_tx_az = linspace(-az_lim, az_lim, grid_num_t_az);
AODstep_az = OMP_grid_tx_az(2) - OMP_grid_tx_az(1);

grid_ARV_r_az = exp(1j*(0:Nr_az-1)'*pi*sin(OMP_grid_rx_az));
grid_ARV_r_el = exp(1j*(0:Nr_el-1)'*pi*sin(OMP_grid_rx_el));

% Planar geometry is [ANT_row1.'; ANT_row2.' cdots, ANT_row_Nr_el.']
grid_ARV_r = kron(grid_ARV_r_el, grid_ARV_r_az);

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
M_BS_burst_az = 2;          % Num. of sector for BS az. plane
M_BS_burst_el = 8;          % Num. of sector for BS el. plane
M_UE_burst_az = 2;          % Num. of sector for UE az. plane
M_UE_burst_el = 2;          % Num. of sector for UE el. plane

% This has to be M = 64
M_burst = [M_BS_burst_az*M_BS_burst_el, M_UE_burst_az*M_UE_burst_el];

% angle grid for sector beams; sector approach gives estimator from grid
BS_az_grid_sec = (az_lim)*linspace(-1+(1/M_BS_burst_az),1-(1/M_BS_burst_az),M_BS_burst_az);
BS_az_grid_step = BS_az_grid_sec(2) - BS_az_grid_sec(1);

BS_el_grid_sec = (el_lim)*linspace(-1+(1/M_BS_burst_el),1-(1/M_BS_burst_el),M_BS_burst_el);
BS_el_grid_step = BS_el_grid_sec(2) - BS_el_grid_sec(1);

UE_az_grid_sec = (az_lim)*linspace(-1+(1/M_UE_burst_az),1-(1/M_UE_burst_az),M_UE_burst_az);
UE_az_grid_step = UE_az_grid_sec(2) - UE_az_grid_sec(1);

UE_el_grid_sec = (el_lim)*linspace(-1+(1/M_UE_burst_el),1-(1/M_UE_burst_el),M_UE_burst_el);
UE_el_grid_step = UE_el_grid_sec(2) - UE_el_grid_sec(1);


% codebook for sector IA approach
W_sec_mat = get_IA_BF_3D(Nr_az, Nr_el,...
                     M_UE_burst_az, M_UE_burst_el,...
                     'sector_FSM_KW', az_lim, el_lim); % Rx BF in IA stage
F_sec_mat = get_IA_BF_3D(Nt_az, Nt_el,...
                     M_BS_burst_az, M_BS_burst_el,...
                     'sector_FSM_KW', az_lim, el_lim); % Tx beamformer in IA stage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       PSS, CP and signal length parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Matrix Initialization in main loop 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug flag used somewhere
debug_flag = 0;           

data_mag_BF = zeros( MCtimes, SNR_num );
data_mag_sec = zeros( MCtimes, SNR_num );
data_mag_dir1 = zeros( MCtimes, SNR_num );
data_mag_dir = zeros( MCtimes, SNR_num );
data_mag_true = zeros( MCtimes, SNR_num );
data_mag_raw = zeros( MCtimes, SNR_num );
best_g_dir = zeros( MCtimes, SNR_num );

% For loop of Monte Carlo Sim. (iid. Realization of Channel, Noise, and Sounding BF)
for MCidx = 1:MCtimes
    
    % Print iteration num.
    clc; fprintf('Monte Carlo Run %d out of %d\n',MCidx, MCtimes);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %        QuaDRiGa Parameter and Setup
    %     It uses QuaDRiGa Tutorial 4.3 as baseline
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Include QuaDRiGa source folder
    addpath('C:\Users\Han\Documents\QuaDriGa_2017.08.01_v2.0.0-664\quadriga_src');

    % Fonts setting
    % set(0,'defaultTextFontSize', 18)                        % Default Font Size
    % set(0,'defaultAxesFontSize', 18)                        % Default Font Size
    % set(0,'defaultAxesFontName','Times')                    % Default Font Type
    % set(0,'defaultTextFontName','Times')                    % Default Font Type
    % set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
    % set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type
    % set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Default Paper Size

    s = qd_simulation_parameters;                           % New simulation parameters
    s.center_frequency = fc;                                % 28 GHz carrier frequency
    s.sample_density = 4;                                   % 4 samples per half-wavelength
    s.use_absolute_delays = 0;                              % Include delay of the LOS path
    s.show_progress_bars = 0;                               % Disable progress bars

    % Antenna array setting
    l = qd_layout( s );                                     % New QuaDRiGa layout
    l.tx_array = qd_arrayant('3gpp-mmw',Nt_el,Nt_az,fc,1,0,0.5,1,1);
    l.rx_array = qd_arrayant('3gpp-mmw',Nr_el,Nr_az,fc,1,0,0.5,1,1);
    % l.rx_array.rotate_pattern(180, 'z');
    % l.rx_array.visualize(1);pause(1)
    l.track = qd_track('linear',1,-pi);

    % BS and UE location 
    UE_dist = 70;                                            % BS/UE distance in x-axis
    UE_height = rand*20;                                     % Height of UE (up to 50m)
    BS_UE_pos_angle = (rand*90-45)/180*pi;
    x_UE = UE_dist * cos(BS_UE_pos_angle);
%     y_UE = UE_dist * sin(BS_UE_pos_angle);
    BS_height = 10;                                          % Height of BS
    y_UE = rand * (2 * UE_dist * tan(45/180*pi))...
        - (UE_dist * tan(45/180*pi));                        % BS/UE dist in y axis(ramdom)
    l.tx_position(:,1) = [0,0,BS_height].';                  % BS position
    l.rx_position(:,1) = [UE_dist,y_UE,UE_height].';         % UE position
    l.set_scenario('mmMAGIC_UMi_NLOS');                      % Set propagation scenario
%     l.visualize;                                           % Plot the layout

    % Call builder to generate channel parameter
    cb = l.init_builder;  
    
    % Some customized setting (for debuging)
    % cb.scenpar.PerClusterDS = 0;                            % Create new builder object
%     cb.scenpar.SF_sigma = 0;                                % 0 dB shadow fading
%     cb.scenpar.KF_mu = 0;                                   % 0 dB K-Factor
%     cb.scenpar.KF_sigma = 0;                                % No KF variation
%     cb.scenpar.SubpathMethod = 'mmMAGIC';
%     cb.plpar = [];                                          % Disable path loss model
    
    % Call builder for small scale fading parameters
    cb.gen_ssf_parameters;                                  % Generate large- and small-scale fading
    
    % Drifting model for channel dynamics (off)
    % s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
    % c = cb.get_channels;                                    % Generate channel coefficients
    % c.individual_delays = 0;                                % Remove per-antenna delays

    % Call builder to generate channel 
    s.use_spherical_waves = 0;                              % Disable drifting
    d = cb.get_channels;                                    % Generate channel coefficients

    % Rearrange since QuaDRiGa uses different az/el order
    % each column in MIMO is [el1; el2; etc]
    % while ours is [az1, az2, etc].'
    chan_coeff = zeros(Nr,Nt,cb.NumClusters,2);
    for c_idx = 1:cb.NumClusters
        chan_MIMO_old = squeeze(d.coeff(:,:,c_idx,1));

        row_order = repmat(1:Nr_el:Nr,1,Nr_el) + kron((0:Nr_el-1),ones(1,Nr_az));
        col_order = repmat(1:Nt_el:Nt,1,Nt_el) + kron((0:Nt_el-1),ones(1,Nt_az));

        chan_MIMO_new1 = chan_MIMO_old(:,col_order);
        chan_MIMO_new2 = chan_MIMO_new1(row_order,:);

        chan_coeff(:,:,c_idx,1) = chan_MIMO_new2;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %        Convert into Channel Taps
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_QDG_f = zeros(Nr,Nt,P);
    H_QDG_t = zeros(Nr,Nt,P);
    frac_delay = 10e-9;
    
    for pp=1:P
        for path_idx = 1:cb.NumClusters

        H_QDG_f(:,:,pp) = H_QDG_f(:,:,pp) +...
            exp(-1j*2*pi*(frac_delay + d.delay(path_idx,1))*(pp-1)/(Ts*P))...
            *squeeze(chan_coeff(:,:,path_idx,1));

        end
    end
    
    for nt_idx = 1:Nt
        for nr_idx = 1:Nr
            H_QDG_t(nr_idx,nt_idx,:) = ifft(H_QDG_f(nr_idx,nt_idx,:));
        end
    end
    
    % Alternatively using QuaDRiGa build-in function to convert to taps
    % Note this does not include any delay fractional to Ts from timing
    % offset except overwrite d.delay
%     h_tap = d.fr(BW,P,1);
%     for nt_idx = 1:Nt
%         for nr_idx = 1:Nr
%     H_QDG_f(nr_idx,nt_idx,:) = squeeze(h_tap(nr_idx,nt_idx,:,1));
%     H_QDG_t(nr_idx,nt_idx,:) = ifft(squeeze(H_QDG_f(nr_idx,nt_idx,:)),[],1);
%         end
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Channel parameter and related (SV model)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch channel_model
        case 'QuaDRiGa'
            % phi & theta are not functional unless in SV; Only for debug 
            phi0_az(MCidx) = cb.AoA(1)+pi;
            if phi0_az(MCidx)>pi
                phi0_az(MCidx) = - (2*pi - (cb.AoA(1) + pi));
            end
            
            phi0_el(MCidx) = cb.EoA(1);
            theta0_az(MCidx) = -cb.AoD(1);
            theta0_el(MCidx) = -cb.EoD(1);
        case 'SV'
            % cluster specific angles
            phi0_az(MCidx) = (rand*90-45)/180*pi;%0.5236;%20/180*pi;%(rand*90-45)/180*pi;
            phi0_el(MCidx) = (rand*30-15)/180*pi;%0.2618;%20/180*pi;%(rand*90-45)/180*pi;
            theta0_az(MCidx) = (rand*90-45)/180*pi;%0.1309;%20/180*pi;%; 0.1309;%
            theta0_el(MCidx) = (rand*30-15)/180*pi;%0.2618;%0;%%;
    end
    
    phi_az = zeros(path_num,1);
    phi_el = zeros(path_num,1);
    
    % AoA of rays 
    phi_az = phi0_az(MCidx) + laprnd(path_num, 1, 0, AOAspread_az);
    phi_el = phi0_el(MCidx) + laprnd(path_num, 1, 0, AOAspread_el);

    % AoD of rays 
    theta_az = zeros(path_num,1);
    theta_el = zeros(path_num,1);

    theta_az = theta0_az(MCidx) + laprnd(path_num, 1, 0, AODspread_az);
    theta_el = theta0_el(MCidx) + laprnd(path_num, 1, 0, AODspread_el);
    
    % delay spread; most simple case
%     ray_delay = 30e-9 + rand(path_num,1)*tauspread;
    
    % delay spread; mmMAGIC specific intra-cluster delay
    c_tau = 23; % intra-cluster delay spread with unit [ns]
    r_tau = 2.03; % parameter in mmMAGIC
    relative_delay_prime = -r_tau*c_tau*log(rand(path_num,1));
    relative_delay = sort(relative_delay_prime) - min(relative_delay_prime);
    delay_pow_scling = exp(-relative_delay*(r_tau-1)./(r_tau*c_tau));
    ray_delay = (relative_delay+10)*1e-9;
    
    % Find closest AoA/AoD in grid (for debug)
    [~,row_az_true] = min(abs(OMP_grid_rx_az - phi0_az(MCidx)));
    [~,row_el_true] = min(abs(OMP_grid_rx_el - phi0_el(MCidx)));
    [~,col_az_true] = min(abs(OMP_grid_tx_az - theta0_az(MCidx)));
    [~,col_el_true] = min(abs(OMP_grid_tx_el - theta0_el(MCidx)));
    
    % Find index in dictionary that correponding to closest AoA/AoD
    col_true = (col_el_true-1) * grid_num_t_az + col_az_true;
    row_true = (row_el_true-1) * grid_num_r_az + row_az_true;
    index_true = (col_true - 1) * (grid_num_r_az * grid_num_r_el) + row_true;

    % Rotate of ray
    tau = rand*(90e-9);
    pathdelay = zeros(path_num,1);
    tap_max = Nc;
    tapdelay = 0:(tap_max-1);
    tau_samp(MCidx) = tau/Ts*2*pi;
    
    % path gain scaled by intra-cluster delay
    g_ray = exp(1j*rand(path_num,1)*2*pi) .* delay_pow_scling;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Wideband channel in Tap domain
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch channel_model 
        case 'SV'
        H_chan_WB = get_H_WB_3D(g_ray.',...
                            ray_delay.',...
                            phi_az.',...
                            phi_el.',...
                            theta_az.',...
                            theta_el.',...
                            1,...                  % cluster number
                            path_num,...           % ray number
                            Nt_az, Nt_el,...
                            Nr_az, Nr_el,...
                            Ts,P);
                                          
        case 'QuaDRiGa'
        H_chan_WB = H_QDG_t * sqrt(10^(Tx_pow_dBW/10));
    end
    
    % Pre-compute some vectors/matrices (when ideal CP-Removal is used)
    for pathindex = 1:path_num
        
        % Spatial response and its derivative over phi
        arx_az(:,pathindex) = exp(1j * pi * (0:Nr_az-1)' * sin(phi_az(pathindex)))/sqrt(Nr_az);
        arx_el(:,pathindex) = exp(1j * pi * (0:Nr_el-1)' * sin(phi_el(pathindex)))/sqrt(Nr_el);
        arx(:,pathindex) = kron( arx_el(:,pathindex), arx_az(:,pathindex) );
        
        % Spatial response and its derivative over theta
        atx_az(:,pathindex) = exp(1j * pi * (0:Nt_az-1)' * sin(theta_az(pathindex)))/sqrt(Nt_az);
        atx_el(:,pathindex) = exp(1j * pi * (0:Nt_el-1)' * sin(theta_el(pathindex)))/sqrt(Nt_el);
        atx(:,pathindex) = kron( atx_el(:,pathindex), atx_az(:,pathindex) );
        
        % Delay response and its derivative over tau
        fvec(:,pathindex) = exp(-1j * (0:P-1)' * tau_samp(MCidx) / P);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       IA PN Sounding BF and sensing dictionary
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Quasi-omni beam from random steering mtx
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(M);
    
    % Sensing dictionary - post-BF response for each possible AoA/AoD pair
    % It has M^2 rows and only M of them are useful
    Measure_mat = kron(transpose(F)*conj(grid_ARV_t),W'*grid_ARV_r);
    
    % Only each pair of m-th column of F and W are useful
    select_row = zeros(1,M);
    for ii=1:M
        select_row(ii) = (ii-1) * M + ii;
    end
    Measure_mat_new = Measure_mat(select_row,:);
    
    % Precompute normalization
    for cc=1:grid_num_r_az*grid_num_t_az
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       CFO, Phase Error and related
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Weights in estimating CFO
    for dd=1:dict_num
        index_new = [abs(Measure_mat_new(1:M-1,dd)),abs(Measure_mat_new(2:M,dd))];
        CFO_weight = min(index_new,[],2);
        CFO_weight_norm = CFO_weight/sum(CFO_weight);
        CFO_select(:,dd) = CFO_weight_norm>1/M;
    end

    % ------ Phase noise process (When Ideal CP-Removal is Used) -----------
    PN_seq = zeros(Nb * M, 1);
    PN_seq(1) = 1;
    for ll=1:length(PN_seq)
        PN_seq(ll+1) = PN_seq(ll).*exp(1j * randn * PN_sigma);
    end
    
    % ------ About CFO and its derivative (When Ideal CP-Removal is Used) --------
    qvec = exp(1j * (0:P-1)' * eF);
    
    % Received signals (When Ideal CP-Removal is Used)
    symb = [seq;1]; %exp(1j*rand(P,1)*2*pi);
    tau_num = 500;
    delay_cand = linspace(0,250,tau_num)*1e-9/Ts*2*pi;
    for tt=1:tau_num
        delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P).*symb);
    end
    
    % ------- Precompute for CS Search (When Ideal CP-Removal is Used) -------------
    sig_rx = zeros(P*M, 1);
    for ll=1:path_num
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            indexPN = (mm-1)*Nb+1:(mm-1)*Nb+P;
            sig_rx(index) = sig_rx(index) + (g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec(:,ll).*symb)) * exp(1j*Nb*(mm-1)*eF)).*PN_seq(indexPN);
        end
    end
    
    awgn = (randn(P*M,1)+1j*randn(P*M,1))/sqrt(2);
    
    % ------- Precompute for Sector Search (When Ideal CP-Removal is Used) -------------
    sig_rx_sec = zeros(P*M, 1);
    for ll=1:path_num
        for mm=1:M
            
            mm_BS = floor( (mm-1)/(M_UE_burst_az * M_UE_burst_el) ) + 1;
            mm_UE = mm - (mm_BS - 1) * (M_UE_burst_az * M_UE_burst_el);
            
            index = (mm-1)*P+1:mm*P;
            indexPN = (mm-1)*Nb+1:(mm-1)*Nb+P;
            sig_rx_sec(index) = sig_rx_sec(index) + (g_ray(ll) * (W_sec_mat(:,mm_UE)'*arx(:,ll)) * conj(F_sec_mat(:,mm_BS)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec(:,ll).*symb)) * exp(1j*Nb*(mm-1)*eF)).*PN_seq(indexPN);
        end
    end
    
    % ------- Precompute for Sector Search (When Nonideal CP-Removal is Used) -------------
    sig_rx_sec2 = zeros(P*M, 1);
    STO = randi(50) + 50; % timing offset as the sample number
    
    % received sample number 
    Rx_sig_length = burst_N * M + ZC_N - 1 + STO; % signal length after ZC correlation;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     IA Waveform (Sector Approach, when Nonideal CP-Removal is Used)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for tap_index = 1:tap_max

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
                g_effective = (w_vec'*squeeze(H_chan_WB(:,:,tap_index))*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
%             index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
%             Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
        Rx_sig0_sec(:,tap_index) = g_save_debug.' .* Tx_sig_CP;
    end
    
    
    % ----- Summation over all delay tap for freq selective sim --------
    Rx_sig_sec = zeros(burst_length*M,1);
    for tap_index = 1:tap_max
        timewindow0 = (1+tapdelay(tap_index)):burst_length*M;
        timewindow1 = 1:(burst_length*M-tapdelay(tap_index));
        Rx_sig_sec(timewindow0,1) = Rx_sig_sec(timewindow0,1) + Rx_sig0_sec(timewindow1,tap_index);
    end
    
    % ------ Phase noise process (PN can be confusing) -----------
    PE_sec = zeros(burst_length*M-1, 1);
    PE_sec(1) = 1;
    for ll=1:length(PE_sec)
        PE_sec(ll+1) = PE_sec(ll).*exp(1j * randn * PN_sigma);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     IA waveform (PN approach, when Nonideal CP-Removal is Used)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % mtx to take moving average over Nc samples in Energy Detection
    Nc_acc_mtx = toeplitz([1,zeros(1,burst_length-Nc)]',[ones(1,Nc),zeros(1,burst_length-Nc)]);
    
    for tap_index = 1:tap_max
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
                g_effective = (w_vec'*squeeze(H_chan_WB(:,:,tap_index))*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
%             index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
%             Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
        Rx_sig0_PN(:,tap_index) = g_save_debug.' .* Tx_sig_CP;
    end
    
    % ----- summation over all delay tap for freq selective sim --------
    Rx_sig_PN = zeros(burst_length*M,1);
    for tap_index = 1:tap_max
        timewindow0 = (1+tapdelay(tap_index)):burst_length*M;
        timewindow1 = 1:(burst_length*M-tapdelay(tap_index));
        Rx_sig_PN(timewindow0,1) = Rx_sig_PN(timewindow0,1) + Rx_sig0_PN(timewindow1,tap_index);
    end
    
    % ------ Phase noise process (PN can be confusing) -----------
    PE_CS = zeros(burst_length*M-1, 1);
    PE_CS(1) = 1;
    for ll=1:length(PE_CS)
        PE_CS(ll+1) = PE_CS(ll).*exp(1j * randn * PN_sigma);
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
        
        % noise power for each SNR
%         noise_pow = 10^(-SNR_range(ss)/10);

        % Absolute noise power (100MHz equivalent BW)
        % Also everything should be in dBW!!! subtract 30dB
        noise_pow = 10^( (-174 + 10*log10(BW) + NF - 30)/10 );
        awgn_CP = noise_CP * sqrt(noise_pow);
        
        % case 'sector beamformer'
        Rx_sig_H1_sec = Rx_sig_sec.*...
            PE_sec .* exp(1j*CFO_samp*(0:length(Rx_sig_sec)-1).') + awgn_CP ;
        Rx_sig_H0_sec = awgn_CP;
        
        % case 'PN beamformer'
        Rx_sig_H1 = Rx_sig_PN.*PE_CS.*exp(1j*CFO_samp*(0:length(Rx_sig_PN)-1).') + awgn_CP ;
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
        
        % Moving average of Nc samples to include energy of all paths
        ave_Nc_H0 = Nc_acc_mtx * post_corr_ED_H0;
        ave_Nc_H1 = Nc_acc_mtx * post_corr_ED_H1;
        
        [peak_pow_H1(MCidx) peakindex_H1(ss)] = max(ave_Nc_H1(1:STO_max));
        peak_pow_H0(MCidx) = max(ave_Nc_H0(1:STO_max));

        % ----- Single-Peak Detection w. directional beams---------
        % Practical scenario where peak location is unknown
        peak_pow_H1_sec(MCidx) = max(corr_out_H1_STO_sec);
        peak_pow_H0_sec(MCidx) = max(corr_out_H0_STO_sec);
        
        % decision
        decision_CS(MCidx) = peak_pow_H1(MCidx) > CS_TH;
        decision_sec(MCidx) = peak_pow_H1_sec(MCidx) > sec_TH;

        
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
            
            % This line give ideal CP removal
%             index_from_rx_sig = (STO+1) + CP + ((mm-1)*Nb:(mm-1)*Nb+P-1);
            
            % This line used estimated STO to remove CP
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
        %     Coarse BF Training (Sector Approach)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         sig_sec_ave = mean(abs(reshape(sig_rx_sec + awgn * sqrt(sigman2), P, M)).^2,1);
        sig_sec_ave = mean(abs(reshape(sig_noisy_sec, P, M)).^2,1);
        [~,best_sec_idx] = max(sig_sec_ave);
        
        % From SS burst index to get sector index in UE, BS_az, and BS_el
        best_sec_BS_idx = floor((best_sec_idx-1)/(M_UE_burst_az*M_UE_burst_el))+1;
        best_sec_UE_idx = best_sec_idx - (best_sec_BS_idx-1)*(M_UE_burst_az*M_UE_burst_el);
        
        best_sec_UE_el_idx = floor((best_sec_UE_idx-1)/M_UE_burst_az)+1;
        best_sec_UE_az_idx = best_sec_UE_idx - (best_sec_UE_el_idx-1)*M_UE_burst_az;        
        best_sec_BS_el_idx = floor((best_sec_BS_idx-1)/M_BS_burst_az)+1;
        best_sec_BS_az_idx = best_sec_BS_idx - (best_sec_BS_el_idx-1)*M_BS_burst_az;
        
        % From sector index to index in the angle estimation
        sec_UE_az_est(MCidx,ss) = UE_az_grid_sec(best_sec_UE_az_idx);
        sec_UE_el_est(MCidx,ss) = UE_el_grid_sec(best_sec_UE_el_idx);
        sec_BS_az_est(MCidx,ss) = BS_az_grid_sec(best_sec_BS_az_idx);
        sec_BS_el_est(MCidx,ss) = BS_el_grid_sec(best_sec_BS_el_idx);
        
        % Error Evaluation
        AOA_az_error_sec(MCidx,ss) = abs(sec_UE_az_est(MCidx,ss) - phi0_az(MCidx));
        AOA_el_error_sec(MCidx,ss) = abs(sec_UE_el_est(MCidx,ss) - phi0_el(MCidx));
        AOD_az_error_sec(MCidx,ss) = abs(sec_BS_az_est(MCidx,ss) - theta0_az(MCidx));
        AOD_el_error_sec(MCidx,ss) = abs(sec_BS_el_est(MCidx,ss) - theta0_el(MCidx));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     #1 Stage Fine BF Training w. CSI-RS (Sector Approach)
        %
        %   CSI-RS frame is not sim-ed since it is not emphasis of paper
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fine_idx = 1;
        BS_az_fine1_num = 3;
        BS_el_fine1_num = 3;
        BS_fine1_num = BS_az_fine1_num * BS_el_fine1_num;
        
        % This number is picked in arbitary way; Can definitely be better
        BS_el_dir1_step = BS_el_grid_step/2;
        BS_az_dir1_step = BS_az_grid_step/2;
        
        UE_az_fine1_num = 3;
        UE_el_fine1_num = 3;
        UE_fine1_num = UE_az_fine1_num * UE_el_fine1_num;
        
        % This number is picked in arbitary way; Can definitely be better
        UE_el_dir1_step = UE_el_grid_step/2;
        UE_az_dir1_step = UE_az_grid_step/2;
        
        RSSI1 = zeros(UE_fine1_num * BS_fine1_num,1);

        for BS_el_fine1_idx = 1:BS_el_fine1_num
            BS_dir1_el = sec_BS_el_est(MCidx,ss) + (BS_el_fine1_idx-1) * BS_el_dir1_step;
            BS_dir1_el_vec = exp(1j*pi*(0:(Nt_el-1)).'*sin(BS_dir1_el))/sqrt(Nt_el);

            for BS_az_fine1_idx = 1:BS_az_fine1_num
                BS_dir1_az = sec_BS_az_est(MCidx,ss) + (BS_az_fine1_idx-1) * BS_az_dir1_step;
                BS_dir1_az_vec = exp(1j*pi*(0:(Nt_az-1)).'*sin(BS_dir1_az))/sqrt(Nt_az);
                
                for UE_el_fine1_idx = 1:UE_el_fine1_num
                    UE_dir1_el = sec_UE_el_est(MCidx,ss) + (UE_el_fine1_idx-1) * UE_el_dir1_step;
                    UE_dir1_el_vec = exp(1j*pi*(0:(Nr_el-1)).'*sin(UE_dir1_el))/sqrt(Nr_el);

                    for UE_az_fine1_idx = 1:UE_az_fine1_num
                        UE_dir1_az = sec_UE_az_est(MCidx,ss) + (UE_az_fine1_idx-1) * UE_az_dir1_step;
                        UE_dir1_az_vec = exp(1j*pi*(0:(Nr_az-1)).'*sin(UE_dir1_az))/sqrt(Nr_az); 
                
                        BS_dir1_vec = kron( BS_dir1_el_vec, BS_dir1_az_vec );
                        UE_dir1_vec = kron( UE_dir1_el_vec, UE_dir1_az_vec );
                        
                        % Evaluating RSSI
                        switch channel_model
                            case 'QuaDRiGa'
                                for ll = 1:path_num
                                    RSSI1(fine_idx) = RSSI1(fine_idx) + ...
                                        abs(UE_dir1_vec'*squeeze(chan_coeff(:,:,ll,1))*BS_dir1_vec)^2;
                                end
                            case 'SV'
                                for ll = 1:path_num
                                    RSSI1(fine_idx) = RSSI1(fine_idx) + g_ray(ll) * (UE_dir1_vec'*arx(:,ll)) ...
                                    * conj(BS_dir1_vec'*atx(:,ll));
                                end
                        end
                        fine_idx = fine_idx + 1;
                        
                    end
                end
            end
        end
        
        [best_g_dir1(MCidx,ss), best_fine1_idx] = max(abs(RSSI1));
        
        best_BS_dir1_idx = floor((best_fine1_idx-1)/UE_fine1_num) + 1;
        best_UE_dir1_idx = best_fine1_idx - (best_BS_dir1_idx-1) * UE_fine1_num;
        
        best_BS_dir1_el_idx = floor((best_BS_dir1_idx-1)/BS_az_fine1_num) + 1;
        best_BS_dir1_az_idx = best_BS_dir1_idx - (best_BS_dir1_el_idx-1) * BS_az_fine1_num;
        
        best_UE_dir1_el_idx = floor((best_UE_dir1_idx-1)/UE_az_fine1_num) + 1;
        best_UE_dir1_az_idx = best_UE_dir1_idx - (best_UE_dir1_el_idx-1) * UE_az_fine1_num;
        
        best_UE_dir1_el(MCidx,ss) = sec_UE_el_est(MCidx,ss) +...
                                    (best_UE_dir1_el_idx-1) * UE_el_dir1_step;
        best_UE_dir1_az(MCidx,ss) = sec_UE_az_est(MCidx,ss) +...
                                    (best_UE_dir1_az_idx-1) * UE_az_dir1_step;
        best_BS_dir1_el(MCidx,ss) = sec_BS_el_est(MCidx,ss) +...
                                    (best_BS_dir1_el_idx-1) * BS_el_dir1_step;
        best_BS_dir1_az(MCidx,ss) = sec_BS_az_est(MCidx,ss) +...
                                    (best_BS_dir1_az_idx-1) * BS_az_dir1_step;


        % Error Evaluation
        AOA_az_error_dir1(MCidx,ss) = abs(best_UE_dir1_az(MCidx,ss) - phi0_az(MCidx));
        AOA_el_error_dir1(MCidx,ss) = abs(best_UE_dir1_el(MCidx,ss) - phi0_el(MCidx));
        AOD_az_error_dir1(MCidx,ss) = abs(best_BS_dir1_az(MCidx,ss) - theta0_az(MCidx));
        AOD_el_error_dir1(MCidx,ss) = abs(best_BS_dir1_el(MCidx,ss) - theta0_el(MCidx));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     #2 Stage Fine BF Training w. CSI-RS (Sector Approach)
        %
        %   CSI-RS frame is not sim-ed since it is not emphasis of paper
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fine_idx = 1;
        BS_az_fine_num = 5;
        BS_el_fine_num = 5;
        BS_fine_num = BS_az_fine_num * BS_el_fine_num;
        
        % This number is picked in arbitary way; Can definitely be better
        BS_el_dir_step = BS_el_grid_step/4;
        BS_az_dir_step = BS_az_grid_step/4;
        
        UE_az_fine_num = 5;
        UE_el_fine_num = 5;
        UE_fine_num = UE_az_fine_num * UE_el_fine_num;
        
        % This number is picked in arbitary way; Can definitely be better
        UE_el_dir_step = UE_el_grid_step/4;
        UE_az_dir_step = UE_az_grid_step/4;
        
        RSSI = zeros(UE_fine_num * BS_fine_num,1);

        for BS_el_fine_idx = 1:BS_el_fine_num
            BS_dir_el = sec_BS_el_est(MCidx,ss) + (BS_el_fine_idx-3) * BS_el_dir_step;
            BS_dir_el_vec = exp(1j*pi*(0:(Nt_el-1)).'*sin(BS_dir_el))/sqrt(Nt_el);

            for BS_az_fine_idx = 1:BS_az_fine_num
                BS_dir_az = sec_BS_az_est(MCidx,ss) + (BS_az_fine_idx-3) * BS_az_dir_step;
                BS_dir_az_vec = exp(1j*pi*(0:(Nt_az-1)).'*sin(BS_dir_az))/sqrt(Nt_az);
                
                for UE_el_fine_idx = 1:UE_el_fine_num
                    UE_dir_el = sec_UE_el_est(MCidx,ss) + (UE_el_fine_idx-3) * UE_el_dir_step;
                    UE_dir_el_vec = exp(1j*pi*(0:(Nr_el-1)).'*sin(UE_dir_el))/sqrt(Nr_el);

                    for UE_az_fine_idx = 1:UE_az_fine_num
                        UE_dir_az = sec_UE_az_est(MCidx,ss) + (UE_az_fine_idx-3) * UE_az_dir_step;
                        UE_dir_az_vec = exp(1j*pi*(0:(Nr_az-1)).'*sin(UE_dir_az))/sqrt(Nr_az); 
                
                        BS_dir_vec = kron( BS_dir_el_vec, BS_dir_az_vec );
                        UE_dir_vec = kron( UE_dir_el_vec, UE_dir_az_vec );
                        
                        % Evaluating RSSI
                        switch channel_model
                            case 'QuaDRiGa'
                                for ll = 1:path_num
                                    RSSI(fine_idx) = RSSI(fine_idx) + ...
                                        abs(UE_dir_vec'*squeeze(chan_coeff(:,:,ll,1))*BS_dir_vec)^2;
                                end
                            case 'SV'
                                for ll = 1:path_num
                                    RSSI(fine_idx) = RSSI(fine_idx) + g_ray(ll) * (UE_dir_vec'*arx(:,ll)) ...
                                    * conj(BS_dir_vec'*atx(:,ll));
                                end
                        end
                        fine_idx = fine_idx + 1;
                        
                    end
                end
            end
        end
        
        [best_g_dir(MCidx,ss), best_fine_idx] = max(abs(RSSI));
        
        best_BS_dir_idx = floor((best_fine_idx-1)/UE_fine_num) + 1;
        best_UE_dir_idx = best_fine_idx - (best_BS_dir_idx-1) * UE_fine_num;
        
        best_BS_dir_el_idx = floor((best_BS_dir_idx-1)/BS_az_fine_num) + 1;
        best_BS_dir_az_idx = best_BS_dir_idx - (best_BS_dir_el_idx-1) * BS_az_fine_num;
        
        best_UE_dir_el_idx = floor((best_UE_dir_idx-1)/UE_az_fine_num) + 1;
        best_UE_dir_az_idx = best_UE_dir_idx - (best_UE_dir_el_idx-1) * UE_az_fine_num;
        
        best_UE_dir_el(MCidx,ss) = sec_UE_el_est(MCidx,ss) +...
                                    (best_UE_dir_el_idx-3) * UE_el_dir_step;
        best_UE_dir_az(MCidx,ss) = sec_UE_az_est(MCidx,ss) +...
                                    (best_UE_dir_az_idx-3) * UE_az_dir_step;
        best_BS_dir_el(MCidx,ss) = sec_BS_el_est(MCidx,ss) +...
                                    (best_BS_dir_el_idx-3) * BS_el_dir_step;
        best_BS_dir_az(MCidx,ss) = sec_BS_az_est(MCidx,ss) +...
                                    (best_BS_dir_az_idx-3) * BS_az_dir_step;


        % Error Evaluation
        AOA_az_error_dir(MCidx,ss) = abs(best_UE_dir_az(MCidx,ss) - phi0_az(MCidx));
        AOA_el_error_dir(MCidx,ss) = abs(best_UE_dir_el(MCidx,ss) - phi0_el(MCidx));
        AOD_az_error_dir(MCidx,ss) = abs(best_BS_dir_az(MCidx,ss) - theta0_az(MCidx));
        AOD_el_error_dir(MCidx,ss) = abs(best_BS_dir_el(MCidx,ss) - theta0_el(MCidx));
        
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
        bestrow = floor((bestindex_comp(MCidx)-1)/(grid_num_r_az*grid_num_r_el))+1;
        bestcol = bestindex_comp(MCidx)-(bestrow-1)*(grid_num_r_az*grid_num_r_el);
        
        bestcol_el = floor((bestcol-1)/grid_num_r_az)+1;
        bestcol_az = bestcol - (bestcol_el-1)*grid_num_r_az;
        
        bestrow_el = floor((bestrow-1)/grid_num_t_az)+1;
        bestrow_az = bestrow - (bestrow_el-1)*grid_num_t_az;
        
        bestAOA_az(MCidx,ss) = (bestcol_az-1)*AOAstep_az - az_lim;
        bestAOA_el(MCidx,ss) = (bestcol_el-1)*AOAstep_el - el_lim;
        
        bestAOD_az(MCidx,ss) = (bestrow_az-1)*AODstep_az - az_lim;
        bestAOD_el(MCidx,ss) = (bestrow_el-1)*AODstep_el - el_lim;
        
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
        AOA_el_error_comp(MCidx,ss) = abs(bestAOA_el(MCidx,ss) - phi0_el(MCidx));

        AOD_az_error_comp(MCidx,ss) = abs(bestAOD_az(MCidx,ss) - theta0_az(MCidx));
        AOD_el_error_comp(MCidx,ss) = abs(bestAOD_el(MCidx,ss) - theta0_el(MCidx));
        
        CFO_error(MCidx,ss) = abs((CFO_final(bestindex_comp(MCidx))*Nb+2*pi)/Nb - eF);
        
        % BF gain analysis (CS approach)
        w_az_data = exp(1j * pi * (0:Nr_az-1)' * sin(bestAOA_az(MCidx,ss)))/sqrt(Nr_az);
        w_el_data = exp(1j * pi * (0:Nr_el-1)' * sin(bestAOA_el(MCidx,ss)))/sqrt(Nr_el);
        w_data = kron( w_el_data, w_az_data ); % Post-training Rx Beamformer 

        v_az_data = exp(1j * pi * (0:Nt_az-1)' * sin(bestAOD_az(MCidx,ss)))/sqrt(Nt_az);
        v_el_data = exp(1j * pi * (0:Nt_el-1)' * sin(bestAOD_el(MCidx,ss)))/sqrt(Nt_el);
        v_data = kron( v_el_data, v_az_data ); % Post-training Tx Beamformer 
        
        % BF gain analysis (sec approach)
        w_az_data_sec = exp(1j * pi * (0:Nr_az-1)' * sin(sec_UE_az_est(MCidx,ss)))/sqrt(Nr_az);
        w_el_data_sec = exp(1j * pi * (0:Nr_el-1)' * sin(sec_UE_el_est(MCidx,ss)))/sqrt(Nr_el);
        w_data_sec = kron( w_el_data_sec, w_az_data_sec ); % Post-training Rx Beamformer     
        
        v_az_data_sec = exp(1j * pi * (0:Nt_az-1)' * sin(sec_BS_az_est(MCidx,ss)))/sqrt(Nt_az);
        v_el_data_sec = exp(1j * pi * (0:Nt_el-1)' * sin(sec_BS_el_est(MCidx,ss)))/sqrt(Nt_el);
        v_data_sec = kron( v_el_data_sec, v_az_data_sec ); % Post-training Tx Beamformer 
        
        
        % BF gain analysis (sec w. CSI-RS #1 stage)
        w_az_data_dir1 = exp(1j * pi * (0:Nr_az-1)' * sin(best_UE_dir1_az(MCidx,ss)))/sqrt(Nr_az);
        w_el_data_dir1 = exp(1j * pi * (0:Nr_el-1)' * sin(best_UE_dir1_el(MCidx,ss)))/sqrt(Nr_el);
        w_data_dir1 = kron( w_el_data_dir1, w_az_data_dir1 );
        
        v_az_data_dir1 = exp(1j * pi * (0:Nt_az-1)' * sin(best_BS_dir1_az(MCidx,ss)))/sqrt(Nt_az);
        v_el_data_dir1 = exp(1j * pi * (0:Nt_el-1)' * sin(best_BS_dir1_el(MCidx,ss)))/sqrt(Nt_el);
        v_data_dir1 = kron( v_el_data_dir1, v_az_data_dir1 );
        
        % BF gain analysis (sec w. CSI-RS #2 stage)
        w_az_data_dir = exp(1j * pi * (0:Nr_az-1)' * sin(best_UE_dir_az(MCidx,ss)))/sqrt(Nr_az);
        w_el_data_dir = exp(1j * pi * (0:Nr_el-1)' * sin(best_UE_dir_el(MCidx,ss)))/sqrt(Nr_el);
        w_data_dir = kron( w_el_data_dir, w_az_data_dir );
        
        v_az_data_dir = exp(1j * pi * (0:Nt_az-1)' * sin(best_BS_dir_az(MCidx,ss)))/sqrt(Nt_az);
        v_el_data_dir = exp(1j * pi * (0:Nt_el-1)' * sin(best_BS_dir_el(MCidx,ss)))/sqrt(Nt_el);
        v_data_dir = kron( v_el_data_dir, v_az_data_dir );
        
        
        % BF gain analysis (true cluster-specific steering)
        w_az_data_true = exp(1j * pi * (0:Nr_az-1)' * sin(phi0_az(MCidx)))/sqrt(Nr_az);
        w_el_data_true = exp(1j * pi * (0:Nr_el-1)' * sin(phi0_el(MCidx)))/sqrt(Nr_el);
        w_data_true = kron( w_el_data_true, w_az_data_true );
        
        v_az_data_true = exp(1j * pi * (0:Nt_az-1)' * sin(theta0_az(MCidx)))/sqrt(Nt_az);
        v_el_data_true = exp(1j * pi * (0:Nt_el-1)' * sin(theta0_el(MCidx)))/sqrt(Nt_el);
        v_data_true = kron( v_el_data_true, v_az_data_true );
        
        % Gain Evaluation
        switch channel_model
            case 'SV'
                for ll=1:path_num
                    data_mag_BF(MCidx,ss) = data_mag_BF(MCidx,ss) +...
                        g_ray(ll) * (w_data'*arx(:,ll)) * conj(v_data'*atx(:,ll));
                    data_mag_sec(MCidx,ss) = data_mag_sec(MCidx,ss) +...
                        g_ray(ll) * (w_data_sec'*arx(:,ll)) * conj(v_data_sec'*atx(:,ll));
                    data_mag_dir1(MCidx,ss) = data_mag_dir1(MCidx,ss) +...
                        g_ray(ll) * (w_data_dir1'*arx(:,ll)) * conj(v_data_dir1'*atx(:,ll));
                    data_mag_dir(MCidx,ss) = data_mag_dir(MCidx,ss) +...
                        g_ray(ll) * (w_data_dir'*arx(:,ll)) * conj(v_data_dir'*atx(:,ll));
                    data_mag_true(MCidx,ss) = data_mag_true(MCidx,ss) +...
                        g_ray(ll) * (w_data_true'*arx(:,ll)) * conj(v_data_true'*atx(:,ll));
                    data_mag_raw(MCidx,ss) = data_mag_raw(MCidx,ss) + g_ray(ll);
                end
            case 'QuaDRiGa'
                for ll=1:cb.NumClusters
                    data_mag_BF(MCidx,ss) = data_mag_BF(MCidx,ss) +...
                        abs(w_data'*squeeze(chan_coeff(:,:,ll,1))*v_data)^2;
                    data_mag_sec(MCidx,ss) = data_mag_sec(MCidx,ss) +...
                        abs(w_data_sec'*squeeze(chan_coeff(:,:,ll,1))*v_data_sec)^2;
                    data_mag_dir1(MCidx,ss) = data_mag_dir1(MCidx,ss) +...
                        abs(w_data_dir1'*squeeze(chan_coeff(:,:,ll,1))*v_data_dir1)^2;
                    data_mag_dir(MCidx,ss) = data_mag_dir(MCidx,ss) +...
                        abs(w_data_dir'*squeeze(chan_coeff(:,:,ll,1))*v_data_dir)^2;
                    data_mag_true(MCidx,ss) = data_mag_true(MCidx,ss) +...
                        abs(w_data_true'*squeeze(chan_coeff(:,:,ll,1))*v_data_true)^2;
                    data_mag_raw(MCidx,ss) = data_mag_raw(MCidx,ss) +...
                        norm(squeeze(chan_coeff(:,:,ll,1)),'fro')^2/(Nt*Nr);
                end
                
        end % end of channel model switch
        
    end % end of SNR sweep
end % end of Monte Carlo loop
%% Alignment evaluation

critical_width = 55; % 3dB beam width 105/N [deg]

for ss=1:SNR_num
% AOA_az_align_comp_mean(ss) = sum((AOA_az_error_comp(:,ss)/pi*180)<(critical_width/Nr_az),1)/MCtimes;
% AOA_el_align_comp_mean(ss) = sum((AOA_el_error_comp(:,ss)/pi*180)<(critical_width/Nr_el),1)/MCtimes;
% AOD_az_align_comp_mean(ss) = sum((AOD_az_error_comp(:,ss)/pi*180)<(critical_width/Nt_az),1)/MCtimes;
% AOD_el_align_comp_mean(ss) = sum((AOD_el_error_comp(:,ss)/pi*180)<(critical_width/Nt_el),1)/MCtimes;

AOA_az_align_dir_mean(ss) = sum((AOA_az_error_dir(:,ss)/pi*180)<(critical_width/Nr_az),1)/MCtimes;
AOA_el_align_dir_mean(ss) = sum((AOA_el_error_dir(:,ss)/pi*180)<(critical_width/Nr_el),1)/MCtimes;
AOD_az_align_dir_mean(ss) = sum((AOD_az_error_dir(:,ss)/pi*180)<(critical_width/Nt_az),1)/MCtimes;
AOD_el_align_dir_mean(ss) = sum((AOD_el_error_dir(:,ss)/pi*180)<(critical_width/Nt_el),1)/MCtimes;


Align_comp_mean(ss) = sum(((AOD_az_error_comp(:,ss)/pi*180)<(critical_width/Nt_az)&...
                           (AOD_el_error_comp(:,ss)/pi*180)<(critical_width/Nt_el)&...
                           (AOA_az_error_comp(:,ss)/pi*180)<(critical_width/Nr_az)&...
                           (AOA_el_error_comp(:,ss)/pi*180)<(critical_width/Nr_el)),1)/MCtimes;
                       
Align_sec_mean(ss) = sum((( AOD_az_error_sec(:,ss)/pi*180)<(critical_width/Nt_az)&...
                           (AOD_el_error_sec(:,ss)/pi*180)<(critical_width/Nt_el)&...
                           (AOA_az_error_sec(:,ss)/pi*180)<(critical_width/Nr_az)&...
                           (AOA_el_error_sec(:,ss)/pi*180)<(critical_width/Nr_el)),1)/MCtimes;

end

figure
% plot(SNR_range,AOA_az_align_comp_mean,'-s','linewidth',2);hold on
% plot(SNR_range,AOA_el_align_comp_mean,'-s','linewidth',2);hold on
% plot(SNR_range,AOD_az_align_comp_mean,'-s','linewidth',2);hold on
% plot(SNR_range,AOD_el_align_comp_mean,'-s','linewidth',2);hold on

plot(SNR_range,AOA_az_align_dir_mean,'-s','linewidth',2);hold on
plot(SNR_range,AOA_el_align_dir_mean,'-s','linewidth',2);hold on
plot(SNR_range,AOD_az_align_dir_mean,'-s','linewidth',2);hold on
plot(SNR_range,AOD_el_align_dir_mean,'-s','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Misalignment Rate')
legend('AoA (az)','AoA (el)','AoD (az)','AOD (el)')

% plot(SNR_range,AOAalign_comp_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_comp_mean,'--','linewidth',2);hold on

figure
plot(SNR_range,Align_comp_mean,'-o','linewidth',2);hold on
plot(SNR_range,Align_sec_mean,'--x','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Misalignment Rate')
legend('CSIA, mmMAGIC UMi LOS','DIA, mmMAGIC UMi LOS')
%% BF gain CDF
BW_data = 400e6;
noise_pow_dBm = -174 + 10*log10(BW_data) + NF - 30; % Data band noise pow in [dBW]

switch channel_model
    case 'SV'
    [b_PN,a_PN] = ecdf(20*log10(sqrt(Nt*Nr)*abs(data_mag_BF)./abs(data_mag_raw)));
    [b_sec,a_sec] = ecdf(20*log10(sqrt(Nt*Nr)*abs(data_mag_sec)./abs(data_mag_raw)));
    [b_dir,a_dir] = ecdf(20*log10(sqrt(Nt*Nr)*abs(data_mag_dir)./abs(data_mag_raw)));
    [b_true,a_true] = ecdf(20*log10(sqrt(Nt*Nr)*abs(data_mag_true)./abs(data_mag_raw)));
    case 'QuaDRiGa'
    [b_PN,a_PN] = ecdf(-noise_pow_dBm + 10*log10(abs(data_mag_BF(find(decision_CS)))));
    [b_sec,a_sec] = ecdf(-noise_pow_dBm + 10*log10(abs(data_mag_sec(find(decision_sec)))));
    [b_dir1,a_dir1] = ecdf(-noise_pow_dBm + 10*log10(abs(data_mag_dir1(find(decision_sec)))));
    [b_dir,a_dir] = ecdf(-noise_pow_dBm + 10*log10(abs(data_mag_dir(find(decision_sec)))));
    [b_true,a_true] = ecdf(-noise_pow_dBm + 10*log10(abs(data_mag_true)));
end

figure
plot(a_PN,b_PN,'linewidth',2);hold on
plot(a_sec,b_sec,'linewidth',2);hold on
plot(a_dir1,b_dir1,'linewidth',2);hold on
plot(a_dir,b_dir,'linewidth',2);hold on
plot(a_true,b_true,'linewidth',2);hold on

grid on
% xlim([-20,40])
xlabel('Data Phase Post-BF SNR [dB]')
ylabel('Prob(Gain<abscissa)')
legend('CSIA w/o RS','DIA  w/o RS','DIA  w/ 1 RS','DIA  w/ 2 RS','"Genie" Steering')
