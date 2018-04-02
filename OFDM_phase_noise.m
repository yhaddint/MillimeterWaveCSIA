% test whether it is possible to get reasonable "initial guess" of channle
% parameter
%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;close all
rng(2)
plot_ellipse = 0; % plot geolocations of Tx/Rx/Scatterers
print_stat = 0; % print channel parameter summary

%-------------------------------------
% System Parameters
%-------------------------------------
awgn_pow = 1; % Receiver thermal noise power;
tx_outpower = 1; % Used to control transmission power 
M = 2; % Number of unique data streams
Nt = 16; % Number of Tx antennas (ULA)
Nr = 16; % Number of Tx antennas (ULA)
Nrf_DA = Nt; % Number of RF-chains in Digital Architecture
Nrf_FH = 2; % Number of RF-chains in Fully-connected Hybrid Architecture
fc = 28e9; % Carrier frequency
cluster_num = 2; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
BW = 1000e6; % Bandwidth of data symbols
Ts = 1/BW; % Sample duration
Nfft = 512; % Number of FFT in OFDM
OFDM_blk_num = 1e1; % Number of OFDM blocks (each has Nfft symb) in simulation
noise_pow = 1; % AWGN power
MCtimes = 10; % Number of Monte Carlo simulations

%-------------------------------------
% Phase Noise Specification
% FOM = L(f0,df) + 20log10(df/f0)+10log10(P_VCO)
%-------------------------------------
VCO_FOM = -114 + 20*log10(1e6/28e9) + 10*log10(27);
P_VCO = 27; % mW as unit
VCO_c = 10^(VCO_FOM/10)/P_VCO; % parameter c in PN specs
PN_sigma2 = VCO_c*4*pi^2*fc^2*Ts; % Time domain variance of PN Wiener process
PN_sweep_num = 10;
PN_range0 = linspace(-10,20,PN_sweep_num);
PN_range = PN_sigma2* (10.^(PN_range0/10));


%-------------------------------------
% Results initialization
%-------------------------------------
IpN_PN = zeros(OFDM_blk_num, PN_sweep_num, MCtimes);
SINR_PN = zeros(OFDM_blk_num, PN_sweep_num, MCtimes);


%%
%------------- for loop of Monte Carlo simulation over MCtimes run -------------------
for MCindex = 1:MCtimes
    
    % Baseband Waveform of OFDM symbols;
    % It provides a (M * Nfft by sig_length) matrix with M unique streams
    sig_unitpow = ((randi(2,Nfft*M, OFDM_blk_num)*2-3)+1j*(randi(2,Nfft*M, OFDM_blk_num)*2-3))/sqrt(2);
    
    % Channel parameter generation (strongest two)
    [raygain,...
     raydelay,...
     rayAOA,...
     rayAOD ]= get_chan_parameter(plot_ellipse,...
                                  print_stat,...
                                  cluster_num,...
                                  ray_num);
    
    
    % MIMO channel generation and normalization
    H_freq0 = get_H_freq( raygain,...
                       raydelay,...
                       rayAOA,...
                       rayAOD,...
                       cluster_num,...
                       ray_num,...
                       Nt, Nr,...
                       fc, Nfft, Ts);
    norm_factor = sqrt(mean(mean(mean(abs(H_freq0).^2))));
    % norm_factor = 1;
    H_freq0 = H_freq0 / norm_factor;
    
    for pp=1:PN_sweep_num % Simulate different PN power
        
        PN_sigma = sqrt(PN_range(pp)); % Pick PN power and take square root
%         PN_sigma = 0; % turn-off phase noise in debug

        % Signal zero-initilizations
        array_sig_DA = zeros(Nt, OFDM_blk_num);
%         array_sig_FH = zeros(Nt, sig_length);
        PN_DA_next_blk = ones(Nrf_DA,1); % Used to connect PN between OFDM block
        
        for ll=1:OFDM_blk_num
            
            % Generate PN sequence within OFDM block
            PN_seq_DA = zeros(Nrf_DA, Nfft);
            PN_seq_DA(:,1) = PN_DA_next_blk;
            for nn=1:Nfft-1
                PN_seq_DA(:,nn+1) = PN_seq_DA(:,nn).*exp(1j*randn(Nrf_DA,1)*PN_sigma);
            end
            PN_seq_DA(:,nn) = PN_DA_next_blk;
            
            % Compute beamformer in digital architecture based on SVD
            [precoder_DA, combiner_DA] = get_WB_DA_beamformer( H_freq0, M );
            
            % Run digital precoding algorithm in DA
            array_sig_DA_raw(:,ll) = run_DA_WB_Precoding( precoder_DA, sig_unitpow(:,ll), M );
            
            % Apply phase noise in digitally precoder symbols in DA
            array_sig_DA_PN(:,ll) = get_OFDM_ICI( array_sig_DA_raw(:,ll), PN_seq_DA );
            
        end
          array_sig_DA = array_sig_DA_PN;

%         array_sig_DA = array_sig_DA_PN./norm(array_sig_DA_PN,'fro')*sqrt(OFDM_blk_num*tx_outpower);
        
%         PN_seq_FH = zeros(Nrf_FH, sig_length);
%         PN_seq_FH(:,1) = 1;
%         for ll=1:sig_length
%             array_sig_FH_PN(:,ll) = F_RF * ((F_BB * sig_unitpow(ll,:).').* PN_seq_FH(:,ll));
%             PN_seq_FH(:,ll+1) = PN_seq_FH(:,ll).*exp(1j*randn(Nrf_FH,1)*PN_sigma);
%         end
%         array_sig_FH = array_sig_FH_PN./norm(array_sig_FH_PN,'fro')*sqrt(sig_length*tx_outpower);
        
        
        %% Signal Distortion at Specific Angle (stream 1 & Angle 1)
        noise_pow = awgn_pow;
        AWGN = (randn(Nr*Nfft,OFDM_blk_num)+1j*randn(Nr*Nfft,OFDM_blk_num))/sqrt(2)*sqrt(noise_pow);
        AWGN(:,1) = zeros(Nr*Nfft,1); % Removed AWGN in first OFDM blk (noiseless pilots)
%         AWGN(:,2:end) = ones(Nr*Nfft,OFDM_blk_num-1)*1e-10; % Trun-off AWGN to debug

%         sig_pow_norm = sig_unitpow./norm(sig_unitpow,'fro')*sqrt(OFDM_blk_num);
        
        for ll = 1:OFDM_blk_num
            % Zero initialization
            rx_sig_DA = zeros(Nr*Nfft,1);
            
            % Apply channel to transmitted signal at each subcarrier
            for kk=1:Nfft
                sigtx_index = (kk-1)*Nt+1:kk*Nt;
                sigrx_index = (kk-1)*Nr+1:kk*Nr;
                rx_sig_DA(sigrx_index) = ...
                    squeeze(H_freq0(:,:,kk)) * array_sig_DA(sigtx_index, ll) + AWGN(sigrx_index, ll);
            end

            % Apply Rx combiner to all subcarriers
            rx_comb_sig_DA(:,ll) = run_DA_WB_Combining( combiner_DA, rx_sig_DA, M );
        end
        
        % Using the first OFDM block as training stage (assuemd to be noiseless)
        gain_hat_DA = rx_comb_sig_DA(:,1)./sig_unitpow(:,1);
        
        % Evaluate SINR of symbols in the rest of OFDM blocks
        for ll=1:OFDM_blk_num
            if ll>1
                IpN_DA(ll,pp,MCindex) = norm(gain_hat_DA .* sig_unitpow(:,ll) - rx_comb_sig_DA(:,ll))^2;
                SINR_DA(ll,pp,MCindex) = norm(gain_hat_DA)^2/(IpN_DA(ll,pp,MCindex));
            end
        end

    %             rx_sig_FH = ((combining_mtx(:,mm)'* H_freq_SC) * array_sig_FH).' + AWGN;
    %             %  Normalize the Beamforming Gain (real gain)
    %             %  ------- min E||sig_org - sig_rx * alpha||^2 -----------
    %             gain_hat_FH = pinv(sig_pow_norm) * rx_sig_FH;
    %             gain_sig_FH(pp,MCindex) = abs(gain_hat_FH)^2;
    %             IpN_FH(mm,pp,MCindex) = norm(gain_hat_FH * sig_pow_norm - rx_sig_FH)^2/OFDM_blk_num;
    %             SINR_FH(mm,pp,MCindex) = abs(gain_hat_FH)^2/(IpN_FH(mm,pp,MCindex));
    end
end
%% Take mean among MC realizations
for ll=1:OFDM_blk_num
    for bb=1:length(PN_range)
        SINR_mean_DA(ll,bb) = mean(squeeze(SINR_DA(ll,bb,:)));
%         SINR_mean_FH(mm,bb) = mean(squeeze(SINR_FH(mm,bb,:)));
    end
end

%% SINR figure
figure
semilogx(PN_range, 10*log10(SINR_mean_DA(2,:)),'-o','linewidth',2);hold on
semilogx(PN_range, 10*log10(SINR_mean_DA(10,:)),'-o','linewidth',2);hold on

% semilogx(PN_range, 10*log10(SINR_mean_FH(1,:)),'-o','linewidth',2);
grid on
xlabel('PN Power (\sigma^2)')
ylabel('SINR (dB)')
legend('DA 2','DA 10')