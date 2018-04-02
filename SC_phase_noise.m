% SINR performance of single carrier waveform in DA and SA when phase noise
% is present

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
Nt = 32; % Number of Tx antennas (ULA)
Nr = 16; % Number of Tx antennas (ULA)
Nrf_DA = Nt; % Number of RF-chains in Digital Architecture
Nrf_FH = 2; % Number of RF-chains in Fully-connected Hybrid Architecture
fc = 28e9; % Carrier frequency
cluster_num = 2; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
BW = 1000e6; % Bandwidth of data symbols
Ts = 1/BW; % Sample duration
Nfft = 512; % Number of FFT in OFDM
sig_length = 1e3; % Number of samples required to generate
MCtimes = 20; % Number of Monte Carlo simulations

%-------------------------------------
% Phase Noise Specification
% FOM = L(f0,df) + 20log10(df/f0)+10log10(P_VCO)
%-------------------------------------
VCO_FOM = -114 + 20*log10(1e6/28e9) + 10*log10(27);
P_VCO = 27; % mW as unit
VCO_c = 10^(VCO_FOM/10)/P_VCO; % parameter c in PN specs
PN_sigma2 = VCO_c*4*pi^2*fc^2*Ts; % Time domain variance of PN Wiener process
PN_sweep_num = 10;
PN_range0 = linspace(10,40,PN_sweep_num);
PN_range = PN_sigma2* (10.^(PN_range0/10));

%-------------------------------------
% Results initialization
%-------------------------------------
IpN_PN = zeros(M, PN_sweep_num, MCtimes);
SINR_PN = zeros(M, PN_sweep_num, MCtimes);


%%
%------------- for loop of Monte Carlo simulation over MCtimes run -------------------
for MCindex = 1:MCtimes
    
    % Baseband Waveform; It provides a (sig_length by M) matrix with M unique streams
    QAM_level = 16; % Use X-QAM modulation
    upsam = 5; % Oversampling ratio
    sig_unitpow = get_waveform(QAM_level, sig_length, upsam, M);
    
    % Channel parameter generation (strongest two)
    [raygain,...
     raydelay,...
     rayAOA,...
     rayAOD ]= get_chan_parameter(plot_ellipse, print_stat, cluster_num, ray_num);
    
    
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
    H_freq_SC = squeeze(H_freq0(:,:,1)); % For SC waveform only the first carrier is used
    
    
    % Precoder and combiner matrix in digital architecture
    [U, Sigma, V] = svd(H_freq_SC);
    precoding_mtx = V(:,1:M);
    combining_mtx = U(:,1:M);
    
    % Precoder in hybrid architecture
    [ F_RF, F_BB ] = get_SC_hybrid_beamformer(precoding_mtx,Nt,Nrf_FH,32);
    Hybrid_Loss = norm(F_RF * F_BB - precoding_mtx,'fro')/norm(precoding_mtx,'fro');
    
    
    for pp=1:PN_sweep_num % Simulate different PN power
        
        PN_sigma = sqrt(PN_range(pp)); % Pick PN power and take square root

        % Signal initilizations
        array_sig_DA = zeros(Nt, sig_length);
        array_sig_FH = zeros(Nt, sig_length);
        
        PN_seq_DA = zeros(Nrf_DA, sig_length);
        PN_seq_DA(:,1) = 1;
        for ll=1:sig_length
            array_sig_DA_PN(:,ll) = (precoding_mtx * sig_unitpow(ll,:).').* PN_seq_DA(:,ll);
            PN_seq_DA(:,ll+1) = PN_seq_DA(:,ll).*exp(1j*randn(Nrf_DA,1)*PN_sigma);
        end
        array_sig_DA = array_sig_DA_PN./norm(array_sig_DA_PN,'fro')*sqrt(sig_length*tx_outpower);
        
        PN_seq_FH = zeros(Nrf_FH, sig_length);
        PN_seq_FH(:,1) = 1;
        for ll=1:sig_length
            array_sig_FH_PN(:,ll) = F_RF * ((F_BB * sig_unitpow(ll,:).').* PN_seq_FH(:,ll));
            PN_seq_FH(:,ll+1) = PN_seq_FH(:,ll).*exp(1j*randn(Nrf_FH,1)*PN_sigma);
        end
        array_sig_FH = array_sig_FH_PN./norm(array_sig_FH_PN,'fro')*sqrt(sig_length*tx_outpower);
        

        %% Signal Distortion at Specific Angle (stream 1 & Angle 1)
        noise_pow = norm(H_freq_SC,'fro')^2/Nt * awgn_pow;
        AWGN = (randn(sig_length,1)+1j*randn(sig_length,1))/sqrt(2)*sqrt(noise_pow);
        
        for mm=1:M
            sig_pow_norm = sig_unitpow(:,mm)./norm(sig_unitpow(:,mm))*sqrt(sig_length);
            
            rx_sig_DA = ((combining_mtx(:,mm)'* H_freq_SC) * array_sig_DA).' + AWGN;
            %  Normalize the Beamforming Gain (real gain)
            %  ------- min E||sig_org - sig_rx * alpha||^2 -----------
            gain_hat_DA = pinv(sig_pow_norm) * rx_sig_DA;
            gain_sig_DA(pp,MCindex) = abs(gain_hat_DA)^2;
            IpN_DA(mm,pp,MCindex) = norm(gain_hat_DA * sig_pow_norm - rx_sig_DA)^2/sig_length;
            SINR_DA(mm,pp,MCindex) = abs(gain_hat_DA)^2/(IpN_DA(mm,pp,MCindex));
            
            rx_sig_FH = ((combining_mtx(:,mm)'* H_freq_SC) * array_sig_FH).' + AWGN;
            %  Normalize the Beamforming Gain (real gain)
            %  ------- min E||sig_org - sig_rx * alpha||^2 -----------
            gain_hat_FH = pinv(sig_pow_norm) * rx_sig_FH;
            gain_sig_FH(pp,MCindex) = abs(gain_hat_FH)^2;
            IpN_FH(mm,pp,MCindex) = norm(gain_hat_FH * sig_pow_norm - rx_sig_FH)^2/sig_length;
            SINR_FH(mm,pp,MCindex) = abs(gain_hat_FH)^2/(IpN_FH(mm,pp,MCindex));
        end
    end
end
%% Take mean among MC realizations
for mm=1:1
    for bb=1:length(PN_range)
        SINR_mean_DA(mm,bb) = mean(squeeze(SINR_DA(mm,bb,:)));
        SINR_mean_FH(mm,bb) = mean(squeeze(SINR_FH(mm,bb,:)));
    end
end

%% SINR figure
figure
semilogx(PN_range, 10*log10(SINR_mean_DA(1,:)),'-o','linewidth',2);hold on
semilogx(PN_range, 10*log10(SINR_mean_FH(1,:)),'-o','linewidth',2);
grid on
xlabel('PN Power (\sigma^2)')
ylabel('SINR (dB)')
legend('DA','FH')