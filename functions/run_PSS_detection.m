function [ peak_pow_H1, peak_pow_H0 ] = run_PSS_detection( SNR_dB, STO )
%RUN_PSS_DETECTION Summary of this function goes here
%   Detailed explanation goes here
    
    % control parameters ---------
    if STO==0
        perfectSTO=1;
    else
        perfectSTO=0;
    end
    
    % ----- system parameters --------
    M = 20; % Number of ZC burst
    noise_pow = 10^(-SNR_dB/10);
    ZC_root = 29; % ZC root, a prime value with ZC length
    ZC_N = 255; % ZC sequence length
    burst_N = ZC_N*2; % number of sample in each burst (2OFDM symbol for now)
    Rx_sig_length = burst_N * M + ZC_N - 1; % signal length after ZC correlation;
    Nt = 32; % Number of antenna in tx
    Nr = 16; % Number of antenna in Rx

    % ------ channel parameter ------
    print_stat = 0;
    cluster_num = 2;
    ray_num = 20;
    sigma_delay_spread = 0;
    centroid_AOA = 'random';
    sigma_AOA_spread = 0;
    centroid_AOD = 'random';
    sigma_AOD_spread = 0;

    % ------ ZC sequence generation ---------
    seq = lteZadoffChuSeq(ZC_root,ZC_N); % ZC symbol mapped into OFDM subcarriers
    seq_1DC = ifft(seq)*sqrt(ZC_N+1); % Time domain signal used to ZC detection & STO estimation; DC subcarrier is not null
    burst_sig = [seq_1DC; zeros(burst_N-ZC_N,1)]; % each burst has one ZC and something else (SSS/other control info)
    Tx_sig = repmat(burst_sig,M,1);
    burst_length = burst_N; % Number of samples in one ZC burst
    Tx_sig_length = length(Tx_sig); % Number of samples in M ZC burst
    ZC_t_domain = conj(flipud(seq_1DC));  % ZC sequence used for correlation in Rx



    % ------- Random parameter realizations -----   
    [ raygain,...
      raydelay,...
      ray_AOA_azim,...
      ray_AOD_azim ] = get_chan_parameter_nogeo(print_stat,...
                                          cluster_num,...
                                          ray_num,...
                                          sigma_delay_spread,...
                                          centroid_AOA,...
                                          sigma_AOA_spread,...
                                          centroid_AOD,...
                                          sigma_AOD_spread); % Channel parameters generation
    W_mat = get_IA_BF(Nr, M,'PN');% Rx beamformer in IA stage
    V_mat = get_IA_BF(Nt, M,'PN'); % Tx beamformer in IA stage
    
    % ----- Initializations of vec, mat -----
    Rx_sig = zeros(Tx_sig_length, 1); % received signal in t domain
    corr_out_H1 = zeros(Rx_sig_length + burst_N - ZC_N + 1, 1); % pad a zero at the end
    corr_out_H0 = zeros(Rx_sig_length + burst_N - ZC_N + 1,1); % pad a zero at the end

    
    % ------- Channel Generation --------
    H_chan = get_H_NB(raygain,...
                      ray_AOA_azim,...
                      ray_AOD_azim,...
                      cluster_num,...
                      ray_num,...
                      Nt, Nr); % Generate discrete time domain frequency-flat channel
    H_chan0 = H_chan./norm(H_chan,'fro')*sqrt(Nt*Nr);
    
    % ----- received signal generation ------
    precoder_index_old = 0;
    combiner_index_old = 0;
    for nn=1:Tx_sig_length
        precoder_index = floor( (nn-1) / burst_length )+1;
        combiner_index = floor( (nn + STO - 1) / burst_length )+1;
        if (precoder_index ~= precoder_index_old) || (combiner_index ~= combiner_index_old)
            w_vec = W_mat(:,combiner_index);
            v_vec = V_mat(:,precoder_index);
            g_effective = (w_vec'*H_chan*v_vec);
            precoder_index_old = precoder_index;
            combiner_index_old = combiner_index;
        end
        Rx_sig(nn) = g_effective * Tx_sig(nn);
    end
    
    % ------- AWGN -------
    noise = (randn(Tx_sig_length,1)+1j*randn(Tx_sig_length,1))/sqrt(2)*sqrt(noise_pow);
    Rx_sig_H1 = Rx_sig + noise;
    Rx_sig_H0 = noise;
    
    % ------ T Domain ZC Correlation -------
    corr_out_H1(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H1)).^2; % corr rx t-domain sig with ZC
    corr_out_H0(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H0)).^2; % corr rx t-domain sig with ZC
    
    % ----- Multi-Peak Detection ---------
    if perfectSTO % Concept scenario where peak locations is know
        peak_pow_H1 = sum(corr_out_H1(ZC_N:burst_length:end));
        peak_pow_H0 = sum(corr_out_H0(ZC_N:burst_length:end));
    else % Practical scenario where peak location is unknown
        peak_pow_H1 = max(sum(reshape(corr_out_H1,burst_length,M+1),2));
        peak_pow_H0 = max(sum(reshape(corr_out_H0,burst_length,M+1),2));
    end

end

