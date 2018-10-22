function [  peak_pow_H1,...
            peak_pow_H0,...
            peakindex_H1 ] = run_PSS_detection_CFO( SNR_range,...
                                                STOtype,...
                                                STOinfo,...
                                                BFtype,...
                                                M_burst,...
                                                CFO )
%RUN_PSS_DETECTION Summary of this function goes here
%[  peak_pow_H1,...
%   peak_pow_H0,...
%   peakindex_H1 ] = run_PSS_detection_CFO( SNR_range,...
%                                           STOtype,...
%                                           STOinfo,...
%                                           BFtype,...
%                                           M_burst,...
%                                           CFO )
% ------ Input explainations:------------
% SNR_range: SNR range in dB scale
% STOtype: 'zero' or 'random'
% STOinfo: 0 or 1 to indicate whether system has timing info
% BFtype: use 'PN', 'directional', or 'sector'
% M_burst: number of burst in discovery
% CFO: Actual CFO value in the unit of ppm
    
    % Convert CFO ppm into phase rotate in discrete time 
    CFO_samp = CFO*28e3/57.6e6*2*pi; %with unit rad per sample
    
    % ------ SNR parameters ---------
    SNR_num = length(SNR_range);

    % ---- Results vector dimension initialization -----
    peak_pow_H1 = zeros(SNR_num,1);
    peak_pow_H0 = zeros(SNR_num,1);
    
    % ----- system parameters --------
    Nt = 256; % Number of antenna in tx
    Nr = 32; % Number of antenna in Rx
    CP = 8; % cyclic prefix
    Nc = 4; % maximum multipath delay in [Ts]
    
    switch BFtype
        case 'directional'
            M = Nt*Nr;
        case 'PN'
            M = M_burst(1)*M_burst(2);
        case 'sector'
            M = M_burst(1)*M_burst(2);
    end
    ZC_root = 29; % ZC root, a coprime number with ZC length
    ZC_N = 127; % ZC sequence length
%     burst_N = ZC_N*2; % number of sample in each burst (2OFDM symbol for now)
    burst_N = ZC_N*2 + CP; % number of sample in each burst (2OFDM symbol for now)
    Rx_sig_length = burst_N * M + ZC_N - 1; % signal length after ZC correlation;


    % ----- STO ------
    switch STOtype
        case 'zero'
            STO = 0;
        case 'random'
            STO = randi(Rx_sig_length);
    end

    % ------ channel parameter ------
    print_stat = 0;
    cluster_num = 3;
    ray_num = 1; %20 % rays within a cluster
    sigma_delay_spread = 0;
    centroid_AOA = 'random';
    sigma_AOA_spread = 0;
    centroid_AOD = 'random';
    sigma_AOD_spread = 0;

    % ------ ZC sequence generation ---------
    seq = lteZadoffChuSeq(ZC_root,ZC_N); % ZC symbol mapped into OFDM subcarriers
    seq_1DC = ifft(seq)*sqrt(ZC_N+1); % Time domain signal used to ZC detection & STO estimation; DC subcarrier is not null
%     burst_sig = [seq_1DC; zeros(burst_N-ZC_N,1)]; % each burst has one ZC and something else (SSS/other control info)
    burst_sig = [seq_1DC(end-CP+1:end);...
                 seq_1DC;...
                 zeros(burst_N-(ZC_N+CP),1)]; % each burst has one ZC (CPed) and something else (SSS/other control info)
    Tx_sig = repmat(burst_sig,M,1);
%     burst_length = burst_N; % Number of samples in one ZC burst
    burst_length = length(burst_sig);
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
      delay_cand = randperm(Nc-1);
      pathdelay = [0,delay_cand(1:cluster_num-1)];
      % Multipath in sample delay domain
      raygain = sqrt(1/cluster_num)*exp(1j*rand(cluster_num,1)*2*pi);
      % I don't want random raygain in this paper; I need equal gain path
    
    switch BFtype
        case 'PN'
            V_mat = get_IA_BF(Nt, M, BFtype); % Tx beamformer in IA stage
            W_mat = get_IA_BF(Nr, M, BFtype); % Rx beamformer in IA stage
        case 'directional'
            V_mat = get_IA_BF(Nt, Nt, BFtype); % Tx beamformer in IA stage
            W_mat = get_IA_BF(Nr, Nr, BFtype); % Rx beamformer in IA stage
        case 'sector'
            V_mat = get_IA_BF(Nt, M_burst(1), BFtype); % Tx beamformer in IA stage
            W_mat = get_IA_BF(Nr, M_burst(2), BFtype); % Rx beamformer in IA stage
            
    end % end of switch
    
    % ----- Initializations of vec, mat -----
    Rx_sig = zeros(Tx_sig_length, 1); % received signal in t domain
    corr_out_H1 = zeros(Rx_sig_length + burst_N - ZC_N + 1, 1); % pad a zero at the end
    corr_out_H0 = zeros(Rx_sig_length + burst_N - ZC_N + 1,1); % pad a zero at the end
    
    % generate NB channel and received signal for each delay tap
    Rx_sig0 = zeros(burst_length*M,cluster_num);
    for path_index = 1:cluster_num
        % ------- Channel Generation --------
        H_chan = get_H_NB(raygain(path_index),...
                          ray_AOA_azim(path_index),...
                          ray_AOD_azim(path_index),...
                          1,...
                          ray_num,...
                          Nt, Nr); % Generate discrete time domain frequency-flat channel
        H_chan0 = H_chan./norm(H_chan,'fro')*sqrt(Nt*Nr/cluster_num);

        % ----- received signal generation ------
        precoder_index_old = 0;
        combiner_index_old = 0;
        for nn=1:Tx_sig_length
            switch BFtype 
                case 'PN' % w/ PN BF, BF changes every burst_length samples 
                precoder_index = floor( (nn-1) / burst_length )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,M)+1;

                case 'directional' % w/ directional beam (steering vector) 
                precoder_index = floor( (nn-1) / (burst_length*Nr) )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,Nr)+1;

                case 'sector' % w/ sector beam 
                precoder_index = floor( (nn-1) / (burst_length*M_burst(2)) )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,M_burst(2))+1;
            end
            if (precoder_index ~= precoder_index_old) || (combiner_index ~= combiner_index_old)
    %             fprintf('precoder index %d, ',precoder_index);
    %             fprintf('combiner index %d\n',combiner_index); 
                w_vec = W_mat(:,combiner_index);
                v_vec = V_mat(:,precoder_index);
                g_effective = (w_vec'*H_chan0*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
            index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
            Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
    end
    
    % ----- summation over all delay tap --------
    Rx_sig = zeros(burst_length*M,1);
    for path_index = 1:cluster_num
        timewindow0 = (1+pathdelay(path_index)):burst_length*M;
        timewindow1 = 1:(burst_length*M-pathdelay(path_index));
        Rx_sig(timewindow0,1) = Rx_sig(timewindow0,1) + Rx_sig0(timewindow1,path_index);
    end
    
    % ------- AWGN -------
    noise = (randn(Tx_sig_length,1)+1j*randn(Tx_sig_length,1))/sqrt(2);
    
    % ------- sweep all SNR values -------
    for ss = 1:SNR_num
        noise_pow = 10^(-SNR_range(ss)/10);
        awgn = noise * sqrt(noise_pow);
        Rx_sig_H1 = Rx_sig.*exp(1j*CFO_samp*(0:length(Rx_sig)-1).') + awgn ;
        Rx_sig_H0 = awgn;

        % ------ T Domain ZC Correlation -------
        corr_out_H1(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H1)/ZC_N).^2; % corr rx t-domain sig with ZC
        corr_out_H0(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H0)/ZC_N).^2; % corr rx t-domain sig with ZC

        
        switch BFtype
            % ----- Multi-Peak Detection ---------
            case 'PN'
                if STOinfo % Concept scenario where peak locations is know
                    peak_window0 = [];
                    for ncindex = 1:Nc
                        peak_window0 = [peak_window0;(ZC_N+CP+ncindex-1:burst_length:burst_length*M)'];
                    end
                    peak_window = sort(peak_window0,'ascend');
%                     peak_pow_H1_t1(ss) = sum(corr_out_H1(ZC_N+CP:burst_length:end))/M;
                    peak_pow_H1(ss) = sum(corr_out_H1(peak_window))/M;
                    peak_pow_H0(ss) = sum(corr_out_H0(peak_window))/M;
                    peakindex_H1(ss) = 0;
                else % Practical scenario where peak location is unknown
                    [peak_pow_H1(ss) peakindex_H1(ss)] = max(mean(reshape(corr_out_H1,burst_length,M+1),2));
                    peak_pow_H0(ss) = max(mean(reshape(corr_out_H0,burst_length,M+1),2));
                end

            % ----- Single-Peak Detection ---------
            case 'directional'
                peak_pow_H1(ss) = max(corr_out_H1);
                peak_pow_H0(ss) = max(corr_out_H0);
                
            case 'sector'
                if STOinfo % Concept scenario where peak locations is know
                    peak_pow_H1(ss) = max(corr_out_H1(ZC_N:burst_length:end));
                    peak_pow_H0(ss) = max(corr_out_H0(ZC_N:burst_length:end));
                    peakindex_H1(ss) = 0;
                else % Practical scenario where peak location is unknown
                    peak_pow_H1(ss) = max(corr_out_H1);
                    peak_pow_H0(ss) = max(corr_out_H0);
                    peakindex_H1(ss) = 0;
%                     [peak_pow_H1(ss) peakindex_H1(ss)] = max(mean(reshape(corr_out_H1,burst_length,M+1),2));
%                     peak_pow_H0(ss) = max(mean(reshape(corr_out_H0,burst_length,M+1),2));

                end
        end % end of switch
    end % end of SNR sweeping

end

