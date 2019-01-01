function [  peak_pow_H1,...
            peak_pow_H0,...
            peakindex_H1 ] = run_PSS_detection_CFO( SNR_range,...
                                                STOtype,...
                                                STOinfo,...
                                                BFtype,...
                                                M_burst,...
                                                CFO,...
                                                Tx_sec_codebook,...
                                                Rx_sec_codebook)
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
% STOinfo: 0 or 1 to indicate whether system has timing info, different
%          algs. are used in two cases.
% BFtype: use 'PN', 'directional', or 'sector'
% M_burst: number of burst in discovery
% CFO: Actual CFO value in the unit of ppm
    
    % Convert CFO ppm into phase rotate in discrete time 
    CFO_samp = CFO*28e3/(2*28.8e6)*2*pi; %with unit rad per sample
    
    % ------ SNR parameters ---------
    SNR_num = length(SNR_range);

    % ---- Results vector dimension initialization -----
    peak_pow_H1 = zeros(SNR_num,1);
    peak_pow_H0 = zeros(SNR_num,1);
    
    % ----- system parameters --------
    Nt = 256; % Number of antenna in tx (ULA)
    Nr = 32; % Number of antenna in Rx (ULA)
    CP = 8; % cyclic prefix in [samples]
    Nc = 4; % maximum multipath delay in [samples]
    SC = 256; % number of subcarrier in IA SS burst
    STO_max = 1000; % range of maximum integer offset
    OFDM_sym_num = 4;
    
    switch BFtype
        case 'directional'
            M = Nt*Nr;
        case 'PN'
            M = M_burst(1)*M_burst(2);
        case 'sector_LS'
            M = M_burst(1)*M_burst(2);
        case 'sector_FSM_KW'
            M = M_burst(1)*M_burst(2);
    end
    
    % ---- signal length parameters --------
    ZC_root = 29; % ZC root, a coprime number with ZC length
    ZC_N = 127; % ZC sequence length
    burst_N = SC*OFDM_sym_num; % number of sample in each burst (4OFDM symbol for now)
%     burst_N = ZC_N*OFDM_sym_num + CP; % number of sample in each burst (4OFDM symbol for now)
    
    % ----- STO ------
    switch STOtype
        case 'zero'
            STO = 0;
        case 'random'
            STO = randi(burst_N * M + ZC_N - 1);
            % I've studied scenario with unbounded STO; It works but I'm not
            % going to put it in the paper.
        otherwise
            STO = STOtype;
    end

    % received sample number 
    Rx_sig_length = burst_N * M + ZC_N - 1 + STO; % signal length after ZC correlation;

    
    % ------ channel parameter ------
    print_stat = 0; % print AoA/AoD statistics
    cluster_num = 4; % 2 or 3 or 4 (not 5!)
    ray_num = 1; %20 % rays within a cluster
    sigma_delay_spread = 0;
    centroid_AOA = 'random';
    sigma_AOA_spread = 0; % clustered sparse channel is left for future
    centroid_AOD = 'random';
    sigma_AOD_spread = 0;

    % ------ ZC sequence generation ---------
    seq = lteZadoffChuSeq(ZC_root,ZC_N); % ZC symbol mapped into OFDM subcarriers
    seq_1DC = ifft(seq)*sqrt(ZC_N+1); % Time domain signal used to ZC detection & STO estimation; DC subcarrier is not null
%     burst_sig = [seq_1DC; zeros(burst_N-ZC_N,1)]; % each burst has one ZC and something else (SSS/other control info)
    burst_sig = [seq_1DC(end-CP+1:end);...
                 seq_1DC;...
                 zeros(burst_N-(ZC_N+CP),1)]; % each burst has one ZC (CP-ed) and something else (SSS/other control info)
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
        case 'sector_LS'
            V_mat = Tx_sec_codebook;
            W_mat = Rx_sec_codebook;
        case 'sector_FSM_KW'
            V_mat = Tx_sec_codebook;
            W_mat = Rx_sec_codebook;
            
    end % end of switch
    
  
    % generate NB channel and received signal for each delay tap
    Nc_acc_mtx = toeplitz([1,zeros(1,burst_length-Nc)]',[ones(1,Nc),zeros(1,burst_length-Nc)]);
    Rx_sig0 = zeros(burst_length*M,cluster_num);
    for path_index = 1:cluster_num
        % ------- Channel Generation --------
        H_chan = get_H_NB(raygain(path_index),...
                          ray_AOA_azim(path_index),...
                          ray_AOD_azim(path_index),...
                          1,...
                          ray_num,...
                          Nt, Nr); % Generate discrete time domain frequency-flat channel
        H_chan0 = H_chan./norm(H_chan,'fro')*sqrt(Nt*Nr/cluster_num); % H per multipath

        % ----- received signal generation ------
        precoder_index_old = 0;
        combiner_index_old = 0;
        for nn=1:Tx_sig_length
            
% I've turned off switch since it saves time in large scale sim; get it
% back later!
%             switch BFtype 
%                 case 'PN' % w/ PN BF, BF changes every burst_length samples 
                precoder_index = floor( (nn-1) / burst_length )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,M)+1;

%                 case 'directional' % w/ directional beam (steering vector) 
%                 precoder_index = floor( (nn-1) / (burst_length*Nr) )+1;
%                 combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
%                 combiner_index = mod(combiner_index_raw-1,Nr)+1;

%                 case 'sector_LS' % w/ sector beam 
%                 precoder_index = floor( (nn-1) / (burst_length*M_burst(2)) )+1;
%                 combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
%                 combiner_index = mod(combiner_index_raw-1,M_burst(2))+1;
                
%                 case 'sector_FSM_KW' % w/ sector beam 
%                 precoder_index = floor( (nn-1) / (burst_length*M_burst(2)) )+1;
%                 combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
%                 combiner_index = mod(combiner_index_raw-1,M_burst(2))+1;
%             end
            if (precoder_index ~= precoder_index_old) || (combiner_index ~= combiner_index_old)
    %             fprintf('precoder index %d, ',precoder_index);
    %             fprintf('combiner index %d\n',combiner_index); 
                w_vec = W_mat(:,combiner_index);
                v_vec = V_mat(:,precoder_index);
                g_effective = (w_vec'*H_chan0*v_vec);
                precoder_index_old = precoder_index;
                combiner_index_old = combiner_index;
            end
%             index_debug(:,nn) = [precoder_index;combiner_index];
            g_save_debug(nn) = g_effective;
%             Rx_sig0(nn,path_index) = g_effective * Tx_sig(nn);
        end % end of sample sweeping
        Rx_sig0(:,path_index) = g_save_debug.' .* Tx_sig;
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
    noise_at_STO = (randn(STO,1)+1j*randn(STO,1))/sqrt(2);
    
    % ----- Initializations of vec, mat -----
%     Rx_sig = zeros(Tx_sig_length, 1); % received signal in t domain
    corr_out_H1 = zeros(Rx_sig_length + burst_N - ZC_N + 1, 1); % pad a zero at the end
    corr_out_H0 = zeros(Rx_sig_length + burst_N - ZC_N + 1,1); % pad a zero at the end
  
    % ------- sweep all SNR values -------
    for ss = 1:SNR_num
        noise_pow = 10^(-SNR_range(ss)/10);
        awgn = noise * sqrt(noise_pow);
        Rx_sig_H1 = Rx_sig.*exp(1j*CFO_samp*(0:length(Rx_sig)-1).') + awgn ;
        Rx_sig_H0 = awgn;

        % ------ T Domain ZC Correlation -------
        if STO==0
            corr_out_H1(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H1)/ZC_N).^2; % corr rx t-domain sig with ZC
            corr_out_H0(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H0)/ZC_N).^2; % corr rx t-domain sig with ZC
        else
            Rx_sig_H0_wSTO = [noise_at_STO * sqrt(noise_pow); Rx_sig_H0];
            Rx_sig_H1_wSTO = [noise_at_STO * sqrt(noise_pow); Rx_sig_H1];
            corr_out_H1_STO = abs(conv(ZC_t_domain,Rx_sig_H1_wSTO)/ZC_N).^2; % corr rx t-domain sig with ZC
            corr_out_H0_STO = abs(conv(ZC_t_domain,Rx_sig_H0_wSTO)/ZC_N).^2; % corr rx t-domain sig with ZC
        end

        
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
                end

            % ----- Single-Peak Detection with directional beams---------
            case 'directional'
                peak_pow_H1(ss) = max(corr_out_H1);
                peak_pow_H0(ss) = max(corr_out_H0);
            
            % ----- Single-Peak Detection with directional beams---------
            otherwise
                if STOinfo % Concept scenario where peak locations is know
                    peak_pow_H1(ss) = max(corr_out_H1(ZC_N+CP:burst_length:end));
                    peak_pow_H0(ss) = max(corr_out_H0(ZC_N+CP:burst_length:end));
                    peakindex_H1(ss) = 0;
                else % Practical scenario where peak location is unknown
                    peak_pow_H1(ss) = max(corr_out_H1_STO);
                    peak_pow_H0(ss) = max(corr_out_H0_STO);
                    peakindex_H1(ss) = 0;
%                     [peak_pow_H1(ss) peakindex_H1(ss)] = max(mean(reshape(corr_out_H1,burst_length,M+1),2));
%                     peak_pow_H0(ss) = max(mean(reshape(corr_out_H0,burst_length,M+1),2));

                end
        end % end of switch
    end % end of SNR sweeping

end

