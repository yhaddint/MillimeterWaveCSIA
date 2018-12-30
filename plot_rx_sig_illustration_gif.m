clear;clc;
%
%%
% ------ script control parameters -------
rng(3)
plot_cdf = 0;

% parameter for evaluation of theoretical equations
Nt = 128;
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
SNR_dB = -10;
%% Monte Carlo evaluations
STO = 170; %[170 (3us), 568(10us), 960 (critical BF mismatch)]; Type of sample timing offset: value or 'zero'/'random'
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


% received sample number 
Rx_sig_length = burst_N * M + ZC_N - 1 + STO; % signal length after ZC correlation;


% ------ channel parameter ------
print_stat = 0; % print AoA/AoD statistics
cluster_num = 2; % 2 or 3 or 4 (not 5!)
ray_num = 1; %20 % rays within a cluster
sigma_delay_spread = 0;
centroid_AOA = [23.18,56.66]/180*pi; % one particular example
% centroid_AOA = [-5]/180*pi; % one particular example

sigma_AOA_spread = 0; % clustered sparse channel is left for future
% centroid_AOD = [-5]/180*pi; % one particular example
centroid_AOD = [6.10,-37.46]/180*pi; % one particular example

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
        switch BFtype 
            case 'PN' % w/ PN BF, BF changes every burst_length samples 
                precoder_index = floor( (nn-1) / burst_length )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,M)+1;

            case 'directional' % w/ directional beam (steering vector) 
                precoder_index = floor( (nn-1) / (burst_length*Nr) )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,Nr)+1;

            case 'sector_LS' % w/ sector beam 
                precoder_index = floor( (nn-1) / (burst_length*M_burst(2)) )+1;
                combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
                combiner_index = mod(combiner_index_raw-1,M_burst(2))+1;

            case 'sector_FSM_KW' % w/ sector beam 
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

% ------- setup noise power -------
noise_pow = 10^(-SNR_dB/10);
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
            peak_pow_H1 = sum(corr_out_H1(peak_window))/M;
            peak_pow_H0 = sum(corr_out_H0(peak_window))/M;
            peakindex_H1 = 0;
        else % Practical scenario where peak location is unknown
            post_corr_ED_H1 = sum(reshape([corr_out_H1_STO(ZC_N+CP:end);...
                zeros(burst_length*(M+2)-length(corr_out_H1_STO(ZC_N+CP:end)),1)],...
                burst_length,M+2),2)/M;

            post_corr_ED_H0 = sum(reshape([corr_out_H0_STO(ZC_N+CP:end);...
                zeros(burst_length*(M+2)-length(corr_out_H0_STO(ZC_N+CP:end)),1)],...    
                burst_length,M+2),2)/M;


            ave_Nc_H0 = Nc_acc_mtx*post_corr_ED_H0;
            ave_Nc_H1 = Nc_acc_mtx*post_corr_ED_H1;
            [peak_pow_H1, peakindex_H1] = max(ave_Nc_H1(1:STO_max));
            peak_pow_H0 = max(ave_Nc_H0(1:STO_max));
        end

    % ----- Single-Peak Detection with directional beams---------
    case 'directional'
        peak_pow_H1 = max(corr_out_H1);
        peak_pow_H0 = max(corr_out_H0);

    % ----- Single-Peak Detection with directional beams---------
    otherwise
        if STOinfo % Concept scenario where peak locations is know
            peak_pow_H1 = max(corr_out_H1(ZC_N+CP:burst_length:end));
            peak_pow_H0 = max(corr_out_H0(ZC_N+CP:burst_length:end));
            peakindex_H1 = 0;
        else % Practical scenario where peak location is unknown
            peak_pow_H1 = max(corr_out_H1_STO);
            peak_pow_H0 = max(corr_out_H0_STO);
            peakindex_H1 = 0;
%                     [peak_pow_H1(ss) peakindex_H1(ss)] = max(mean(reshape(corr_out_H1,burst_length,M+1),2));
%                     peak_pow_H0(ss) = max(mean(reshape(corr_out_H0,burst_length,M+1),2));

        end
end % end of switch


%% plot FSM-KW codebook
N_element_tx = Nt;
M_burst_tx = 64;%M_burst(1); % T/Rx partial, so it's smaller than 64
W_mat_tx = get_IA_BF(N_element_tx, M_burst_tx, 'PN'); % Tx beamformer in IA stage
angle_range = 90;
FF_tx = zeros(N_element_tx, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF_tx(:,kk) = exp(1j*pi*(0:N_element_tx-1).'*sin((kk - angle_range -1 )/180*pi));
end

N_element_rx = Nr;
M_burst_rx = 64;%M_burst(2); % T/Rx partial, so it's smaller than 64
W_mat_rx = get_IA_BF(N_element_rx, M_burst_rx, 'PN'); % Tx beamformer in IA stage
angle_range = 90;
FF_rx = zeros(N_element_rx, (angle_range*2+1));
for kk = 1:(angle_range*2+1)
    FF_rx(:,kk) = exp(1j*pi*(0:N_element_rx-1).'*sin((kk - angle_range -1 )/180*pi));
end
% for mm=1:M_burst
%     pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
%     % ------ beam pattern test -------
%     angle_test = 0;
%     % xdata = linspace(-pi/2,pi/2,181);
%     xdata = linspace(-pi/2,pi/2,181);
%     %
%     figure(100)
%     subplot(132)
%     polarplot(xdata+pi,pattern,'linewidth',2);hold on
%     set(gca,'FontSize', font_size)
%     title('FSM-Sec. BF [33]')
%     rlim([-20,20])
%     % thetalim([-90 90])
%     thetalim([90,270])
%     ax = gca;
%     % ax.ThetaTick = -90:22.5:90;
%     ax.ThetaTick = 90:22.5:270;
%     ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
%     grid on
% end



% %% plot PN codebook
% N_element = 256;
% M_burst = 16; % T/Rx partial, so it's smaller than 64
% W_mat = get_IA_BF(N_element, M_burst, 'PN'); % Tx beamformer in IA stage
% 
% angle_range = 90;
% FF = zeros(N_element, (angle_range*2+1));
% for kk = 1:(angle_range*2+1)
%     FF(:,kk) = exp(1j*pi*(0:N_element-1).'*sin((kk - angle_range -1 )/180*pi));
% end
% for mm=1:M_burst
%     pattern = 20*log10(abs(FF'*W_mat(:,mm))+1e-1);
%     % ------ beam pattern test -------
%     angle_test = 0;
%     % xdata = linspace(-pi/2,pi/2,181);
%     xdata = linspace(-pi/2,pi/2,181);
%     %
%     figure(100)
%     subplot(133)
%     polarplot(xdata+pi,pattern,'linewidth',2);hold on
%     set(gca,'FontSize', font_size)
%     title('PN BF (prop.)')
%     rlim([-20,20])
%     % thetalim([-90 90])
%     thetalim([90,270])
%     ax = gca;
%     % ax.ThetaTick = -90:22.5:90;
%     ax.ThetaTick = 90:22.5:270;
%     ax.ThetaTickLabels = {'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'};
%     grid on
% end

%% plot
blue = [   0    0.4470    0.7410];
figure

nImages = 64;
n_range = floor(linspace(0,length(Rx_sig_H1),nImages));

fig = figure;
for idx = 1:nImages
    
%     sector_index = idx;
    sector_index = floor((idx-1)/1)+1;
    sector_idx_tx = idx;%floor((sector_index-1)/M_burst_rx)+1;
    sector_idx_rx = idx;%sector_index - (sector_idx_tx-1)*M_burst_rx;
    % ---- subplot --------
    
    pattern_tx = 20*log10(abs(FF_tx'*W_mat_tx(:,sector_idx_tx))+1e-1);
    xdata = linspace(-pi/2,pi/2,181);
    subplot(221)
    polarplot(xdata+pi*0,pattern_tx,'linewidth',2);
    rlim([-20,20])
    thetalim([-90 90])
    title(['Tx Sector = ' num2str( sector_idx_tx) ])

%     thetalim([90,270])
    
    % ---- subplot --------
    pattern_rx = 20*log10(abs(FF_rx'*W_mat_rx(:,sector_idx_rx))+1e-1);
    xdata = linspace(-pi/2,pi/2,181);
    subplot(222)
    polarplot(xdata+pi,pattern_rx,'linewidth',2);
    rlim([-20,20])
    thetalim([90 270])
    title(['Rx Sector = ' num2str( sector_idx_rx) ])

%     thetalim([90,270])
    
    % ---- subplot --------
    data_range = 1:n_range(idx);
    subplot(223)
    plot((1:length(Rx_sig_H1(data_range)))/burst_N,abs(Rx_sig_H1(data_range)),'color',blue);hold on
    grid on
    xlabel('SS Burst Index')
    ylabel('Rx Signal')
    xlim([0,64])
    ylim([0,15])
    xticks(0:8:64)
    title(['Rx Signal till Burst Index = ' num2str( sector_index) ])

    
    data_range = ZC_N+CP + (1:n_range(idx));
    subplot(224)
    plot((1:length(corr_out_H1_STO(data_range)))/burst_N,abs(corr_out_H1_STO(data_range)),'color',blue);hold on
    grid on
    xlabel('SS Burst Index')
    ylabel('PSS Corr. Output')
    xlim([0,64])
    ylim([0,25])
    xticks(0:8:64)
    title(['PSS Corr. out till Burst Index = ' num2str( sector_index) ])
    
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
close;

% figure;
% for idx = 1:nImages
%     subplot(4,4,idx)
%     imshow(im{idx});
% end

filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.15);
    end
end



