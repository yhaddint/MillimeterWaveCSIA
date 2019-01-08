function [CS_TH, sec_TH] = get_Genie_Detection_TH()
%GET_GENIE_DETECTION_TH Summary of this function goes here
%   Detailed explanation goes here

BW = 57.6e6;
NF = 4;
NB = 512;
M = 64;
Tx_sig_length = NB*M;
noise_pow = 10^( (-174 + 10*log10(BW) + NF - 30)/10 );
MCtimes = 5e3;
PFA = 0.01;
CP = 32;
Nc = 31;
burst_length = 512;
STO_max = 470;

% ZC sequence
ZC_root = 29;                               % ZC root, a coprime number with ZC length
ZC_N = 127;                                 % ZC sequence length
seq = lteZadoffChuSeq(ZC_root,ZC_N);        % Generate ZC sequence
seq_1DC = ifft(seq)*sqrt(ZC_N+1);           % T-domain sig. for ZC corr in detction & STO estimation; DC subcarrier is not null
ZC_t_domain = conj(flipud(seq_1DC));        % ZC sequence used for correlation in Rx


Nc_acc_mtx = toeplitz([1,zeros(1,burst_length-Nc)]',[ones(1,Nc),zeros(1,burst_length-Nc)]);

for MCidx = 1:MCtimes
    
    STO = randi(50)+50;
    noise_at_STO = (randn(STO,1)+1j*randn(STO,1))/sqrt(2);
    noise_CP = (randn(Tx_sig_length,1)+1j*randn(Tx_sig_length,1))/sqrt(2);

    awgn_CP = noise_CP * sqrt(noise_pow);

    % case 'sector beamformer'
    Rx_sig_H0_sec = awgn_CP;

    % case 'PN beamformer'
    Rx_sig_H0 = awgn_CP;

    % ------ T Domain ZC Correlation -------

    Rx_sig_H0_wSTO = [noise_at_STO * sqrt(noise_pow); Rx_sig_H0];
    corr_out_H0_STO = abs(conv(ZC_t_domain,Rx_sig_H0_wSTO)/ZC_N).^2; % corr rx t-domain sig with ZC

    Rx_sig_H0_wSTO_sec = [noise_at_STO * sqrt(noise_pow); Rx_sig_H0_sec];
    corr_out_H0_STO_sec = abs(conv(ZC_t_domain,Rx_sig_H0_wSTO_sec)/ZC_N).^2; % corr rx t-domain sig with ZC

    % ----- Multi-Peak Detection w. PN beam ---------
    % Practical scenario where peak location is unknown

    post_corr_ED_H0 = sum(reshape([corr_out_H0_STO(ZC_N+CP:end);...
        zeros(burst_length*(M+2)-length(corr_out_H0_STO(ZC_N+CP:end)),1)],...    
        burst_length,M+2),2)/M;

    % Moving average of Nc samples to include energy of all paths
    ave_Nc_H0 = Nc_acc_mtx * post_corr_ED_H0;

    peak_pow_H0(MCidx) = max(ave_Nc_H0(1:STO_max));

    % ----- Single-Peak Detection w. directional beams---------
    % Practical scenario where peak location is unknown
    peak_pow_H0_sec(MCidx) = max(corr_out_H0_STO_sec);
end

CS_H0 = sort(peak_pow_H0,'descend');
CS_TH = CS_H0(MCtimes*PFA);

sec_H0 = sort(peak_pow_H0_sec,'descend');
sec_TH = sec_H0(MCtimes*PFA);



end

