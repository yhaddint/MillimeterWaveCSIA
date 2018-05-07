clear;clc;

% ------ script control parameter -------
rng(2)
STOinfo = 0; % Whether system know STO perfectly
MCtimes = 1e2; % Monte Carlo simulation

% ----- system parameters --------
M = 20; % Number of ZC burst
ZC_root = 29; % ZC root, a prime value with ZC length
ZC_N = 255; % ZC sequence length
burst_N = ZC_N*2; % number of sample in each burst (2OFDM symbol for now)
noise_pow = 1e-1;
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

%%
% ------------ MC iterations --------------
for MCindex = 1:MCtimes
    
    % ------- Random parameter realizations -----
    STO = 0;%randi(M*sig_length); % Sample timing offset
    
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
    for nn=1:Tx_sig_length
        precoder_index = floor( (nn-1) / burst_length )+1;
        combiner_index_raw = floor( (nn + STO - 1) / burst_length )+1;
        combiner_index = mod(combiner_index_raw-1,M)+1;
        w_vec = W_mat(:,combiner_index);
        v_vec = V_mat(:,precoder_index);
        Rx_sig(nn) = (w_vec'*H_chan*v_vec) * Tx_sig(nn);
    end
    
    % ------- AWGN -------
    noise = (randn(Tx_sig_length,1)+1j*randn(Tx_sig_length,1))/sqrt(2)*sqrt(noise_pow);
    Rx_sig_H1 = Rx_sig + noise;
    Rx_sig_H0 = noise;
    
    % ------ T Domain ZC Correlation -------
    corr_out_H1(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H1)).^2; % corr rx t-domain sig with ZC
    corr_out_H0(1:Rx_sig_length) = abs(conv(ZC_t_domain,Rx_sig_H0)).^2; % corr rx t-domain sig with ZC
    
    % ----- Multi-Peak Detection ---------
    if STOinfo % Concept scenario where peak locations is know
        peak_pow_H1(MCindex) = sum(corr_out_H1(ZC_N:burst_length:end));
        peak_pow_H0(MCindex) = sum(corr_out_H0(ZC_N:burst_length:end));
    else % Practical scenario where peak location is unknown
        peak_pow_H1(MCindex) = max(sum(reshape(corr_out_H1,burst_length,M+1),2));
        peak_pow_H0(MCindex) = max(sum(reshape(corr_out_H0,burst_length,M+1),2));
    end
end
%%
figure
[a,b] = ksdensity(peak_pow_H1);hold on;
plot(b,a);grid on
[a,b] = ksdensity(peak_pow_H0);hold on;
plot(b,a);grid on
legend('H1','H0')

%% The following is used for theoretical analysis
env = exp(1j*rand(3,1)*2*pi);
figure
subplot(211)
plot(abs(conv(circshift([Rx_1DC*env(1);Rx_1DC*env(2);Rx_1DC*env(3)],20),conj(flipud(seq_1DC)))./ZC_N));
hold on
grid on
title('With STO')
subplot(212)
plot(abs(conv([Rx_1DC*env(1);Rx_1DC*env(2);Rx_1DC*env(3)],conj(flipud(seq_1DC)))./ZC_N));
grid on
title('No STO')

%%
rx_sig = repmat(Rx_1DC,M,1);
tx_bf_gain = zeros(M,1);
rx_bf_gain = zeros(M,1);
for mm=1:M
    tx_bf_gain(mm) = randn+1j*randn;
    rx_bf_gain(mm) = randn+1j*randn;
end

tx_bf_gain_vec = kron(tx_bf_gain,ones(sig_length,1));
rx_bf_gain_vec = kron(rx_bf_gain,ones(sig_length,1));
rx_bf_gain_vec_STO = circshift(rx_bf_gain_vec,STO);

rx_sig_STO = rx_sig.*tx_bf_gain_vec.*rx_bf_gain_vec_STO;
rx_sig_noSTO = rx_sig.*tx_bf_gain_vec.*rx_bf_gain_vec;
%%
figure
subplot(211)
plot(abs(conv(rx_sig_noSTO,conj(flipud(seq_1DC)))./ZC_N));
grid on
title('No Timing Offset')
subplot(212)
plot(abs(conv(rx_sig_STO,conj(flipud(seq_1DC)))./ZC_N));
grid on
title('With Timing Offset')
%% theoretical power
MCtimes = 5e4;
M = 10;
Nr = 32;
for MCindex = 1:MCtimes
    STO_fraction = rand;
    phi = rand * (pi/3) - (pi*2/3);
    arx = exp(1j*pi*(0:Nr-1).'*sin(phi));
    alpha = arx'*exp(1j*rand(Nr,2)*2*pi)/Nr;
    corr_peak_STO(MCindex) = abs(STO_fraction * alpha(1) + (1-STO_fraction)*alpha(2))^2;
    corr_peak(MCindex) = abs(alpha(1))^2;
end
figure
[a,b] = ecdf(10*log10(mean(reshape(corr_peak,M,MCtimes/M),1)));
plot(b,a);hold on
[a,b] = ecdf(10*log10(mean(reshape(corr_peak_STO,M,MCtimes/M),1)));
plot(b,a);hold on
legend('No STO','With STO')
grid on
%%
M = 20;
MCtimes = 1e3;
for MCindex = 1:MCtimes
    alphai = randn(M,1)+1j*randn(M,1);
    alphaj = randn(M,1)+1j*randn(M,1);
    theta = rand;
    pow_STO(MCindex) = sum(abs(alphai*theta + alphaj*(1-theta)).^2);
    pow_noSTO(MCindex) = sum(abs(alphai).^2);
end
figure
[a,b] = ecdf(pow_STO);hold on;
plot(b,a);grid on
[a,b] = ecdf(pow_noSTO);hold on;
plot(b,a);grid on
legend('STO','No STO')


