clear;clc;
M = 20;
ZC_root = 29;
ZC_N = 255;
STO = randi(M*1*ZC_N);

seq = lteZadoffChuSeq(ZC_root,ZC_N);
% seq_0DC = ifft([seq(1:63);seq(65:end)]);
seq_1DC = ifft(seq);
Rx_1DC = [seq_1DC];
sig_length = length(Rx_1DC);
%%
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


