clear;clc;
N = 512;
DFT = dftmtx(N);
PN_sigma2 = 1e-1;
MCtimes = 1000;
for MCindex = 1:MCtimes
PN_seq = zeros(N,1);
PN_seq(1) = 1;
for nn=1:N-1
    PN_seq(nn+1) = PN_seq(nn)*exp(1j*randn*sqrt(PN_sigma2));
end
%%
result(:,MCindex) = DFT'* (PN_seq.*DFT(:,1))/N;
end

figure
semilogy(mean(abs(result),2));
