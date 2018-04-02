% Script to test get_ICI_matrix 

Nfft = 16;
Nrf = 2;
sig_vec = ones(Nfft*Nrf,1);
PN_seq = ones(Nfft, Nrf);
[ sig_out ] = get_OFDM_ICI( sig_vec, PN_seq );