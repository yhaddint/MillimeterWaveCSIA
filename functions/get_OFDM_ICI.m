function [ sig_out ] = get_OFDM_ICI( sig_vec, PN_seq )
%GET_ICI_MATRIX Summary of this function goes here
%   This function provides ICI matrix defined as slides 01/24/2018
%   [ Q_mtx ] = get_ICI_matrix( PN_seq )
%   IP: sig_vec is Nfft*Nrf by 1 vector
%   IP: PN_seq is Nfft by Nrf matrix with time domain PN in Nrf RF-chain
%   OP: sig_out is output signal after modeling ICI from PN
    
    Nfft = size(PN_seq,1);
    Nrf = size(PN_seq,2);
    DFT = dftmtx(Nfft);
    sig_mtx = reshape(sig_vec, Nrf, Nfft);
    sig_out = zeros(Nrf*Nfft,1);
    
    for pp = 1:Nfft
        for nn=1:Nrf
            Q(nn,pp) = (PN_seq(:,nn).' * DFT(:,pp))/Nfft;
            % Each column has diagonal components of Q_p in slides
        end
    end
    
    Q_new = Q;
    for p1 = 1:Nfft
        sigout_index = (p1-1)*Nrf+1:p1*Nrf; % Index of M IP symbols at kth subcarriers
        sig_out(sigout_index) = sum(Q_new .* sig_mtx, 2);
        
        % Circulate o move into second block row of Q matrix in the slide
        Q_temp = Q_new;
        Q_new = [Q_temp(:,end), Q_temp(:,1:end-1)];
    end
end

