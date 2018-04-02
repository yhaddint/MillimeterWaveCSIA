function [ precoder_DA, combiner_DA ] = get_WB_DA_beamformer( chan, M )
%GET_WB_DA_BEAMFORMER Summary of this function goes here
%   Precoder and combiner in fully digital architecture.
%   This version does not include water-filling power allocation
%   [ precoder_DA, combiner_DA ] = DA_WB_Precoding( chan, M )
%   IP: chan is a Nr by Nt by Nfft generalized matrix for channel matrix
%   IP: M is scalar, multiplexing level
%   OP: precoder_DA is M*Nt by Nfft matrix, the structured precoding matrix
%   OP: combiner_DA is M*Nr by Nfft matrix, the structured precoding matrix

    Nr = size(chan,1); % Number of Tx antennas 
    Nt = size(chan,2); % Number of Tx antennas 
    Nfft = size(chan,3); % Number of subcarriers
    precoder_DA = zeros(M*Nt, Nfft);
    combiner_DA = zeros(M*Nr, Nfft);
    
    for kk=1:Nfft
        [U, Sigma, V] = svd(squeeze(chan(:,:,kk))); % Standard SVD approach
        precoder_DA(:,kk) = reshape(V(:,1:M),M*Nt,1);
        combiner_DA(:,kk) = reshape(U(:,1:M),M*Nr,1);

    end
    
end

