function [ sig_out ] = run_DA_WB_Precoding( precoder_DA, sig_vec, M )
%DA_WB_PRECODING Summary of this function goes here
%   Precoding approach in fully digital architecture.
%   This version does not include water-filling power allocation
%   [ sig_out ] = DA_WB_Precoding( chan, sig_vec, M )
%   IP: chan is a Nr by Nt by Nfft generalized matrix for channel matrix
%   IP: sig_vec is a Nfft*M by 1 vector; block divided by subcarriers
%   IP: M is scalar, multiplexing level
%   OP: sig_out is Nfft*Nt by 1 vector, the precoded symbols

    Nt = size(precoder_DA,1)/M; % Number of Tx antennas 
    Nfft = size(precoder_DA,2); % Number of subcarriers
    sig_out = zeros(Nfft*Nt,1);
    
    for kk=1:Nfft
        sigin_index = (kk-1)*M+1:kk*M; % Index of M IP symbols at kth subcarriers
        sigout_index = (kk-1)*Nt+1:kk*Nt; % Index of Nt OP symbols at kth subcarriers       
        sig_out(sigout_index) = reshape(precoder_DA(:,kk),Nt,M) * sig_vec(sigin_index);
    end
end

