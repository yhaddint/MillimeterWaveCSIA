function [ sig_out ] = run_DA_WB_Combining( combiner_DA, sig_vec, M )
%FH_WB_COMBINING Summary of this function goes here
%   Combining approach in fully digital architecture.
%   [ sig_out ] = DA_WB_Combining( chan, sig_vec, M )
%   IP: combiner_DA is a Nr*M by Nfft generalized matrix for combiner matrix
%   IP: sig_vec is a Nfft*Nr by 1 vector; block divided by subcarriers
%   IP: M is scalar, multiplexing level
%   OP: sig_out is Nfft*M by 1 vector, the post-combining symbols

    Nr = size(combiner_DA,1)/M; % Number of Tx antennas 
    Nfft = size(combiner_DA,2); % Number of subcarriers
    sig_out = zeros(M*Nfft,1);
    
    for kk=1:Nfft
        sigin_index = (kk-1)*Nr+1:kk*Nr; % Index of Nr IP symbols at kth subcarriers
        sigout_index = (kk-1)*M+1:kk*M; % Index of M OP symbols at kth subcarriers

        sig_out(sigout_index) = (reshape(combiner_DA(:,kk),Nr,M))' * sig_vec(sigin_index);
    end
end

