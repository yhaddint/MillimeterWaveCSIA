function [ sig_out ] = get_rapp_NL( sig_in,Vsat,P )
%GET_RAPP_NL Summary of this function goes here
%   Detailed explanation goes here
    
    sig_ratio = abs(sig_in)/Vsat;
    sig_out = sig_in./((1+sig_ratio.^(2*P)).^(1/2/P));

end

