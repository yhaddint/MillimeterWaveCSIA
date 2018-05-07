function [ y,x ] = get_emp_maxrv_cdf(mu_sig,sigma_sig,mu_noise,sigma_noise )
%GET_EMP_MAXRV_CDF Summary of this function goes here
%   Detailed explanation goes here
    MCtimes = 1e4;
    for MCindex = 1:MCtimes
        rv1 = randn * sigma_sig + mu_sig;
        rv2 = randn * sigma_noise + mu_noise;
        rv(MCindex) = max(rv1,rv2);
    end
    [y,x] = ecdf(rv);

end

