% script to simulate and evaluate probability of n* = max (n_i) where n_i
% is i.i.d normal r.v.
% There is a useful reference "Extremes" by Rober Wolpert, Duke Univ. 

clear;clc;
mu = 10;
sigma = 1;
N = 500;
MCtimes = 5e4;

% --- simulation results ------
for MCindex = 1:MCtimes

    result_max(MCindex) = max(randn(N,1)*sigma+mu);
%     result(MCindex) = abs(randn(1,1)).^2;
end

figure
% [a,b] = ksdensity(result);
% plot(b,a);hold on
% [a,b] = ksdensity(result_max);
[a,b] = ksdensity(result_max);
plot(b,a);hold on
grid on


% The theoretical value from reference is
% mu_max = mu - sigma * invnormcdf(1/N)
% sigma_max = -sigma / invnormcdf(1/N)
mu_max = mu - sigma*(-qfuncinv(1/N));
sigma_max = -sigma/(-qfuncinv(1/N));
x = linspace(mu_max-4*sigma_max,mu_max+4*sigma_max,1e3);
y = evpdf(-x,-mu_max,sigma_max);
plot(x,y)
legend('sim.','approx.')