clear;clc;

% parameters
N_range = [100,200,400,800];

% Original Gaussian parameter
mu0 = 5;
sigma0 = 1;
MCtimes = 1e5;
Pfa = 0.05;

% color
plotcolor = [0    0.4470    0.7410;
            0.8500    0.3250    0.0980;
            0.9290    0.6940    0.1250;
            0.4940    0.1840    0.5560];

%% MC sim - emprical dist.
for nn=1:length(N_range)
    
    % pick parameter N
    N = N_range(nn);
    
    % Pick maximum from N Gaussian R.V.
    for MCindex=1:MCtimes
        emp(MCindex) = max(randn(N,1)*sigma0+mu0);
    end
    
    % emprical cdf
    [b1(:,nn),a1(:,nn)] = ecdf(emp);  
    
end

%% Gumbel distribution - theoretical dist.
for nn=1:length(N_range)
    
    % pick parameter N
    N = N_range(nn);
    
    % parameter for Gumbel Dist.
    mu = mu0 + sigma0*qfuncinv(1/N);
    sigma = sigma0/qfuncinv(1/N);
    beta = sigma*sqrt(6)/pi;
    
    % theoretical cdf
    xdata(:,nn) = linspace(0,max(emp),100);
    ydata(:,nn) = exp(-exp(-(xdata(:,nn)-mu)/beta));
    TH(nn) = mu-beta*log(-log(1-Pfa));
end

%% plots
figure(1)
for nn=1:length(N_range)
    semilogy(a1(:,nn),1-b1(:,nn),'--','linewidth',2,'color',plotcolor(nn,:));hold on
    semilogy(xdata(:,nn),1-ydata(:,nn),'-','linewidth',2,'color',plotcolor(nn,:));hold on
end
grid on
ylim([0.001,1])
xlim([min(emp),max(emp)])
xlabel('X')
ylabel('CCDF P(x>X)')
yticks([0.001,0.01,0.1,1])
yticklabels({'0.001','0.01','0.1','1'})
% legend('N=200 (Sim.)','N=200 (Theo.)')

legendtext = [];
for nn = 1:2*length(N_range)
    nindex = ceil(nn/2);
    if mod(nn,2)==1
        legendtext = [legendtext;'N=',num2str(N_range(nindex)),' (Simu.)'];
    else
        legendtext = [legendtext;'N=',num2str(N_range(nindex)),' (Appx.)'];
    end
end
legend(legendtext)
