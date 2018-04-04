clear;clc;

% ------ script control parameters -------
rng(3)
MCtimes = 5e2; % Monte Carlo simulation
Pfa = 0.05;
SNR_data_range = -25:5:0;
SNR_range = SNR_data_range + 15;
SNR_num = length(SNR_range);

%% Monte Carlo evaluations

STO = 'random'; % Type of sample timing offset
BFtype = 'PN'; % Type of beamformer in IA
STOinfo = 1; % Assuming perfect knowledge of peak
M_burst = [256, 256]; % Number of bursts in IA; For directional use [M_Tx_BF,M_Rx_BF] for beams in BS and UE

% ------------ MC iterations (each has all SNRs)--------------
for MCindex = 1:MCtimes
    clc
    fprintf('Iteration %d:\n',MCindex);
    [ peak_pow_H1(:,MCindex),...
      peak_pow_H0(:,MCindex) ] = run_PSS_detection( SNR_range,...
                                                    STO,...
                                                    STOinfo,...
                                                    BFtype,...
                                                    M_burst);
end

% --------- Detection based on emprical threshold --------
for ss = 1:SNR_num
    temp_H0 = sort(peak_pow_H0(ss,:),'descend');
    TH(ss) = temp_H0(MCtimes*Pfa);
    Pm(ss) = sum(peak_pow_H1(ss,:)<TH(ss))/MCtimes;
end


%% Figure
figure
semilogy(SNR_data_range,Pm);hold on
grid on
xlabel('SNR (dB)')
ylabel('Miss Detection of PSS')
%%

% figure
% [a,b] = ksdensity(peak_pow_H1(1,:));
% plot(b,a);hold on
% [a,b] = ksdensity(peak_pow_H0(1,:));
% plot(b,a);hold on
% grid on



