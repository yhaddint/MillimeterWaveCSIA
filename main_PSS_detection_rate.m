clear;clc;

% ------ script control parameter -------
rng(2)
MCtimes = 1e2; % Monte Carlo simulation
Pfa = 0.05;
SNR_data_range = -25:5:10;
SNR_range = SNR_data_range + 15;
SNR_num = length(SNR_range);

%%
for ss = 1:SNR_num
    SNR_dB = SNR_range(ss);
    STO = 0;
    % ------------ MC iterations --------------
    for MCindex = 1:MCtimes
        clc
        fprintf('SNR =  %d dB\n',SNR_dB);
        fprintf('Iteration %d:\n',MCindex);
        [ peak_pow_H1(MCindex), peak_pow_H0(MCindex) ] = run_PSS_detection( SNR_dB, STO );
    end

    % --------- Detection based on emprical threshold --------
    temp_H0 = sort(peak_pow_H0,'descend');
    TH = temp_H0(MCtimes*Pfa);
    Pm(ss) = sum(peak_pow_H1<TH)/MCtimes;
end

%% Figure
figure
semilogy(SNR_data_range,Pm);hold on
grid on
xlabel('SNR (dB)')
ylabel('Miss Detection of PSS')





