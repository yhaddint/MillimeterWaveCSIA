%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
path_num = 4; % Num of rays in a cluster
Nr = 16; % Number of antenna in Rx
Nt = 32;
M_range = [34:2:44,48:4:88,94:6:118,126:15:180];
% M_range = [32,48,64];
% M = 64; % Length of training
MCtimes = 1e3; % Num of Monte Carlo Sim.
AOAspread2 = 0;
AOAspread = 0;
AODspread2 = 0;
AODspread = 0;
SNR_num = 9;
SNR_data_range = linspace(-40,-20,SNR_num);
SNR_range = SNR_data_range + 30; % There is 1000x diff. in IA and data AWGN BW
P = 1; % subcarrier. it is one in NB case

%-------- dictionary generation -------------
cand_num_r = 41;
cand_num_t = 61;
dict_num = cand_num_r*cand_num_t;

cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Nt-1)'*pi*sin(cand_angle_t));

AOA_error = zeros(MCtimes,SNR_num,length(M_range));
AOD_error = zeros(MCtimes,SNR_num,length(M_range));
AoA_align_i = zeros(MCtimes,SNR_num,length(M_range));
AoD_align_i = zeros(MCtimes,SNR_num,length(M_range));

desired_AoA_gap = (105/180*pi)/(Nr/2);
desired_AoD_gap = (105/180*pi)/(Nt/2);

% zero initialization of true AoA/AoD
phi = zeros(MCtimes, path_num);
theta = zeros(MCtimes, path_num);


% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCindex = 1:MCtimes
    clc
    fprintf('iteration %d:\n',MCindex);
        % AoA of rays with disired seperation
          
    % complex gain with normalization; Equi to exp(1j*rand*2*pi) in L=1 case
    g_ray0 = randn(path_num,1) + 1j*randn(path_num,1);
    g_ray = g_ray0./norm(g_ray0);

    % generate AoA/AoD with desired gap
    phi_temp = rand(1,path_num)*2*pi/3-pi/3;
    while get_gap(phi_temp)<desired_AoA_gap
        phi_temp = rand(1,path_num)*2*pi/3-pi/3;
    end
    
    theta_temp = rand(1,path_num)*2*pi/3-pi/3;
    while get_gap(theta_temp)<desired_AoD_gap
        theta_temp = rand(1,path_num)*2*pi/3-pi/3;
    end
    
    % AoA of paths with disired seperation
    phi(MCindex, :) = phi_temp;
    
    % AoD of paths with disired seperation
    theta(MCindex, :) = theta_temp;
    
    
    for pathindex = 1:path_num
        % Spatial response and its derivative over phi
        arx(:,pathindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(MCindex, pathindex)));%/sqrt(Nr);
        
        % Spatial response and its derivative over theta
        atx(:,pathindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(MCindex, pathindex)));%/sqrt(Nt);   
    end
    
    % Processing. Various measurement (M) are tested
    for M_index=1:length(M_range)
        
        % Pick M from pool
        M=M_range(M_index);

        % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
        % directional beam from angle steering vector
        probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
        W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);

        probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
        F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(M);
        
        % Pre-compute dictionary
        Measure_mat = kron(transpose(F)*conj(cand_ARV_t),W'*cand_ARV_r);
        select_row = zeros(1,M);
        for ii=1:M
            select_row(ii) = (ii-1)*M+ii;
        end
        Measure_mat_new = Measure_mat(select_row,:);
        for cc=1:cand_num_r*cand_num_t
            Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
        end
        
        % CS channel measurement (noiseless)
        sig_rx = zeros(P*M,1);
        
        for mm=1:M
            sig_rx(mm) = 0;
            for ll=1:path_num
                sig_rx(mm) = sig_rx(mm) + g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll));
            end
        end
        
        % AWGN realization (normalized power)
        awgn = (randn(P*M,1)+1j*randn(P*M,1))/sqrt(2);
        
        % For loop for SNR
        for ss = 1:SNR_num

            % SNR and additive awgn (power scaling)
            sigman2 = 10^(-SNR_range(ss)/10);
            sig_noisy = sig_rx + awgn * sqrt(sigman2);

            % Finding angle
            sig_desymb = zeros(M,1);
            for mm=1:M
                indextemp = (mm-1)*P+1:mm*P;
                sig_desymb(mm) = sig_noisy(indextemp);
            end

            for dd=1:dict_num
                score_final(dd) = abs(sig_desymb'*(Measure_mat_new(:,dd)))...
                    /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
            end
            [~,bestindex(MCindex)] = max(abs(score_final));
            bestrow = floor((bestindex(MCindex)-1)/cand_num_r)+1;
            bestcol = bestindex(MCindex)-(bestrow-1)*cand_num_r;
            bestAOA(MCindex,ss) = (bestcol-1)*AOAstep-60*pi/180;
            bestAOD(MCindex,ss) = (bestrow-1)*AODstep-60*pi/180;
            if path_num==1
                AOA_error(MCindex,ss,M_index) = abs(bestAOA(MCindex,ss) - phi(MCindex,1));
                AOD_error(MCindex,ss,M_index) = abs(bestAOD(MCindex,ss) - theta(MCindex,1));
            else
                [AOA_error(MCindex,ss,M_index), AoA_align_i(MCindex,ss,M_index)] =...
                            min(abs(bestAOA(MCindex,ss) - phi(MCindex,:)));
                [AOD_error(MCindex,ss,M_index), AoD_align_i(MCindex,ss,M_index)] =...
                            min(abs(bestAOD(MCindex,ss) - theta(MCindex,:)));
            end
        end 
    end
end
%%
for ss=1:SNR_num
    for M_index = 1:length(M_range)
    AOAalign_mean(ss,M_index) = sum((squeeze(AOA_error(:,ss,M_index))/pi*180)<(105/Nr),1)/MCtimes;
    AODalign_mean(ss,M_index) = sum((squeeze(AOD_error(:,ss,M_index))/pi*180)<(105/Nt),1)/MCtimes;
    end
end

for ss=1:SNR_num
    for M_index = 1:length(M_range)
        DoubleAlign(ss,M_index) = sum(((squeeze(AOA_error(:,ss,M_index))/pi*180)<(105/Nr)).*...
                      ((squeeze(AOA_error(:,ss,M_index))/pi*180)<(105/Nr)).*...
                      (squeeze(AoA_align_i(:,ss,M_index))==squeeze(AoD_align_i(:,ss,M_index))))/MCtimes;
    end
end

figure
subplot(211)
plot(SNR_range,AOAalign_mean,'-o','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
title('AoA Alignment')
subplot(212)
plot(SNR_range,AODalign_mean,'--o','linewidth',2);hold on
grid on
title('AoD Alignment')
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legendtext = [];
for M_index = 1:length(M_range)
    legendtext = [legendtext;'M=',num2str(M_range(M_index))];
end
legend(legendtext)

figure
plot(SNR_range,DoubleAlign,'-o','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legendtext = [];
for M_index = 1:length(M_range)
    legendtext = [legendtext;'M=',num2str(M_range(M_index))];
end
legend(legendtext)

%% RMSE evaluation
portion = 0.90;
for ss=1:SNR_num
    
    input_sig = abs(AOD_error(:,ss)).^2;
    RMSE_theta(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
    
    input_sig = abs(AOA_error(:,ss)).^2;
    RMSE_phi(ss) = sqrt(mean(get_portion(input_sig,portion)))*(1/pi*180);
end
figure
subplot(311)
semilogy(SNR_range,RMSE_delay);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of normalized delay [ns]')
title('delay est')

subplot(312)
semilogy(SNR_range,RMSE_CFO);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of CFO est [Hz]')
title('CFO est')

subplot(313)
semilogy(SNR_range,RMSE_theta);hold on
semilogy(SNR_range,RMSE_phi);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA/AoD [deg]')
title('angle est')

