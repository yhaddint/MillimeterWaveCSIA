%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
path_num = 1; % Num of rays in a cluster
Nr = 8; % Number of antenna in Rx
Nt = 64;
M = 64; % Length of training
MCtimes = 5e1; % Num of Monte Carlo Sim.
AOAspread2 = 0;
AOAspread = 0;
AODspread2 = 0;
AODspread = 0;
SNR_num = 5;
SNR_range = linspace(-20,20,SNR_num);
Ts = 1/(50e6);
Nb = 512;
CFO_ppm = 1; % CFO in ppm
CFO = rand*(28e9/1e6*CFO_ppm); % With unit Hz
eF = CFO*Ts*2*pi; % 
P = 128;
DFT = dftmtx(P);

%-------- dictionary generation -------------
cand_num_r = 61;
cand_num_t = 121;
dict_num = cand_num_r*cand_num_t;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Nt-1)'*pi*sin(cand_angle_t));

% test scenario when assuming CFO is known
phase_error_mat = kron(exp(1j*eF*Nb*(0:M-1)),exp(1j*eF*(0:P-1)'));
phase_error = reshape(phase_error_mat,M*P,1);

% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCindex = 1:MCtimes
    
    clc
    fprintf('iteration %d:\n',MCindex);
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);
    
    Measure_mat = kron(transpose(F)*conj(cand_ARV_t),W'*cand_ARV_r);
    select_row = zeros(1,M);
    for ii=1:M
        select_row(ii) = (ii-1)*M+ii;
    end
    Measure_mat_new = Measure_mat(select_row,:);
    for cc=1:cand_num_r*cand_num_t
        Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
    end

    
%     probe_Tx_BF = ones(Nt,M);
%     F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);   

    % AoA of rays with disired seperation
    phi = zeros(path_num,1);
    phi0 = 0/180*pi;
    phi = phi0 + randn(path_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta = zeros(path_num,1);
    theta0 = 0/180*pi;
    theta = theta0 + randn(path_num,1) * AODspread;

%     % Gain
%     g_cmplx = exp(1j*rand(ray_num,1)*2*pi)/sqrt(ray_num);
%     g = g_cmplx;
    % Rotate of ray
    tau = rand*(90e-9);
    tau_samp(MCindex) = tau/Ts*2*pi;
%     g_ray = (randn+1j*randn)/sqrt(2);
    g_ray = exp(1j*rand*2*pi);

    % Pre-compute some vectors/matrices in FIM
    for pathindex = 1:path_num

        % Spatial response and its derivative over phi
        arx(:,pathindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(pathindex)))/sqrt(Nr);
        Darx(:,pathindex) = 1j * pi * (0:Nr-1)' * cos(phi(pathindex));
        drx(:,pathindex) = Darx(:,pathindex).*arx(:,pathindex);

        % Spatial response and its derivative over theta
        atx(:,pathindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(pathindex)))/sqrt(Nt);
        Datx(:,pathindex) = 1j * pi * (0:Nt-1)' * cos(theta(pathindex));
        dtx(:,pathindex) = Datx(:,pathindex).*atx(:,pathindex);
        
        % Delay response and its derivative over tau
        fvec(:,pathindex) = exp(-1j * (0:P-1)' * tau_samp(MCindex) / P);
        Dfvec(:,pathindex) = -1j * (0:P-1)'/P;
        dfvec(:,pathindex) = Dfvec(:,pathindex).*fvec(:,pathindex);
    end
    
    % About CFO and its derivative
    qvec = exp(1j * (0:P-1)' * eF);
    for mm=1:M
        timeindex = ((mm-1)*Nb+0):((mm-1)*Nb+P-1);
        Dqvec(:,mm) = 1j * (timeindex).';
        dqvec(:,mm) = exp(1j*Nb*(mm-1)*eF)*qvec.*Dqvec(:,mm);
    end
    
    
    % Pre-compute some equations, vecx, vecalpha, vecbeta
    % vecx is derivative of phi_0 in f(r,eta)
    vaa = zeros(M,1);
    vad = zeros(M,1);
    vda = zeros(M,1);
    
    % Zero initialization of vectors for FIM
    vdcfo = zeros(P*M,1);
    vdtheta = zeros(P*M,1);
    vdphi = zeros(P*M,1);
    vdtau = zeros(P*M,1);
    vdalpha = zeros(P*M,1);
    vdbeta = zeros(P*M,1);
    
    % Received signals
    symb = exp(1j*rand(P,1)*2*pi);
    tau_num = 200;
    delay_cand = linspace(0,100,tau_num)*1e-9/Ts*2*pi;
    for tt=1:tau_num
        delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P).*symb);
    end
    
    sig_rx = zeros(P*M,1);
    for ll=1:path_num
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            sig_rx(index) = g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*eF);
        end
    end
    
    awgn = (randn(P*M,1)+1j*randn(P*M,1))/sqrt(2);
    % For loop for SNR

    for ss = 1:SNR_num
        
        % SNR and addint awgn
        sigman2 = 10^(-SNR_range(ss)/10);
        sig_noisy = sig_rx + awgn * sqrt(sigman2);

        % CFO estimation
        sig_ave = mean(reshape(sig_noisy,P,M),2);
        for tt=1:tau_num
            if MCindex==15 & tt==6
                apple=1;
            end
            sig_tone = sig_ave.*conj(delay_mtx(:,tt));
            for pp = 1:P-1
                CFO_hat(pp) = angle(sig_tone(pp+1).*conj(sig_tone(pp)));
            end
            CFO_est = mean(CFO_hat);
            CFO_est1 = angle(mean(exp(1j*CFO_hat)));
%             CFO_est = eF;
            sig_deCFO = sig_ave.*exp(-1j*CFO_est*(0:P-1).');
            score(tt) = abs(sum(sig_deCFO.*conj(delay_mtx(:,tt)))/norm(delay_mtx(:,tt))^2);
        end
        [maxscore,maxindex] = max(score);
        maxindex_est(ss,MCindex) = maxindex;
        delay_est(ss,MCindex) = delay_cand(maxindex);
        
        % Finding angle
        

        for dd=1:dict_num
            if dd==3691
                apple=1;
            end
            sig_cand = kron(Measure_mat_new(:,dd),delay_mtx(:,maxindex));
            score_final(dd) = abs(sig_noisy'*(sig_cand.* phase_error))/norm(sig_cand)^2;
        end
        [~,bestindex_comp(MCindex)] = max(abs(score_final));
        bestrow = floor((bestindex_comp(MCindex)-1)/cand_num_r)+1;
        bestcol = bestindex_comp(MCindex)-(bestrow-1)*cand_num_r;
        bestAOA_nocomp(MCindex,ss) = (bestcol-1)*AOAstep-60*pi/180;
        bestAOD_nocomp(MCindex,ss) = (bestrow-1)*AODstep-60*pi/180;
        
        AOA_error_nocomp(MCindex,ss) = abs(bestAOA_nocomp(MCindex,ss) - phi0);
        AOD_error_nocomp(MCindex,ss) = abs(bestAOD_nocomp(MCindex,ss) - theta0);
%         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);

        % RMSE evaluation from CRLB perspective
%         CRLB_CFO(ss,MCindex) = sqrt(temp(1,1))*(1/Ts/2/pi);
%         CRLB_theta(ss,MCindex) = sqrt(temp(2,2))*(1/pi*180);
%         CRLB_phi(ss,MCindex) = sqrt(temp(3,3))*(1/pi*180);
    end
end
%%
for ss=1:SNR_num
AOAalign_comp_mean(ss) = sum((AOA_error_nocomp(:,ss)/pi*180)<(105/Nr),1)/MCtimes;
AODalign_comp_mean(ss) = sum((AOD_error_nocomp(:,ss)/pi*180)<(105/Nt),1)/MCtimes;
end

figure
plot(SNR_range,AOAalign_comp_mean,'--','linewidth',2);hold on
plot(SNR_range,AODalign_comp_mean,'--','linewidth',2);hold on

grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
%% RMSE evaluation
for ss=1:SNR_num
    RMSE_delay(ss) = sqrt(mean(abs(delay_est(ss,:)-tau_samp).^2))*(Ts/(2*pi)/1e-9);
end
figure
plot(SNR_range,RMSE_delay);hold on
grid on
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of normalized delay [ns]')


%% Plot CRLB of angle est./tracking
% figure
% subplot(211)
% semilogy(SNR_range,mean(CRLB_theta,2));hold on
% semilogy(SNR_range,mean(CRLB_phi,2));hold on
% grid on
% legend('Theta','Phi')
% xlabel('Point-to-Point SNR [dB]')
% ylabel('RMSE of AoA/AoD [deg]')
% 
% subplot(212)
% semilogy(SNR_range,mean(CRLB_CFO,2));hold on
% grid on
% title('CRLB of CFO')
% xlabel('Point-to-Point SNR [dB]')
% ylabel('RMSE of CFO [Hz]')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
%% best should be 3691
for mm=1:10
    temp(mm) = (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll));
end
[Measure_mat_new(1:10,3691) 22.6175*temp.']
