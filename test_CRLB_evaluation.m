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
M = 32; % Length of training
MCtimes = 50; % Num of Monte Carlo Sim.
AOAspread2 = 0;
AOAspread = 0;
AODspread2 = 0;
AODspread = 0;
SNR_num = 50;
SNR_range = linspace(-35,15,SNR_num);
Ts = 1/(50e6);
Nb = 512;
CFO_ppm = 1; % CFO in ppm
CFO = 28e9/1e6*CFO_ppm; % With unit Hz
eF = CFO*Ts*2*pi; % 
P = 128;
DFT = dftmtx(P);

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
    tau = rand * (200e-9);
    tau_samp = mod(tau/Ts*2*pi,2*pi);
    g_ray = (randn+1j*randn)/sqrt(2);

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
        fvec(:,pathindex) = exp(1j * (0:P-1)' * tau_samp);
        Dfvec(:,pathindex) = 1j * (0:P-1)';
        dfvec(:,pathindex) = Dfvec(:,pathindex).*fvec(:,pathindex);
    end
    
    % About CFO and its derivative
    qvec = exp(1j * (0:P-1)' * eF);
    for mm=1:M
        timeindex = ((mm-1)*Nb+0):((mm-1)*Nb+P-1);
        Dqvec(:,mm) = 1j * (timeindex).';
        dqvec(:,mm) = exp(1j*Nb*(mm-1)*CFO)*qvec.*Dqvec(:,mm);
    end
    
    
    % Evaluation of FIM with multiple rays

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
    
    % Precomputation of vectors for FIM
    symb = exp(1j*rand(P,1)*2*pi);
    for ll=1:path_num
        for mm=1:M
            index = (mm-1)*P+1:mm*P;
            vdcfo(index) = g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(dqvec(:,mm))*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
            vdtheta(index) = g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*dtx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
            vdphi(index) = g_ray(ll) * (W(:,mm)'*drx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
            vdtau(index) = g_ray(ll) * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(dfvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
            vdalpha(index) = (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
            vdbeta(index) = 1j * (W(:,mm)'*arx(:,ll)) * conj(F(:,mm)'*atx(:,ll))...
                *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*CFO);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of J_{0,0} part, a 6 by 6 matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_00 = zeros(6,6);

    % Partial of CFO and CFO
    J_00(1,1) = real((vdcfo)' * (vdcfo));
    
    % Partial of CFO and theta
    J_00(1,2) = real((vdcfo)' * (vdtheta));
    J_00(2,1) = J_00(1,2);
    
    % Partial of CFO and phi
    J_00(1,3) = real((vdcfo)' * (vdphi));
    J_00(3,1) = J_00(1,3);
    
    % Partial of CFO and tau
    J_00(1,4) = real((vdcfo)' * (vdtau));
    J_00(4,1) = J_00(1,4);
    
    % Partial of CFO and alpha & beta
    J_00(1,5) = real((vdcfo)' * (vdalpha));
    J_00(5,1) = J_00(1,5);
    
    J_00(1,6) = real((vdcfo)' * (vdbeta));
    J_00(6,1) = J_00(1,6);


    % Partial of theta and theta              
    J_00(2,2) = real((vdtheta)' * (vdtheta));

    % Partial of theta and phi
    J_00(2,3) = real((vdtheta)' * (vdphi));
    J_00(3,2) = J_00(2,3);

    % Partial of theta and tau
    J_00(2,4) = real((vdtheta)' * (vdtau));
    J_00(4,2) = J_00(2,4);
    
    % Partial of theta and gain
    J_00(2,5) = real((vdtheta)' * (vdalpha));
    J_00(5,2) = J_00(2,5);
    
    J_00(2,6) = real((vdtheta)' * (vdbeta));
    J_00(6,2) = J_00(2,6);

    % Partial phi and phi
    J_00(3,3) = real((vdphi)' * (vdphi));
    
    % Partial phi and tau
    J_00(3,4) = real((vdphi)' * (vdtau));
    J_00(4,3) = J_00(3,4);  
    
    % Partial phi and gain
    J_00(3,5) = real((vdphi)' * (vdalpha));
    J_00(5,3) = J_00(3,5);  
    
    J_00(3,6) = real((vdphi)' * (vdbeta));
    J_00(6,3) = J_00(3,6);  
    
    % Partial tau and tau
    J_00(4,4) = real((vdtau)' * (vdtau));
    
    % Partial tau and gain
    J_00(4,5) = real((vdtau)' * (vdalpha));
    J_00(5,4) = J_00(4,5);  
    
    J_00(4,6) = real((vdtau)' * (vdbeta));
    J_00(6,4) = J_00(4,6);

    % Partial gain 
    J_00(5,5) = real((vdalpha)' * (vdalpha));
    J_00(6,6) = real((vdbeta)' * (vdbeta));
    
    
    % For loop for SNR
    for ss = 1:SNR_num
        
        % SNR
        sigman2 = 10^(-SNR_range(ss)/10);
        
        % Evaluate FIM
        J = 2/sigman2 * J_00;
        
        % Evaluation of FIM with single rays
%         temp = inv(J(1:3,1:3));
        temp = inv(J);
        % RMSE evaluation from CRLB perspective
        CRLB_CFO(ss,MCindex) = sqrt(temp(1,1))*(1/Ts/2/pi);
        CRLB_theta(ss,MCindex) = sqrt(temp(2,2))*(1/pi*180);
        CRLB_phi(ss,MCindex) = sqrt(temp(3,3))*(1/pi*180);
    end
end
%% Plot CRLB of angle est./tracking
figure
subplot(211)
semilogy(SNR_range,mean(CRLB_theta,2));hold on
semilogy(SNR_range,mean(CRLB_phi,2));hold on
grid on
legend('Theta','Phi')
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA/AoD [deg]')

subplot(212)
semilogy(SNR_range,mean(CRLB_CFO,2));hold on
grid on
title('CRLB of CFO')
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of CFO [Hz]')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
