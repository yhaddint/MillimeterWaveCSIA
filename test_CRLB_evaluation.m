%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 1; % Num of rays in a cluster
Nr = 8; % Number of antenna in Rx
Nt = 32;
M = 64; % Length of training
MCtimes = 20; % Num of Monte Carlo Sim.
AOAspread2 = 0;
AOAspread = 0;
AODspread2 = 0;
AODspread = 0;
SNR_num = 50;
SNR_range = linspace(-15,30,SNR_num);
Ts = 1/(50e6);
CFO_ppm = 5; % CFO in ppm
CFO = 28e9/1e6*CFO_ppm; % With unit Hz
eF = CFO*Ts*2*pi; % 

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
    phi = zeros(ray_num,1);
    phi0 = 0/180*pi;
    phi = phi0 + randn(ray_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta = zeros(ray_num,1);
    theta0 = 0/180*pi;
    theta = theta0 + randn(ray_num,1) * AODspread;

%     % Gain
%     g_cmplx = exp(1j*rand(ray_num,1)*2*pi)/sqrt(ray_num);
%     g = g_cmplx;
    % Rotate of ray
    tau = rand * (200e-9);
    g_ray = (randn+1j*randn)/sqrt(2);



    % Pre-compute some vectors/matrices in FIM
    for rayindex = 1:ray_num

        % Spatial response and its derivative over phi
        arx(:,rayindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(rayindex)))/sqrt(Nr);
        Darx(:,rayindex) = 1j * pi * (0:Nr-1)' * cos(phi(rayindex));
        drx(:,rayindex) = Darx(:,rayindex).*arx(:,rayindex);

        % Spatial response and its derivative over theta
        atx(:,rayindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(rayindex)))/sqrt(Nt);
        Datx(:,rayindex) = 1j * pi * (0:Nt-1)' * cos(theta(rayindex));
        dtx(:,rayindex) = Datx(:,rayindex).*atx(:,rayindex);
        
        % Delay response and its derivative over tau
        fvec(:,rayindex) = exp(1j * 2*pi * (0:P-1)' * tau(rayindex) * Ts);
        Dfvec(:,rayindex) = 1j * 2*pi * (0:P-1)' * Ts;
        dfvec(:,rayindex) = Dfvec(:,rayindex).*fvec(:,rayindex);
    end
    
    % About CFO and its derivative
    qvec = exp(1j * (0:P-1)' * eF);
    Dqvec = 1j * (p:P-1)';
    dqvec = qvec.*Dqvec;

    % Evaluation of FIM with multiple rays

    % Pre-compute some equations, vecx, vecalpha, vecbeta
    % vecx is derivative of phi_0 in f(r,eta)
    vaa = zeros(M,1);
    vad = zeros(M,1);
    vda = zeros(M,1);

    for rr=1:ray_num
        vaa = vaa + g_ray(rr) * (W'*arx(:,rr)).* conj(F'*atx(:,rr));
        vad = vad + g_ray(rr) * (W'*arx(:,rr)).* conj(F'*dtx(:,rr));
        vda = vda + g_ray(rr) * (W'*drx(:,rr)).* conj(F'*atx(:,rr));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of J_{0,0} part, a 5 by 5 matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_00 = zeros(6,6);

    % Partial of CFO and CFO
    J_00(1,1) = real((vda)' * (vda));
    
    % Partial of CFO and theta
    J_00(1,2) = real((vda)' * (vda));
    J_00(2,1) = J_00(1,2);
    
    % Partial of CFO and phi
    J_00(1,3) = real((vda)' * (vda));
    J_00(3,1) = J_00(1,3);
    
    % Partial of CFO and tau
    J_00(1,4) = real((vda)' * (vda));
    J_00(4,1) = J_00(1,4);
    
    % Partial of CFO and alpha & beta
    J_00(1,5) = real((vda)' * (vda));
    J_00(5,1) = J_00(1,5);
    
    J_00(1,6) = real((vda)' * (vda));
    J_00(6,1) = J_00(1,6);


    % Partial of theta and theta              
    J_00(2,2) = real((vda)' * (vda));

    % Partial of theta and phi
    J_00(2,3) = real((vda)' * (vad));
    J_00(3,2) = J_00(2,3);

    % Partial of theta and tau
    J_00(2,4) = real((vda)' * (1*vaa));
    J_00(4,2) = J_00(2,4);
    
    % Partial of theta and gain
    J_00(2,5) = real((vda)' * (1*vaa));
    J_00(5,2) = J_00(2,5);
    
    J_00(2,6) = real((vda)' * (1*vaa));
    J_00(6,2) = J_00(2,6);

    % Partial phi and phi
    J_00(3,3) = real((vad)' * (vad));
    
    % Partial phi and tau
    J_00(3,4) = real((vda)' * (1*vaa));
    J_00(4,3) = J_00(3,4);  
    
    % Partial phi and gain
    J_00(3,5) = real((vda)' * (1*vaa));
    J_00(5,3) = J_00(3,5);  
    
    J_00(3,6) = real((vda)' * (1*vaa));
    J_00(6,3) = J_00(3,6);  
    
    % Partial tau and tau
    J_00(4,4) = real((vad)' * (vad));
    
    % Partial tau and gain
    J_00(4,5) = real((vda)' * (1*vaa));
    J_00(5,4) = J_00(4,5);  
    
    J_00(4,6) = real((vda)' * (1*vaa));
    J_00(6,4) = J_00(4,6);

    % Partial gain 
    J_00(5,5) = real((vad)' * (1*vaa));
    J_00(6,6) = real((vad)' * (1*vaa));
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J matrix for apriori distribution of dphi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_p = zeros(3+ray_num*3);
    % A-priori of angular spread
    for r1 = 4:3:(ray_num*3+3)
        J_p(r1,r1) = 2/AOAspread;
    end
    for r2 = 5:3:(ray_num*3+3)
        if AODspread==0
            J_p(r2,r2) = 0;
        else
            J_p(r2,r2) = 2/AODspread;
        end
    end
    
    % For loop for SNR
    for ss = 1:SNR_num
        
        % SNR
        sigman2 = 10^(-SNR_range(ss)/10);
        
        % Evaluate FIM
        J = 2/sigman2*[J_00, J_0D; J_D0, J_DD] + J_p;
        
        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        dphiindex = 4:3:(ray_num*3+3);
        dpsiindex = 6:3:(ray_num*3+3);
        usefulindex = reshape([dphiindex;dpsiindex],1,ray_num*2);
        temp = inv(J);
%         temp = inv(J([1,3,usefulindex],[1,3,usefulindex]));
        CRLB_multiple(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);
        
        % Evaluation of FIM with single rays
%         temp = inv(J(1:3,1:3));
        temp = inv([J(1,1),J(1,3);J(3,1),J(3,3)]);
        % RMSE evaluation from CRLB perspective
        CRLB_single(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);

    end
end
%% Plot CRLB of angle est./tracking
figure
semilogy(SNR_range,mean(CRLB_single,2));hold on
semilogy(SNR_range,mean(CRLB_multiple,2));hold on
grid on
legend('Single Ray','Multple Rays')
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA [deg]')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
