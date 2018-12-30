% Alignment performance evaluation

clear;clc;

% ---------- system parameters -----------
rand('seed',3)                                              % random seed
MCtimes = 1e3;                                              % Monte CArlo Runs

CFO_num = 1;
CFO_pool = 0;
% freq_offset = 0e3;                                        % CFO value [Hz]
BW = 57.6e6;                                                % Sample rate [Hz]
Ts = 1/BW;                                                  % Sample Duration [second]
TRN_length = 127;                                           % Pilots number in each sounding
M = 16;                                                    % Num of CS measurement
sample_num = TRN_length*M;                                  % Tot Sample num
PN_sigma = 0^2;                                             % Phase noise power
Ntx = 32;                                                   % transmitter array ULA
Nrx = 8;                                                    % transmitter array ULA

SNR_num = 21;                                                % SNR sweep number
SNR_range = linspace(-40,10,SNR_num);                        % SNR range of interest


%% Pre-compute array responses in all angles

Gt = 33;
Gr = 9;
cand_num_r = Gr;
cand_num_t = Gt;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nrx-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Ntx-1)'*pi*sin(cand_angle_t));

%% MC simulations
for MCidx = 1:MCtimes
    clc;fprintf('Ite %d out of %d\n',MCidx,MCtimes);
    
    % To test random beamformer, the sounder BF weights shuffle from time
    % to time
    if mod(MCidx,20)==1
        steer_vec_tx = ((randi(2,Ntx,M)*2-3)+1j*(randi(2,Ntx,M)*2-3))/sqrt(2*Ntx);
        steer_vec_rx = ((randi(2,Nrx,M)*2-3)+1j*(randi(2,Nrx,M)*2-3))/sqrt(2*Nrx);
        Measure_mat = kron(transpose(steer_vec_tx)*conj(cand_ARV_t),steer_vec_rx'*cand_ARV_r);
        select_row = zeros(1,M);
        for ii=1:M
            select_row(ii) = (ii-1)*M+ii;
        end
        Measure_mat_new = Measure_mat(select_row,:);
        for cc=1:cand_num_r*cand_num_t
            Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
        end
    end
    
    
    % Generate AOA/AOD and g_m and CS measurement
    AOD(MCidx) = 0;%(rand*120-60)/180*pi;
    angle_response_tx = (exp(1j*(0:Ntx-1)*pi*sin(AOD(MCidx)))).';
    AOA(MCidx) = 0;%(rand*120-60)/180*pi;
    angle_response_rx = (exp(1j*(0:Nrx-1)*pi*sin(AOA(MCidx)))).';
    chan_H = angle_response_rx*angle_response_tx';
    y = (steer_vec_tx.'*conj(angle_response_tx)).*...
        (steer_vec_rx'*angle_response_rx);
    
    % Received signal, currently no phase error
    sig_received = (kron(y,ones(TRN_length,1))).';
    noise_normal = (randn(sample_num,1)+1j*randn(sample_num,1))/sqrt(2);

    % For loop of different SNR
    for SNRindex = 1:SNR_num
        
        % Adding noise with different power
        noise_pow = 10^(-SNR_range(SNRindex)/10);
        r_noisy = sig_received.'+noise_normal*sqrt(noise_pow);
        

        % ----- OMP use Complex Value ( no Phase Error/Comp) ----------
        for mm=1:M
            sig(MCidx,mm) = mean(r_noisy((mm-1)*TRN_length+1:mm*TRN_length));
        end
        
        cand_score_nocomp = (Measure_mat_new'*(sig(MCidx,:).'.*abs(y)))./Measure_mat_new_norm';
        [~,bestindex_nocomp(MCidx)] = max(abs(cand_score_nocomp));

        bestrow = floor((bestindex_nocomp(MCidx)-1)/cand_num_r)+1;
        bestcol = bestindex_nocomp(MCidx)-(bestrow-1)*cand_num_r;
        bestAOA_nocomp(MCidx,SNRindex) = (bestcol-1)*AOAstep-60*pi/180;
        bestAOD_nocomp(MCidx,SNRindex) = (bestrow-1)*AODstep-60*pi/180;

        AOA_error_nocomp(MCidx,SNRindex) = abs(bestAOA_nocomp(MCidx,SNRindex) - AOA(MCidx));
        AOD_error_nocomp(MCidx,SNRindex) = abs(bestAOD_nocomp(MCidx,SNRindex) - AOD(MCidx));
%         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);
        

        % ----- OMP using RSSI Only ----------
        for mm=1:M
            sig_post_seq_corr(MCidx,mm) = mean(r_noisy((mm-1)*TRN_length+1:mm*TRN_length));
        end

        cand_score_mag = (abs(Measure_mat_new)'*abs(sig_post_seq_corr(MCidx,:).'))./Measure_mat_new_norm';
        [~,bestindex_mag(MCidx)] = max(abs(cand_score_mag));

        bestrow_mag = floor((bestindex_mag(MCidx)-1)/cand_num_r)+1;
        bestcol_mag = bestindex_mag(MCidx)-(bestrow_mag-1)*cand_num_r;
        bestAOA_mag(MCidx) = (bestcol_mag-1)*AOAstep-60*pi/180;
        bestAOD_mag(MCidx) = (bestrow_mag-1)*AODstep-60*pi/180;

        AOA_error_mag(MCidx,SNRindex) = abs(bestAOA_mag(MCidx) - AOA(MCidx));
        AOD_error_mag(MCidx,SNRindex) = abs(bestAOD_mag(MCidx) - AOD(MCidx));
%         align_counter_mag(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_mag(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_mag(runindex))<0.05);
        
    end % end of SNR for loop
end % end of MC loop
%%
AOAalign_mag_mean = sum((AOA_error_mag/pi*180)<(105/Nrx),1)/MCtimes;
AODalign_mag_mean = sum((AOD_error_mag/pi*180)<(105/Ntx),1)/MCtimes;
align_mag_mean = sum((AOA_error_mag/pi*180)<(105/Nrx)& ...
                     (AOD_error_mag/pi*180)<(105/Ntx), 1)/MCtimes;
                 
AOAalign_nocomp_mean = sum((AOA_error_nocomp/pi*180)<(105/Nrx),1)/MCtimes;
AODalign_nocomp_mean = sum((AOD_error_nocomp/pi*180)<(105/Ntx),1)/MCtimes;
align_nocomp_mean = sum((AOA_error_nocomp/pi*180)<(105/Nrx)& ...
                        (AOD_error_nocomp/pi*180)<(105/Ntx), 1)/MCtimes;
                 
figure
% plot(SNR_range,AOAalign_nocomp_mean,'-','linewidth',2);hold on
% plot(SNR_range,AODalign_nocomp_mean,'-','linewidth',2);hold on
% plot(SNR_range,AOAalign_mag_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_mag_mean,'--','linewidth',2);hold on
plot(SNR_range,align_nocomp_mean,'-o','linewidth',2,'markersize',10);hold on
plot(SNR_range,align_mag_mean,'--x','linewidth',2,'markersize',10);hold on

grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
% legend('AoA, complex alg','AoD, complex alg','AoA, mag alg','AoD, mag alg')
legend(['Coherent, M=' num2str(M)],['Noncoherent, M=' num2str(M)])
%%
% figure(99)
% subplot(211)
% [b,a] = ecdf(AOD_error_nocomp(:,5)/pi*180);
% plot(a,b);hold on
% [b,a] = ecdf(AOD_error_mag(:,5)/pi*180);
% plot(a,b);hold on
% grid on
% xlabel('Estimation Error [deg]')
% ylabel('CDF')
% legend('AoD, complex alg','AoD, mag alg')
% xlim([0,105*8/Ntx])
% 
% subplot(212)
% [b,a] = ecdf(AOA_error_nocomp(:,5)/pi*180);
% plot(a,b);hold on
% [b,a] = ecdf(AOA_error_mag(:,5)/pi*180);
% plot(a,b);hold on
% grid on
% xlabel('Estimation Error [deg]')
% ylabel('CDF')
% xlim([0,105*8/Nrx])
% legend('AoA, complex alg','AoA, mag alg')
