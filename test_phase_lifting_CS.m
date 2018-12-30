% Alignment performance evaluation

clear;clc;

% parameters
rand('seed',1)
CFO_num = 1;
CFO_pool = 0;
% freq_offset = 0e3;
Ts = 1/(2e6);

% SNR_dB_pool = -15:5:15;
% SNR_pool = 10.^(SNR_dB_pool./10);
% noise_pow_dB_pool = 0;
runtimes = 5e1;
TRN_length = 127;
M = 64;
MCS = 15;
sample_num = TRN_length*M;
Nt = 32;
Nr = 8;
SNR_num = 2;
SNR_range = linspace(10,20,SNR_num);

%% dictionary generation

cand_num_r = 9;
cand_num_t = 33;
GtGr = cand_num_r*cand_num_t;
NtNr = Nt*Nr;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Nt-1)'*pi*sin(cand_angle_t));

%% MC simulations

for runindex=1:runtimes
    clc;fprintf('Ite %d out of %d\n',runindex,runtimes);
    if mod(runindex,10)==1
        % dictionary generation
        steer_vec_tx = ((randi(2,Nt,M)*2-3)+1j*(randi(2,Nt,M)*2-3))/sqrt(2*Nt);
        steer_vec_rx = ((randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3))/sqrt(2*Nr);
        Atot_temp = kron(transpose(steer_vec_tx),steer_vec_rx');
        select_row = zeros(1,M);
        for ii=1:M
            select_row(ii) = (ii-1)*M+ii;
        end
        Atot = Atot_temp(select_row,:);
        ACS = (randn(MCS,NtNr) + 1j*randn(MCS,NtNr))*sqrt(1/2/NtNr);
%         APR = Atot* (ACS'*inv(ACS*ACS'));
%         for bb=1:64
%             A = ACS';
%             B = Atot';
%             b = B(:,bb);
%             x = pinv(A)*b;
%             X(:,bb) = x;
%         end
%         APR = X';
        
        APR = (randn(M,MCS) + 1j*randn(M,MCS))*sqrt(1/2/MCS);
        Atot_decomp = APR*ACS;
        
        new_dict = kron(conj(cand_ARV_t),cand_ARV_r); 
        
%         Measure_mat = kron(transpose(steer_vec_tx)*conj(cand_ARV_t),steer_vec_rx'*cand_ARV_r);
%         select_row = zeros(1,M);
%         for ii=1:M
%             select_row(ii) = (ii-1)*M+ii;
%         end
        Measure_mat_new = ACS*new_dict;
        for cc=1:cand_num_r*cand_num_t
            Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
        end
        
        
        
        
    end
    
    % Find True AoA/AoD in grid (for debug)
    [~,row_true] = min(abs(cand_angle_r - 0));
    [~,col_true] = min(abs(cand_angle_t - 0));
    index_true = (col_true-1)*cand_num_r + row_true;
    
%     % generate AOA/AOD and g_m
%     AOD(runindex) = 0;%(rand*120-60)/180*pi;
%     angle_response_tx = transpose(exp(1j*(0:N_t-1)*pi*sin(AOD(runindex))));
%     AOA(runindex) = 0;%(rand*120-60)/180*pi;
%     angle_response_rx = transpose(exp(1j*(0:N_r-1)*pi*sin(AOA(runindex))));
%     chan_H = angle_response_rx*angle_response_tx';
%     y = (steer_vec_tx.'*conj(angle_response_tx)).*(steer_vec_rx'*angle_response_rx);
%     

    
%     figure;plot(abs(y_alt - y));
%     hold on;
%     plot(abs(y))
%     legend('error','y')

    
%     error_norm = norm(Atot - APR * ACS,'fro')
%     sig_norm = norm(Atot,'fro')
%     apple=1;
% %        
%         freq_offset = CFO_pool(1);
%         initial_phase = angle(y(1));
%         x(:,1) = [initial_phase,freq_offset];
%         F = [1, Ts*2*pi; 0, 1];
%         w_PN = (rand(sample_num,1)*2-1)*sqrt(3)*sqrt(PN_sigma);
% 
%         for ii=1:sample_num-1
%             x(:,ii+1) = F * x(:,ii)+[w_PN(ii);0];
%             if mod(ii,TRN_length)==0
%                 temp = ii/TRN_length;
%                 x(1,ii+1) = x(1,ii+1)+angle(y(temp+1))-angle(y(temp));
%             end
%         end
%         
%         env = kron(abs(y),ones(TRN_length,1));
%         sig_received = exp(1j*x(1,:));
%         sig_received_debug = exp(1j*x(1,:)).*env';
        noise_normal = (randn(M,1)+1j*randn(M,1))/sqrt(2);

    
    for SNRindex = 1:SNR_num
        
        noise_pow = 10^(-SNR_range(SNRindex)/10);
        y = APR * (ACS*new_dict(:,3691)) + noise_normal*sqrt(noise_pow);
%         r_noisy = sig_received.'+noise_normal*sqrt(noise_pow);
%         r_noisy_debug = sig_received_debug.'+noise_normal*sqrt(noise_pow);
%         
       
        % ----- OMP use Complex Value without Phase Comp ----------
%         for mm=1:M
%             sig(runindex,mm) = mean(r_noisy((mm-1)*TRN_length+1:mm*TRN_length));
%         end
%         
%         cand_score_nocomp = (Measure_mat_new'*(sig(runindex,:).'.*abs(y)))./Measure_mat_new_norm';
%         [~,bestindex_nocomp(runindex)] = max(abs(cand_score_nocomp));
% 
%         bestrow = floor((bestindex_nocomp(runindex)-1)/cand_num_r)+1;
%         bestcol = bestindex_nocomp(runindex)-(bestrow-1)*cand_num_r;
%         bestAOA_nocomp(runindex,SNRindex) = (bestcol-1)*AOAstep-60*pi/180;
%         bestAOD_nocomp(runindex,SNRindex) = (bestrow-1)*AODstep-60*pi/180;
% 
%         AOA_error_nocomp(runindex,SNRindex) = abs(bestAOA_nocomp(runindex,SNRindex) - AOA(runindex));
%         AOD_error_nocomp(runindex,SNRindex) = abs(bestAOD_nocomp(runindex,SNRindex) - AOD(runindex));
% %         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
% %             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);
%         
        % ----- OMP using Magnitude Only ----------
%         for mm=1:M
%             sig_post_seq_corr(runindex,mm) = mean(r_noisy_debug((mm-1)*TRN_length+1:mm*TRN_length));
%         end
        d_vec = (abs(y)).^2;
%         seldiag = zeros(1024,1);
%         for ii=1:32
%             for jj=1:32
%                 index = (ii-1)*32+jj;
%                 Asdp(:,index) = APR(:,ii).*conj(APR(:,jj));
%                 if ii==jj
%                     seldiag(index) = 1;
%                 end
%             end
%         end
        cvx_begin sdp quiet
            variable Z(MCS,MCS) hermitian semidefinite 
            minimize( trace(Z) )
            subject to
                norm(diag(APR * Z * APR') - d_vec,2) <= 1e-1
        cvx_end
        
        z_new = (ACS*new_dict(:,3691));
        Z_new = z_new*z_new';
        norm(diag(APR * Z_new * APR') - d_vec);
        
        [Umat,Sigma,Vmat] = svd(Z);
        alpha = pinv(Vmat(:,1))*z_new;
        norm(alpha*Vmat(:,1)-z_new)/norm(z_new);
        
        cand_score_nocomp = (Measure_mat_new'*Vmat(:,1))./Measure_mat_new_norm';
        [~,bestindex_nocomp(runindex)] = max(abs(cand_score_nocomp));

        bestrow = floor((bestindex_nocomp(runindex)-1)/cand_num_r)+1;
        bestcol = bestindex_nocomp(runindex)-(bestrow-1)*cand_num_r;
        bestAOA_nocomp(runindex,SNRindex) = (bestcol-1)*AOAstep-60*pi/180;
        bestAOD_nocomp(runindex,SNRindex) = (bestrow-1)*AODstep-60*pi/180;

        AOA_error_nocomp(runindex,SNRindex) = abs(bestAOA_nocomp(runindex,SNRindex) - 0);
        AOD_error_nocomp(runindex,SNRindex) = abs(bestAOD_nocomp(runindex,SNRindex) - 0);
%         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);
        
    end
end
%%
% AOAalign_mag_mean = sum((AOA_error_mag/pi*180)<(105/N_r),1)/runtimes;
% AODalign_mag_mean = sum((AOD_error_mag/pi*180)<(105/N_t),1)/runtimes;
AOAalign_nocomp_mean = sum((AOA_error_nocomp/pi*180)<(105/Nr),1)/runtimes;
AODalign_nocomp_mean = sum((AOD_error_nocomp/pi*180)<(105/Nt),1)/runtimes;
figure
plot(SNR_range,AOAalign_nocomp_mean,'-','linewidth',2);hold on
plot(SNR_range,AODalign_nocomp_mean,'-','linewidth',2);hold on
% plot(SNR_range,AOAalign_mag_mean,'--','linewidth',2);hold on
% plot(SNR_range,AODalign_mag_mean,'--','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legend('AoA, complex alg','AoD, complex alg')
%%
figure(99)
subplot(211)
[b,a] = ecdf(AOD_error_nocomp(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOD_error_mag(:,5)/pi*180);
plot(a,b);hold on
grid on
xlabel('Estimation Error [deg]')
ylabel('CDF')
legend('AoD, complex alg','AoD, mag alg')
xlim([0,105*8/Nt])

subplot(212)
[b,a] = ecdf(AOA_error_nocomp(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOA_error_mag(:,5)/pi*180);
plot(a,b);hold on
grid on
xlabel('Estimation Error [deg]')
ylabel('CDF')
xlim([0,105*8/Nr])
legend('AoA, complex alg','AoA, mag alg')
