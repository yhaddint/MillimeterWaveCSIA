% Alignment performance evaluation

clear;clc;

% parameters
rand('seed',3)
CFO_num = 1;
CFO_pool = 0;
% freq_offset = 0e3;
Ts = 1/(2e6);

% SNR_dB_pool = -15:5:15;
% SNR_pool = 10.^(SNR_dB_pool./10);
% noise_pow_dB_pool = 0;
runtimes = 1e3;
TRN_length = 127;
M = 64;
sample_num = TRN_length*M;
PN_sigma = 0^2;
N_t = 128;
N_r = 16;
SNR_num = 1;
SNR_range = linspace(20,20,SNR_num);


% parameter for sector beams
% theta_range = linspace(-pi/3,pi/3,N_t);
% phi_range = linspace(-pi/3,pi/3,N_r);

G_t = 16;
G_r = 4;
theta_range = linspace(-pi/3,pi/3,G_t);
phi_range = linspace(-pi/3,pi/3,G_r);
M_sector = G_t*G_r;
sample_num_sector = M_sector*TRN_length;
for gg=1:G_t
    A_stopband = 30; % attenuation outside mainlobe (dB)
     temp_tx = get_FSM_KW_codebook( theta_range(gg), pi/G_t, N_t, A_stopband);
     txcodebook(:,gg) = temp_tx./norm(temp_tx);
end

for gg=1:G_r
    A_stopband = 30; % attenuation outside mainlobe (dB)
    temp_rx = get_FSM_KW_codebook( phi_range(gg), pi/G_r, N_r, A_stopband);
    rxcodebook(:,gg) = temp_rx./norm(temp_rx);
end

%% dictionary generation

cand_num_r = 61;
cand_num_t = 121;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:N_r-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:N_t-1)'*pi*sin(cand_angle_t));

%% MC simulations

for runindex=1:runtimes
    if mod(runindex,10)==1
        clc;fprintf('Ite %d out of %d\n',runindex,runtimes);
        % dictionary generation
        steer_vec_tx = ((randi(2,N_t,M)*2-3)+1j*(randi(2,N_t,M)*2-3))/sqrt(2*N_t);
        steer_vec_rx = ((randi(2,N_r,M)*2-3)+1j*(randi(2,N_r,M)*2-3))/sqrt(2*N_r);
        Measure_mat = kron(transpose(steer_vec_tx)*conj(cand_ARV_t),steer_vec_rx'*cand_ARV_r);
        select_row = zeros(1,M);
        for ii=1:M
            select_row(ii) = (ii-1)*M+ii;
        end
        Measure_mat_new = Measure_mat(select_row,:);
        for cc=1:cand_num_r*cand_num_t
            Measure_mat_new_norm(cc) = norm(Measure_mat_new(:,cc),2)^2;
        end

        Measure_mat_new_diff(1,:) = Measure_mat_new(1,:);
        for mm=2:M
            Measure_mat_new_diff(mm,:) = Measure_mat_new(mm,:).*exp(-1j*phase(Measure_mat_new(mm-1,:)));
        end
        
    end
    % generate AOA/AOD and g_m
    AOD(runindex) = 0;%(rand*120-60)/180*pi;
    angle_response_tx = transpose(exp(1j*(0:N_t-1)*pi*sin(AOD(runindex))));
    AOA(runindex) = 0;%(rand*120-60)/180*pi;
    angle_response_rx = transpose(exp(1j*(0:N_r-1)*pi*sin(AOA(runindex))));
    chan_H = angle_response_rx*angle_response_tx';
    y = (steer_vec_tx.'*conj(angle_response_tx)).*(steer_vec_rx'*angle_response_rx);
    
    % received signal from sector beams
    
    for mm=1:M_sector
        row_now = floor((mm-1)/G_t)+1;
        col_now = mm-(row_now-1)*G_t;
        phi = phi_range(row_now);
        theta = theta_range(col_now);
%         gain = (exp(1j*(0:N_r-1).'*pi*sin(phi))'*chan_H*exp(1j*(0:N_t-1).'*pi*sin(theta)))/sqrt(N_r*N_t);
        gain = rxcodebook(:,row_now)'*chan_H*txcodebook(:,col_now);
        sig_received_sector(:,mm) = ones(TRN_length,1)*gain;
    end
    sig_received_sector_vec = reshape(sig_received_sector,TRN_length*M_sector,1);
    noise_sector = (randn(sample_num_sector,1)+1j*randn(sample_num_sector,1))/sqrt(2);
    
    
    % parameter to start EKF
    
        freq_offset = CFO_pool(1);
        initial_phase = angle(y(1));
        x(:,1) = [initial_phase,freq_offset];
        F = [1, Ts*2*pi; 0, 1];
        w_PN = (rand(sample_num,1)*2-1)*sqrt(3)*sqrt(PN_sigma);

        for ii=1:sample_num-1
            x(:,ii+1) = F * x(:,ii)+[w_PN(ii);0];
            if mod(ii,TRN_length)==0
                temp = ii/TRN_length;
                x(1,ii+1) = x(1,ii+1)+angle(y(temp+1))-angle(y(temp));
            end
        end
        
        env = kron(abs(y),ones(TRN_length,1));
        sig_received = exp(1j*x(1,:));
        sig_received_debug = exp(1j*x(1,:)).*env';
        noise_normal = (randn(sample_num,1)+1j*randn(sample_num,1))/sqrt(2);

    
    for SNRindex = 1:SNR_num
        noise_pow = 10^(-SNR_range(SNRindex)/10);
        r_noisy = sig_received.'+noise_normal*sqrt(noise_pow);
        r_noisy_debug = sig_received_debug.'+noise_normal*sqrt(noise_pow);
        
        %------------- Sector beam search -------
        r_noisy_sector = sig_received_sector_vec + noise_sector*sqrt(noise_pow);
        r_slot_ave = mean(reshape(r_noisy_sector,TRN_length,M_sector),1).';
        [~,bestsector] = max(abs(r_slot_ave));
        bestrow = floor((bestsector-1)/G_t)+1;
        bestcol = bestsector - (bestrow-1)*G_t;
        bestAOA_sector(runindex,SNRindex) = phi_range(bestrow);
        bestAOD_sector(runindex,SNRindex) = theta_range(bestcol);
        AOA_error_sector(runindex,SNRindex) = abs(bestAOA_sector(runindex,SNRindex) - AOA(runindex));
        AOD_error_sector(runindex,SNRindex) = abs(bestAOD_sector(runindex,SNRindex) - AOD(runindex));
        
        % ----- OMP use Complex Value without Phase Comp ----------
        for mm=1:M
            sig(runindex,mm) = mean(r_noisy((mm-1)*TRN_length+1:mm*TRN_length));
        end
        
        cand_score_nocomp = (Measure_mat_new'*(sig(runindex,:).'.*abs(y)))./Measure_mat_new_norm';
        [~,bestindex_nocomp(runindex)] = max(abs(cand_score_nocomp));

        bestrow = floor((bestindex_nocomp(runindex)-1)/cand_num_r)+1;
        bestcol = bestindex_nocomp(runindex)-(bestrow-1)*cand_num_r;
        bestAOA_nocomp(runindex,SNRindex) = (bestcol-1)*AOAstep-60*pi/180;
        bestAOD_nocomp(runindex,SNRindex) = (bestrow-1)*AODstep-60*pi/180;

        AOA_error_nocomp(runindex,SNRindex) = abs(bestAOA_nocomp(runindex,SNRindex) - AOA(runindex));
        AOD_error_nocomp(runindex,SNRindex) = abs(bestAOD_nocomp(runindex,SNRindex) - AOD(runindex));
%         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);
        
        % ----- OMP using Magnitude Only ----------
        for mm=1:M
            sig_post_seq_corr(runindex,mm) = mean(r_noisy_debug((mm-1)*TRN_length+1:mm*TRN_length));
        end
        C_mat = Measure_mat_new;
        d_vec = sig_post_seq_corr.';
        cvx_begin
            variable xxx(n)
            minimize( norm(xxx, 2 ) )
            subject to
                C_mat * xxx == d_vec
        cvx_end
        xxx
        cand_score_mag = (abs(Measure_mat_new)'*abs(sig_post_seq_corr(runindex,:).'))./Measure_mat_new_norm';
        [~,bestindex_mag(runindex)] = max(abs(cand_score_mag));

        bestrow_mag = floor((bestindex_mag(runindex)-1)/cand_num_r)+1;
        bestcol_mag = bestindex_mag(runindex)-(bestrow_mag-1)*cand_num_r;
        bestAOA_mag(runindex) = (bestcol_mag-1)*AOAstep-60*pi/180;
        bestAOD_mag(runindex) = (bestrow_mag-1)*AODstep-60*pi/180;

        AOA_error_mag(runindex,SNRindex) = abs(bestAOA_mag(runindex) - AOA(runindex));
        AOD_error_mag(runindex,SNRindex) = abs(bestAOD_mag(runindex) - AOD(runindex));
%         align_counter_mag(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_mag(runindex))<0.1)&&...
%             (abs(AOD(runindex) - bestAOD_mag(runindex))<0.05);
        
    end
end
%%
AOAalign_mag_mean = sum((AOA_error_mag/pi*180)<(105/N_r),1)/runtimes;
AODalign_mag_mean = sum((AOD_error_mag/pi*180)<(105/N_t),1)/runtimes;
AOAalign_nocomp_mean = sum((AOA_error_nocomp/pi*180)<(105/N_r),1)/runtimes;
AODalign_nocomp_mean = sum((AOD_error_nocomp/pi*180)<(105/N_t),1)/runtimes;
AOAalign_sector_mean = sum((AOA_error_sector/pi*180)<(105/N_r),1)/runtimes;
AODalign_sector_mean = sum((AOD_error_sector/pi*180)<(105/N_t),1)/runtimes;
figure
plot(SNR_range,AOAalign_nocomp_mean,'-','linewidth',2);hold on
plot(SNR_range,AODalign_nocomp_mean,'-','linewidth',2);hold on
plot(SNR_range,AOAalign_mag_mean,'--','linewidth',2);hold on
plot(SNR_range,AODalign_mag_mean,'--','linewidth',2);hold on
plot(SNR_range,AOAalign_sector_mean,'-o','linewidth',2);hold on
plot(SNR_range,AODalign_sector_mean,'-o','linewidth',2);hold on
grid on
xlabel('SNR (dB)')
ylabel('Alignment Prob')
legend('AoA, complex alg','AoD, complex alg','AoA, mag alg','AoD, mag alg','AoA, sector','AoD, sector')
%%
figure(99)
subplot(211)
[b,a] = ecdf(AOD_error_nocomp(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOD_error_mag(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOD_error_sector(:,5)/pi*180);
plot(a,b);hold on
grid on
xlabel('Estimation Error [deg]')
ylabel('CDF')
legend('AoD, complex alg','AoD, mag alg','AoD, sector')
xlim([0,105*8/N_t])

subplot(212)
[b,a] = ecdf(AOA_error_nocomp(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOA_error_mag(:,5)/pi*180);
plot(a,b);hold on
[b,a] = ecdf(AOA_error_sector(:,5)/pi*180);
plot(a,b);hold on
grid on
xlabel('Estimation Error [deg]')
ylabel('CDF')
xlim([0,105*8/N_r])
legend('AoA, complex alg','AoA, mag alg','AoA, sector')
%%
% figure;
% subplot(211)
% plot(abs(cand_score_nocomp))
% subplot(212)
% plot(abs(cand_score_mag))
