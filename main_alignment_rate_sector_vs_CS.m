% Alignment performance evaluation

clear;clc;

% parameters
rand('seed',3)
CFO_num = 1;
CFO_pool = 0;
% freq_offset = 0e3;
BW = 57.6e6; % SS burst BW
Ts = 1/BW; % SS burst sample duration
MCtimes = 1e3; % Monte Carlo simulation trials
P = 127; % num of PSS subcarriers
M = 64; % num of bursts
sample_num = P*M;
PN_sigma = 0^2;
Nt = 128; % Tx antenna
Nr = 32; % Rx antenna
SNR_num = 17;
SNR_range = linspace(-30,10,SNR_num);
runSector = 1;
runCS = 0;
runCSmag = 0;


% parameter for sector beams
% theta_range = linspace(-pi/3,pi/3,N_t);
% phi_range = linspace(-pi/3,pi/3,N_r);

Mt = 16; % transmit sector number
Mr = 4; % reciever sector number
theta_range = linspace(-pi/3,pi/3,Mt);
phi_range = linspace(-pi/3,pi/3,Mr);
M_sector = Mt*Mr;
sample_num_sector = M_sector*P;
for gg=1:Mt
    A_stopband = 10; % attenuation outside mainlobe (dB)
     temp_tx = get_FSM_KW_codebook( theta_range(gg), pi/Mt, Nt, A_stopband);
     txcodebook(:,gg) = temp_tx./norm(temp_tx);
end

for gg=1:Mr
    A_stopband = 10; % attenuation outside mainlobe (dB)
    temp_rx = get_FSM_KW_codebook( phi_range(gg), pi/Mr, Nr, A_stopband);
    rxcodebook(:,gg) = temp_rx./norm(temp_rx);
end

%% dictionary generation

cand_num_r = 17;
cand_num_t = 65;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(1j*(0:Nt-1)'*pi*sin(cand_angle_t));

%% MC simulations

for MCindex=1:MCtimes
    clc
    fprintf('Iteration %d:\n',MCindex);
    
    
    % ------- Generate Channel AOA/AOD and gain ----------
    AOD(MCindex) = (rand*90-45)/180*pi;
    angle_response_tx = transpose(exp(1j*(0:Nt-1)*pi*sin(AOD(MCindex))));
    AOA(MCindex) = (rand*90-45)/180*pi;
    angle_response_rx = transpose(exp(1j*(0:Nr-1)*pi*sin(AOA(MCindex))));
    chan_H = angle_response_rx*angle_response_tx';
    
    
    % ------- CS dictionary generation ----------
    if mod(MCindex,10)==1
        % dictionary generation
        steer_vec_tx = ((randi(2,Nt,M)*2-3)+1j*(randi(2,Nt,M)*2-3))/sqrt(2*Nt);
        steer_vec_rx = ((randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3))/sqrt(2*Nr);
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
    
    % ------- Received signal from PN BF -----------
    y = (steer_vec_tx.'*conj(angle_response_tx)).*(steer_vec_rx'*angle_response_rx);
    
    
    % ------ Received signal from sector beams -----------
    for mm=1:M_sector
        row_now = floor((mm-1)/Mt)+1;
        col_now = mm-(row_now-1)*Mt;
        phi = phi_range(row_now);
        theta = theta_range(col_now);
%         gain = (exp(1j*(0:N_r-1).'*pi*sin(phi))'*chan_H*exp(1j*(0:N_t-1).'*pi*sin(theta)))/sqrt(N_r*N_t);
        gain = rxcodebook(:,row_now)'*chan_H*txcodebook(:,col_now);
        sig_received_sector(:,mm) = ones(P,1)*gain;
    end
    sig_received_sector_vec = reshape(sig_received_sector,P*M_sector,1);
    noise_sector = (randn(sample_num_sector,1)+1j*randn(sample_num_sector,1))/sqrt(2);
    
    
    % ------- Parameter to start EKF (old) ---------
    freq_offset = CFO_pool(1);
    initial_phase = angle(y(1));
    x(:,1) = [initial_phase,freq_offset];
    F = [1, Ts*2*pi; 0, 1];
    w_PN = (rand(sample_num,1)*2-1)*sqrt(3)*sqrt(PN_sigma);

    for ii=1:sample_num-1
        x(:,ii+1) = F * x(:,ii)+[w_PN(ii);0];
        if mod(ii,P)==0
            temp = ii/P;
            x(1,ii+1) = x(1,ii+1)+angle(y(temp+1))-angle(y(temp));
        end
    end

    env = kron(abs(y),ones(P,1));
    sig_received = exp(1j*x(1,:));
    sig_received_debug = exp(1j*x(1,:)).*env';
    noise_normal = (randn(sample_num,1)+1j*randn(sample_num,1))/sqrt(2);

    % --------- Test in different SNR setting ----------
    for SNRindex = 1:SNR_num
        noise_pow = 10^(-SNR_range(SNRindex)/10);
        r_noisy = sig_received.'+noise_normal*sqrt(noise_pow);
        r_noisy_debug = sig_received_debug.'+noise_normal*sqrt(noise_pow);
        
        %------------- Sector Beam Search -------
        if runSector
            r_noisy_sector = sig_received_sector_vec + noise_sector*sqrt(noise_pow);
            r_slot_ave = mean(reshape(r_noisy_sector,P,M_sector),1).';
            [~,bestsector] = max(abs(r_slot_ave));
            bestrow = floor((bestsector-1)/Mt)+1;
            bestcol = bestsector - (bestrow-1)*Mt;
            bestAOA_sector(MCindex,SNRindex) = phi_range(bestrow);
            bestAOD_sector(MCindex,SNRindex) = theta_range(bestcol);
            AOA_error_sector(MCindex,SNRindex) = abs(bestAOA_sector(MCindex,SNRindex) - AOA(MCindex));
            AOD_error_sector(MCindex,SNRindex) = abs(bestAOD_sector(MCindex,SNRindex) - AOD(MCindex));
        end
        
        % ----- OMP use Complex Value without Phase Comp ----------
        if runCS
            for mm=1:M
                sig(MCindex,mm) = mean(r_noisy((mm-1)*P+1:mm*P));
            end

            cand_score_nocomp = (Measure_mat_new'*(sig(MCindex,:).'.*abs(y)))./Measure_mat_new_norm';
            [~,bestindex_nocomp(MCindex)] = max(abs(cand_score_nocomp));

            bestrow = floor((bestindex_nocomp(MCindex)-1)/cand_num_r)+1;
            bestcol = bestindex_nocomp(MCindex)-(bestrow-1)*cand_num_r;
            bestAOA_nocomp(MCindex,SNRindex) = (bestcol-1)*AOAstep-60*pi/180;
            bestAOD_nocomp(MCindex,SNRindex) = (bestrow-1)*AODstep-60*pi/180;

            AOA_error_nocomp(MCindex,SNRindex) = abs(bestAOA_nocomp(MCindex,SNRindex) - AOA(MCindex));
            AOD_error_nocomp(MCindex,SNRindex) = abs(bestAOD_nocomp(MCindex,SNRindex) - AOD(MCindex));
    %         align_counter_nocomp(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_nocomp(runindex))<0.1)&&...
    %             (abs(AOD(runindex) - bestAOD_nocomp(runindex))<0.05);
        end
        
        % ----- OMP using Magnitude Only ----------
        if runCSmag
            for mm=1:M
                sig_post_seq_corr(MCindex,mm) = mean(r_noisy_debug((mm-1)*P+1:mm*P));
            end
            cand_score_mag = (abs(Measure_mat_new)'*abs(sig_post_seq_corr(MCindex,:).'))./Measure_mat_new_norm';
            [~,bestindex_mag(MCindex)] = max(abs(cand_score_mag));

            bestrow_mag = floor((bestindex_mag(MCindex)-1)/cand_num_r)+1;
            bestcol_mag = bestindex_mag(MCindex)-(bestrow_mag-1)*cand_num_r;
            bestAOA_mag(MCindex) = (bestcol_mag-1)*AOAstep-60*pi/180;
            bestAOD_mag(MCindex) = (bestrow_mag-1)*AODstep-60*pi/180;

            AOA_error_mag(MCindex,SNRindex) = abs(bestAOA_mag(MCindex) - AOA(MCindex));
            AOD_error_mag(MCindex,SNRindex) = abs(bestAOD_mag(MCindex) - AOD(MCindex));
    %         align_counter_mag(runindex,CFOindex) = (abs(AOA(runindex) - bestAOA_mag(runindex))<0.1)&&...
    %             (abs(AOD(runindex) - bestAOD_mag(runindex))<0.05);
        end
    end
end
%%
AoA_crit = 8*105/Nr;
AoD_crit = 8*105/Nt;

figure
if runCSmag
    AOAalign_mag_mean = sum((AOA_error_mag/pi*180)<(105/Nr),1)/MCtimes;
    AODalign_mag_mean = sum((AOD_error_mag/pi*180)<(105/Nt),1)/MCtimes;
    
    plot(SNR_range,AOAalign_nocomp_mean,'-','linewidth',2);hold on
    plot(SNR_range,AODalign_nocomp_mean,'-','linewidth',2);hold on
end

if runCS
    AOAalign_nocomp_mean = sum((AOA_error_nocomp/pi*180)<(105/Nr),1)/MCtimes;
    AODalign_nocomp_mean = sum((AOD_error_nocomp/pi*180)<(105/Nt),1)/MCtimes;
    
    plot(SNR_range,AOAalign_mag_mean,'--','linewidth',2);hold on
    plot(SNR_range,AODalign_mag_mean,'--','linewidth',2);hold on
end

if runSector
    AOAalign_sector_mean = sum((AOA_error_sector/pi*180)<AoA_crit,1)/MCtimes;
    AODalign_sector_mean = sum((AOD_error_sector/pi*180)<AoD_crit,1)/MCtimes;
    Align_sector_mean = sum(((AOA_error_sector/pi*180)<AoA_crit)...
                           &((AOD_error_sector/pi*180)<AoD_crit),1)/MCtimes;
    
    semilogy(SNR_range,1-AOAalign_sector_mean,'-o','linewidth',2);hold on
    semilogy(SNR_range,1-AODalign_sector_mean,'-o','linewidth',2);hold on
    semilogy(SNR_range,1-Align_sector_mean,'-o','linewidth',2);hold on

end

grid on
xlabel('SNR [dB]')
ylabel('Mis-Alignment Rate')
ylim([0.01,1])
xlim([-30,10])
% legend('AoA, complex alg','AoD, complex alg','AoA, mag alg','AoD, mag alg','AoA, sector','AoD, sector')
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
xlim([0,105*8/Nt])

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
xlim([0,105*8/Nr])
legend('AoA, complex alg','AoA, mag alg','AoA, sector')
%%
% figure;
% subplot(211)
% plot(abs(cand_score_nocomp))
% subplot(212)
% plot(abs(cand_score_mag))
