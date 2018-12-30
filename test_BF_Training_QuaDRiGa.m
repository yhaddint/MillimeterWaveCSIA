% My test file for multi-antenna channel generation

%% Channel model setup and coefficient generation
% First, we parameterize the channel model. We start with the basic simulation parameters. For the
% desired output, we need two additional options: we want to evaluate absolute delays and we need to
% get all 20 sub-paths. Normally, the sub-paths are added already in the channel builder.   

clear all
close all
addpath('C:\Users\Han\Documents\QuaDriGa_2017.08.01_v2.0.0-664\quadriga_src')
rng(4);                                                 % random seed

set(0,'defaultTextFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontName','Times')                    % Default Font Type
set(0,'defaultTextFontName','Times')                    % Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type
set(0,'DefaultFigurePaperSize',[14.5 7.3])              % Default Paper Size

s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 28e9;                            % 2.53 GHz carrier frequency
s.sample_density = 4;                                   % 4 samples per half-wavelength
s.use_absolute_delays = 0;                              % Include delay of the LOS path
s.show_progress_bars = 0;                               % Disable progress bars

%%
% Second, we define a user track. Here we choose a linear track with a length of 30 m. The track
% start 20 m east of the transmitter and runs in east direction, thus linearly increasing the
% distance from the receiver.   
Ntx = 32;
Nrx = 8;
l = qd_layout( s );                                     % New QuaDRiGa layout
l.tx_array = qd_arrayant('3gpp-mmw',1,Ntx,28e9,1,0,0.5,1,1);
l.rx_array = qd_arrayant('3gpp-mmw',1,Nrx,28e9,1,0,0.5,1,1);
% l.rx_array.rotate_pattern(180, 'z');

% l.rx_array.visualize(1);pause(1)

l.track = qd_track('linear',1,-pi);

% l.no_tx = 16;                                           % Transmitter antenna in BS
% l.no_rx = 4;                                            % Receiver antenna in UE

l.tx_position(:,1) = [0,0,25].';                        % 25 m BE height
l.rx_position(:,1) = [20,20/sqrt(3),25].';                        % 25 m BE height
l.set_scenario('mmMAGIC_UMi_LOS');                    % Set propagation scenario
% l.set_scenario('3GPP_38.901_UMi_LOS');                    % Set propagation scenario
% l.visualize;                                            % Plot the layout

%%
% Now, we generate the LSPs. We set the shadow fading and K-factor to 1 and disable the path loss
% model. 

cb = l.init_builder;  
% cb.scenpar.PerClusterDS = 0;                            % Create new builder object
cb.scenpar.SF_sigma = 0;                                % 0 dB shadow fading
cb.scenpar.KF_mu = 0;                                   % 0 dB K-Factor
cb.scenpar.KF_sigma = 0;                                % No KF variation
cb.scenpar.SubpathMethod = 'mmMAGIC';
% cb.scenpar.SubpathMethod = 'legacy';
cb.plpar = [];                                          % Disable path loss model
cb.gen_ssf_parameters;                                  % Generate large- and small-scale fading

%%
% Now, we generate the channel coefficients. The first run uses the drifting module, the second run
% disables it. Note that drifting needs significantly more computing resources. In some scenarios it
% might thus be useful to disable the feature to get quicker simulation results.   

% s.use_spherical_waves = 1;                              % Enable drifting (=spherical waves)
% c = cb.get_channels;                                    % Generate channel coefficients
% c.individual_delays = 0;                                % Remove per-antenna delays

s.use_spherical_waves = 0;                              % Disable drifting
d = cb.get_channels;                                    % Generate channel coefficients
%% freq. response and time response
if 0
    h = d.fr(57.6e6,256,1);
    h = squeeze(h(1,1,:,1));
    pdp = 20*log10(abs(ifft(h,[],1)));
    figure
    subplot(211)
    plot(linspace(0,256/57.6e6*1e9,256),pdp)
    grid on
    xlabel('Delay [ns]')
    ylabel('Power [dB]')
    xlim([-20,500])
    subplot(212)
    stem(cb.taus(1)/1e-9,10*log10(cb.pow(1))+30);hold on
    stem(cb.taus(2:21)/1e-9,10*log10(cb.pow(2:21))+30);hold on
    stem(cb.taus(22:41)/1e-9,10*log10(cb.pow(22:41))+30)
    grid on
    xlim([-20,500])
    xlabel('Delay [ns]')
    ylabel('Power [dB]')
end
%% test AoA
if 1
    cluster_of_int = cb.NumClusters;
    path_tap = zeros(cluster_of_int,Nrx);
    for kk=1:cluster_of_int
        path_tap(kk,:) = squeeze(d.coeff(:,2,kk,1));
    end

    blue = [0    0.4470    0.7410];
    red = [0.8500    0.3250    0.0980];
    yellow = [0.9290    0.6940    0.1250];
    purple = [0.4940    0.1840    0.5560];

    %  
    angle_range = linspace(-pi/2,pi/2,500);
    angle_num = 500;
    corr_pow = zeros(angle_num,cluster_of_int);
    for angle_idx = 1:angle_num
        AoA_test = angle_range(angle_idx);
        angle_vec = exp(-1j*pi*(0:(Nrx-1)).'*sin(AoA_test));
        for kk=1:cluster_of_int
            corr_pow(angle_idx,kk) = abs(path_tap(kk,:) * angle_vec)^2;
        end
    end

    range = zeros(cb.scenpar.NumSubPaths,cb.scenpar.NumClusters-1);
    for cc = 1 : (cb.scenpar.NumClusters-1)
        range(:,cc) = ((cc-1)*cb.scenpar.NumSubPaths+2) : cc*cb.scenpar.NumSubPaths+1;
    end

    figure
    for kk=1:cluster_of_int
        if sum(range(:,1) == kk)
            subplot(222)
            curvecolor = red;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor)   
        elseif sum(range(:,2) == kk)
            subplot(223)
            curvecolor = yellow;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor);    
    %     elseif sum(range(:,3) == kk)
    %         subplot(224)
    %         curvecolor = purple;
    %         gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
    %         set(gca,'Color',curvecolor)   
        else
            subplot(221)
            curvecolor = blue;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor)    
        end
    end

    for kk=1:cb.scenpar.NumClusters
    subplot(2,2,kk)
    xlabel('angle [deg]')
    ylabel('pow [dB]')
    grid on
    ylim([-40,50])
    xlim([-90,90])
    end

    AoA_wrap = 180 + cb.AoA/pi*180;
    AoA_wrap(AoA_wrap>180) = AoA_wrap(AoA_wrap>180)-360;
    figure;stem(AoA_wrap,10*log10(cb.pow)+40)
    grid on
    % xlim([0,360])
end
%% test AoD
if 1
    cluster_of_int = cb.NumClusters;
    path_tap = zeros(cluster_of_int,Ntx);
    for kk=1:cluster_of_int
        path_tap(kk,:) = squeeze(d.coeff(2,:,kk,1));
    end

    blue = [0    0.4470    0.7410];
    red = [0.8500    0.3250    0.0980];
    yellow = [0.9290    0.6940    0.1250];
    purple = [0.4940    0.1840    0.5560];

    %  
    angle_range = linspace(-pi/2,pi/2,500);
    angle_num = 500;
    corr_pow = zeros(angle_num,cluster_of_int);
    for angle_idx = 1:angle_num
        AoD_test = angle_range(angle_idx);
        angle_vec = exp(-1j*pi*(0:(Ntx-1)).'*sin(AoD_test));
        for kk=1:cluster_of_int
            corr_pow(angle_idx,kk) = abs(path_tap(kk,:) * angle_vec)^2;
        end
    end

    range = zeros(cb.scenpar.NumSubPaths,cb.scenpar.NumClusters-1);
    for cc = 1 : (cb.scenpar.NumClusters-1)
        range(:,cc) = ((cc-1)*cb.scenpar.NumSubPaths+2) : cc*cb.scenpar.NumSubPaths+1;
    end

    figure
    for kk=1:cluster_of_int
        if sum(range(:,1) == kk)
            subplot(222)
            curvecolor = red;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor)   
        elseif sum(range(:,2) == kk)
            subplot(223)
            curvecolor = yellow;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor);    
%         elseif sum(range(:,3) == kk)
%             subplot(224)
%             curvecolor = purple;
%             gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
%             set(gca,'Color',curvecolor)   
        else
            subplot(221)
            curvecolor = blue;
            gca = plot(angle_range/pi*180, 10*log10(corr_pow(:,kk)));hold on
            set(gca,'Color',curvecolor)    
        end
    end

    for kk=1:cb.scenpar.NumClusters
    subplot(2,2,kk)
    xlabel('angle [deg]')
    ylabel('pow [dB]')
    grid on
    ylim([-40,50])
    xlim([-60,60])
    end

    figure;stem(cb.AoD/pi*180,10*log10(cb.pow)+40)
    grid on
    xlim([-90,90])
end
%%
if 1
    path_tap = squeeze(d.coeff(:,:,1,1));
    angle_num = 100;
    angle_range = linspace(-pi/2,pi/2,angle_num);

    corr_pow = zeros(angle_num,angle_num);
    for AoDangle_idx = 1:angle_num
        AoD_test = angle_range(AoDangle_idx);
        AoDangle_vec = exp(-1j*pi*(0:(Ntx-1)).'*sin(AoD_test));
        for AoAangle_idx = 1:angle_num
            AoA_test = angle_range(AoAangle_idx);
            AoAangle_vec = exp(1j*pi*(0:(Nrx-1)).'*sin(AoA_test));
            corr_pow(AoDangle_idx,AoAangle_idx) = ...
                abs(AoAangle_vec' * path_tap * AoDangle_vec)^2;
        end
    end
    
    [AoDgrid, AoAgrid] = meshgrid(angle_range,angle_range);
    mesh(AoDgrid/pi*180,AoAgrid/pi*180,corr_pow)
    xlabel('AoD [deg]')
    ylabel('AoA [deg]')
end
%% Received signal using QuaDRiGa channels

path_num = (cb.scenpar.NumClusters-1)*cb.scenpar.NumSubPaths+1; % Num of paths in totl
Nr = Nrx; % Number of antenna in Rx
Nt = Ntx;
M = 64; % Length of training
MCtimes = 1; % Num of Monte Carlo Sim.
BW = 57.6e6; % IA bandiwdth
Ts = 1/BW; % Sample duration
Nb = 512; % Sample per SS burst
CFO_ppm = 0; % CFO in ppm
CFO = (28e9/1e6*CFO_ppm); % With unit Hz
eF = CFO*Ts*2*pi; % 
P = 128;
DFT = dftmtx(P);
to_est_CFO=1;
% max_ite_num = 1e3;
refine_CFO = 0; % refine CFO when coarse estimation of AoA/AoD (long long time!!!)

%-------- dictionary generation -------------
cand_num_r = 9;
cand_num_t = 33;
dict_num = cand_num_r*cand_num_t;

cand_y = zeros(M,cand_num_r*cand_num_t);
cand_angle_r = linspace(-pi*60/180,pi*60/180,cand_num_r);
AOAstep = cand_angle_r(2)-cand_angle_r(1);
cand_angle_t = linspace(-pi*60/180,pi*60/180,cand_num_t);
AODstep = cand_angle_t(2)-cand_angle_t(1);

cand_ARV_r = exp(1j*(0:Nr-1)'*pi*sin(cand_angle_r));
cand_ARV_t = exp(-1j*(0:Nt-1)'*pi*sin(cand_angle_t));

% test scenario when assuming CFO is known
phase_error_mat = kron(exp(1j*eF*Nb*(0:M-1)),exp(1j*eF*(0:P-1)'));
phase_error = reshape(phase_error_mat,M*P,1);


% ---- signal length parameters --------
ZC_root = 29; % ZC root, a coprime number with ZC length
ZC_N = 127; % ZC sequence length
seq = lteZadoffChuSeq(ZC_root,ZC_N);
                
debug_flag = 0;           

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

% weighted average
for dd=1:dict_num
index_new = [abs(Measure_mat_new(1:M-1,dd)),abs(Measure_mat_new(2:M,dd))];
CFO_weight = min(index_new,[],2);
CFO_weight_norm = CFO_weight/sum(CFO_weight);
CFO_select(:,dd) = CFO_weight_norm>1/M;
end


% Find True AoA/AoD in grid (for debug)
[~,row_true] = min(abs(cand_angle_r - (cb.AoA(1)+pi)));
[~,col_true] = min(abs(cand_angle_t - cb.AoD(1)));
index_true = (col_true-1)*cand_num_r + row_true;


% Pre-compute some vectors/matrices in FIM
% for pathindex = 1:path_num
% 
%     % Spatial response and its derivative over phi
%     arx(:,pathindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(pathindex)))/sqrt(Nr);
% 
%     % Spatial response and its derivative over theta
%     atx(:,pathindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(pathindex)))/sqrt(Nt);
% 
%     % Delay response and its derivative over tau
%     fvec(:,pathindex) = exp(-1j * (0:P-1)' * tau_samp(MCindex) / P);
% 
% end

% About CFO and its derivative
qvec = exp(1j * (0:P-1)' * eF);
for mm=1:M
    timeindex = ((mm-1)*Nb+0):((mm-1)*Nb+P-1);
    Dqvec(:,mm) = 1j * (timeindex).';
    dqvec(:,mm) = exp(1j*Nb*(mm-1)*eF)*qvec.*Dqvec(:,mm);
end


% Pre-compute some equations, vecx, vecalpha, vecbeta
% vecx is derivative of phi_0 in f(r,eta)

% Received signals (Using random symbol or ZC sequence)
symb = [seq;1]; %exp(1j*rand(P,1)*2*pi);
tau_num = 500;
delay_cand = linspace(0,300,tau_num)*1e-9/Ts*2*pi;
for tt=1:tau_num
    delay_mtx(:,tt) = DFT'*(exp(-1j * (0:P-1)' * delay_cand(tt) / P).*symb);
end

sig_rx = zeros(P*M,1);
for mm=1:M
    index = (mm-1)*P+1:mm*P;
    for ll=1:(cb.scenpar.NumClusters-1)*cb.scenpar.NumSubPaths+1
        fvec = exp(-1j * 2*pi * (0:P-1)' * d.delay(ll,1)  / (P*Ts));
        sig_rx(index) = sig_rx(index) +...
            (W(:,mm)'*squeeze(d.coeff(:,:,ll,1)) * F(:,mm))...
            *(diag(qvec)*DFT'*(fvec.*symb)) * exp(1j*Nb*(mm-1)*eF);
    end
end
%% Run Proposed Algorithm!
% -------- CFO Estimation and Delay Matching Pursuit ------------
MCindex = 1;
ss = 1;
sig_noisy = sig_rx;
sig_ave = mean(reshape(sig_noisy,P,M),2);
for tt=1:tau_num
%             if MCindex==1 & tt==6
%                 apple=1;
%             end

    % If selected delay candidate is true, sig.*conj(p_kx) is like
    % tone signal
    sig_tone = sig_ave.*conj(delay_mtx(:,tt))./abs(delay_mtx(:,tt)).^2;

%             % Method1 ML estimation of CFO; Kay's paper on freq. est. of tone
    for pp = 1:P-1
        CFO_hat(pp) = angle(sig_tone(pp+1).*conj(sig_tone(pp)));
    end
    CFO_est = sum(CFO_hat)/(P-1);

    % Method2 ML estimation of CFO; Kay's paper on freq. est. of tone
    for pp = 1:P-1
        CFO_hat(pp) = sig_tone(pp+1).*conj(sig_tone(pp));
    end
    CFO_est = angle(sum(CFO_hat)/(P-1));
    CFO_est = 0;

    % Use the ML CFO to evaluate efficacy of delay candidate
    CFO_est1 = angle(mean(exp(1j*CFO_hat)));
%             CFO_est = eF; % For debug; plug-in true CFO
    sig_deCFO = sig_ave.*exp(-1j*CFO_est*(0:P-1).');
    score(tt) = abs(sum(sig_deCFO.*conj(delay_mtx(:,tt)))/norm(delay_mtx(:,tt))^2);
end

% Plug-in Matching Pursuit index for CFO and Delay est.
[maxscore,maxindex] = max(score);
maxindex_est(ss,MCindex) = maxindex;
delay_est(ss,MCindex) = delay_cand(maxindex);

% watch score to debug
%         figure;plot(delay_cand/2/pi*Ts/1e-9,score)

%% ---------------   AoA/AoD Estimation   ----------------
sig_desymb = zeros(M,1);

% Evaluate average sample over each SS burst
for mm=1:M
    indextemp = (mm-1)*P+1:mm*P;
    sig_desymb(mm) = mean(sig_noisy(indextemp).*(conj(delay_mtx(:,maxindex))...
        ./abs(delay_mtx(:,maxindex)))).';
end

% Use true CFO (for debug)
if to_est_CFO==0
    phase_error_mat = exp(1j*eF*Nb*(0:M-1).');
    phase_error = phase_error_mat;
end

% Matching pursuit for all AoA/AoD pair in dictionary
for dd=1:dict_num

    % Debug flag, used by seting AoA=AoD = 0
%             if MCindex==26 & dd==3691
%                 apple=1;
%             end

%             sig_cand = zeros(P*M,1);
%             for mm=1:M
%                 indextemp = (mm-1)*P+1:mm*P;
%                 sig_cand(indextemp) = delay_mtx(:,maxindex)*Measure_mat_new(mm,dd);
%             end
%             sig_cand = kron(Measure_mat_new(:,dd),delay_mtx(:,maxindex));

    % use estimated CFO or perfect CFO to comp phase error

    if to_est_CFO

        sig_burst = sig_desymb.*conj(Measure_mat_new(:,dd))./abs(Measure_mat_new(:,dd));

        % adjust N/A numbers
        sig_burst(isnan(sig_burst))=0;

%                 %  ---------  Quasi-ML Method 1  ---------
%                 CFO_hat_new1 = zeros(M-1,1);
%                 for mm=1:M-1
%                     if CFO_select(mm,dd)
%                     CFO_hat_new1(mm) = angle(sig_burst(mm+1).*conj(sig_burst(mm)));
%                     end
%                 end
%                 CFO_est = sum(CFO_hat_new1)/sum(CFO_select(:,dd))/Nb;

%                 
        %  ---------  Quasi ML Method 2  ---------
        CFO_hat_new2 = zeros(M-1,1);
        for mm=1:M-1
            if CFO_select(mm,dd)
                CFO_hat_new2(mm) = sig_burst(mm+1).*conj(sig_burst(mm));
            end
        end
        CFO_est = angle(sum(CFO_hat_new2(CFO_select(:,dd))./abs(CFO_hat_new2(CFO_select(:,dd))))/sum(CFO_select(:,dd)))/Nb;

        % It seems refine CFO is necessary at beginning
        if refine_CFO
            CFO_range = linspace(CFO_est*0.7,CFO_est*1.4,10);
            score_CFO = zeros(10,1);
        for zz=1:10
%                     CFO_estold = CFO_est;
%                     phase_error_mat = exp(1j*CFO_estold*Nb*(0:M-1).');
%                     x_vec_CFO = Measure_mat_new(:,dd).*phase_error_mat;
%                     y_vec_CFO = sig_desymb;
%                     opt_alpha = 1;
%                     error_CFO_est =  y_vec_CFO - opt_alpha*x_vec_CFO;
%                     sig_deF = Measure_mat_new(:,dd).*(1j*(Nb*(0:M-1).').*exp(1j*CFO_estold*Nb*(0:M-1).'));           
%                     Q_cal = [real(sig_deF);imag(sig_deF)];
%                     tr = [real(error_CFO_est);imag(error_CFO_est)];
%                     deF = pinv(Q_cal) * tr;
%                     CFO_est = CFO_estold + deF;
            phase_est_test = exp(1j*CFO_range(zz)*Nb*(0:M-1).');
            score_CFO(zz) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_est_test));
%                     phase_error_mat = exp(1j*CFO_estold*Nb*(0:M-1).');
        end
        [~,CFO_index] = max(score_CFO);
        CFO_est = CFO_range(CFO_index);
        end

        % ----- true CFO (for debug) ---------
%                 CFO_est = eF;

        %  --------- Simpler approach of estimating CFO  (no good)---------
%                 CFO_est = angle(sum(sig_burst(find(CFO_select(:,dd)+1)))...
%                     *conj(sum(sig_burst(find(CFO_select(:,dd))))))/Nb;

        %  --------- simpler approach 2 (no good) ---------
%                 pad = 511;
%                 x_temp = sig_burst.*conj(Measure_mat_new(:,dd));
%                 CFO_est = max(abs(fft([x_temp;zeros(pad*M,1)])))*(2*pi)/((pad+1)*Nb);

        phase_error_mat = exp(1j*CFO_est*Nb*(0:M-1).');
        phase_error = phase_error_mat;
        CFO_final(dd) = CFO_est;
    end
    score_final(dd) = abs(sig_desymb'*(Measure_mat_new(:,dd).* phase_error))...
        /(Measure_mat_new(:,dd)'*Measure_mat_new(:,dd));
end
[~,bestindex_comp(MCindex)] = max(abs(score_final));

%         bestindex_comp(MCindex) = index_true;% debug. comment in main script
bestrow = floor((bestindex_comp(MCindex)-1)/cand_num_r)+1;
bestcol = bestindex_comp(MCindex)-(bestrow-1)*cand_num_r;
bestAOA(MCindex,ss) = (bestcol-1)*AOAstep-60*pi/180;
bestAOD(MCindex,ss) = (bestrow-1)*AODstep-60*pi/180;

% Plot score for debug
if 1
    figure;plot(abs(score_final));hold on;plot(index_true,abs(score_final(index_true)),'o','linewidth',2,'markersize',10);
    grid on
    legend('scores','true AoA/AoD pair')

    CFO_est_range = linspace(0.0010,0.0020,200);
    for xx=1:200
        CFO_est = CFO_est_range(xx) ;
    phase_error_mat = exp(1j*CFO_est*Nb*(0:M-1).');
%             phase_error_mat = exp(1j*eF*Nb*(0:M-1).');
    phase_error = phase_error_mat;
    score_test_CFO(xx) = abs(sig_desymb'*(Measure_mat_new(:,dd).*phase_error));
    end
    figure
    plot(CFO_est_range,score_test_CFO)

end


