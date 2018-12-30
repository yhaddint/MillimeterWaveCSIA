% This script evaluates beam pattern of different beamformer used in IA

clear;clc;
% --- codebook construction -----
Nt = 64;
A_stopband = 30; % attenuation outside mainlobe (dB)
MCtimes = 5e2;
PSquan_range = 2:10;
SINR_noncomp = zeros(length(PSquan_range),MCtimes);
SINR_ncomp = zeros(length(PSquan_range),MCtimes);
for MCindex = 1:MCtimes
    theta1 = rand*10+30;
    theta2 = -theta1;

    vec1 = get_FSM_KW_codebook( theta1/180*pi, 5/180*pi, Nt, A_stopband);
    vec2 = get_FSM_KW_codebook( theta2/180*pi, 5/180*pi, Nt, A_stopband);
%     vec1 = exp(1j*pi*(0:Nt-1).'*sin(theta1/180*pi));
%     vec2 = exp(1j*pi*(0:Nt-1).'*sin(theta2/180*pi));
    % vec = vec1 + vec2;
    sigman2 = 0.01;

    % vec1_norm = vec1./norm(vec);
    % vec2_norm = vec2./norm(vec);

        % ------- test hybrid approximation-------
    for pp=1:length(PSquan_range)
        phase_bits = PSquan_range(pp);
        vec1_hf = get_hybrid_approx( vec1, phase_bits );
        vec2_hf = get_hybrid_approx( vec2, phase_bits );

        vec1_norm = vec1_hf./norm(vec1_hf);
        vec2_norm = vec2_hf./norm(vec2_hf);

        atx1 = exp(1j*pi*(0:Nt-1).'*sin(theta1/180*pi));
        atx2 = exp(1j*pi*(0:Nt-1).'*sin(theta2/180*pi));

        H_eff = [atx1'*vec1_norm, atx1'*vec2_norm; atx2'*vec1_norm, atx2'*vec2_norm];
        ZF_mtx = pinv(H_eff);
        ZF_mtx_norm = ZF_mtx./norm(ZF_mtx,'fro');

        H_eff_nocomp = H_eff*[sqrt(0.5),0;0,sqrt(0.5)];
        H_eff_ZF = H_eff*ZF_mtx_norm;

        SINR_noncomp(pp,MCindex) = (abs(H_eff_nocomp(1,1))^2/(abs(H_eff_nocomp(1,2))^2+sigman2)+abs(H_eff_nocomp(2,2))^2/(abs(H_eff_nocomp(2,1))^2+sigman2))/2;
        SINR_comp(pp,MCindex) = (abs(H_eff_ZF(1,1))^2/(abs(H_eff_ZF(1,2))^2+sigman2)+abs(H_eff_ZF(2,2))^2/(abs(H_eff_ZF(2,1))^2+sigman2))/2;
    end
end
%%
figure
plot(PSquan_range,10*log10(mean(SINR_noncomp,2)),'-o','linewidth',2);hold on
plot(PSquan_range,10*log10(mean(SINR_comp,2)),'-o','linewidth',2);hold on
grid on
xlabel('PS Quantization Bits')
ylabel('SINR (dB)')
legend('PS Error Non-Aware','PS Error Aware');


