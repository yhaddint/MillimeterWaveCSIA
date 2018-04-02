function [ F_RF, F_BB ] = get_SC_hybrid_beamformer( F_opt, Nt, Nrf, RFCB_size )
%GET_SC_HYBRID_BEAMFORMER Summary of this function goes here
%   Detailed explanation goes here
    At = zeros(Nt,RFCB_size);
    angle_candidate = linspace(-pi/2,pi/2,RFCB_size);
    for aa = 1:RFCB_size
        At(:,aa) = exp(1j*pi*(0:Nt-1).'*sin(angle_candidate(aa)));
    end

    F_res = F_opt;
    F_RF = zeros(Nt,Nrf);
    for rr=1:Nrf
        Psi = At'*F_res;
        [~,aindex] = max(diag(Psi*Psi'));
        F_RF(:,rr) = At(:,aindex);
        F_BB = pinv(F_RF(:,1:rr)) * F_opt;
        F_res0 = F_opt - F_RF(:,1:rr) * F_BB;
        F_res = F_res0./norm(F_res0,'fro');
    end

end

