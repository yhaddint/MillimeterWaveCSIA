function [ v_approx ] = get_hybrid_approx( v, phase_bits )
%GET_HYBRID_APPROX Summary of this function goes here
%   Detailed explanation goes here
    v_abs = abs(v);
    v_angle = phase(v);
    vmax  = max(v_abs);
    angle1 = v_angle - acos(v_abs./(2*vmax));
    angle2 = v_angle + acos(v_abs./(2*vmax));
    angle1_quan = get_phase_quan( wrapTo2Pi(angle1)-pi, phase_bits )+pi;
    angle2_quan = get_phase_quan( wrapTo2Pi(angle2)-pi, phase_bits )+pi;;
    v_approx = vmax*(exp(1j*angle1_quan)+exp(1j*angle2_quan));

end

