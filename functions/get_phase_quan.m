function [ phase_out ] = get_phase_quan( phase_in, phase_bits )
%GET_PHASE_QUAN Summary of this function goes here
%   Detailed explanation goes here
    
    PS_bits_num = phase_bits;
    PS_scale = 2^(PS_bits_num-1);
    phase_out = (round(phase_in/pi*PS_scale+0.5)-0.5)/PS_scale*pi;
end

