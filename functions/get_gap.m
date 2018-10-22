function [gap] = get_gap(inputvec)
%GET_GAP Summary of this function goes here
%   Detailed explanation goes here
ele_num = length(inputvec);
gap = 2*pi;
for i1=1:ele_num
    for i2=i1+1:ele_num
        gap = min(abs(inputvec(i1) - inputvec(i2)),gap);
    end
end
end

