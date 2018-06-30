function [ output ] = get_portion( input, portion )
%GET_PORTION Summary of this function goes here
%   Detailed explanation goes here
sortinput = sort(input,'ascend');
output = sortinput(1:floor(length(input)*portion));
end

