function [ BF0 ] = get_IA_BF( ant_num, M, type )
%GET_IA_BF Summary of this function goes here
%   Detailed explanation goes here
    switch type
        case 'PN'
            BF = (randi(2,ant_num,M)*2-3) + 1j * (randi(2,ant_num, M)*2-3);
        case 'directional'
            angle_sweep = linspace(-pi/2,pi/2,M);
            for mm=1:ant_num
                BF(:,mm) = exp(1j*pi*(0:ant_num01).'*sin(angle_sweep(mm)));
            end
        otherwise
            fprintf('Error in IA beamformer type');
    end
    BF0 = BF./norm(BF,'fro')*sqrt(M);
    

end

