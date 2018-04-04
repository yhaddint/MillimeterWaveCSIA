function [ BF0 ] = get_IA_BF( ant_num, M, type )
%GET_IA_BF Summary of this function goes here
%   Detailed explanation goes here
    switch type
        
        % Quasi-Omni PN beams
        case 'PN'
            BF = (randi(2,ant_num,M)*2-3) + 1j * (randi(2,ant_num, M)*2-3);
        
        % Sector beam codebook design via LS 
        case 'sector'
            beam_width = pi/M;
            angle_sweep = (pi/2)*linspace(-1+(1/M),1-(1/M),M); % Sector beams with center evenly in -pi/2 to pi/2
            
            
            angle_range = 90;
            pointing_range = angle_sweep;
            desired_pattern = zeros(angle_range*2+1,M);
            for mm=1:M
                center = fix(pointing_range(mm)/pi*180 + angle_range) + 1;
                left_boundary = center - floor((beam_width/pi*180-1)/2);
                right_boundary = center + floor((beam_width/pi*180-1)/2);
                actual_width = right_boundary - left_boundary + 1;
                index_desired = (left_boundary:right_boundary);
                desired_pattern(index_desired,mm) = ones(actual_width,1);
            end
            for kk = 1:(angle_range*2+1)
                FF(:,kk) = exp(1j*pi*(0:ant_num-1).'*sin((kk - angle_range -1 )/180*pi));
            end

            for mm=1:M
                 BF(:,mm) = pinv(FF')*desired_pattern(:,mm);
            end
            
        % Narrow (104/N deg) beam codebook design
        case 'directional'
            angle_sweep = (pi/2)*linspace(-1+(1/M),1-(1/M),M); % Sector beams with center evenly in -pi/2 to pi/2
            for mm=1:ant_num
                BF(:,mm) = exp(1j*pi*(0:ant_num-1).'*sin(angle_sweep(mm)));
            end
            
        otherwise
            fprintf('Error in IA beamformer type');
    end
    
    % Normalization for unit tranmission (reception) power
    BF0 = BF./norm(BF,'fro')*sqrt(M);
    

end

