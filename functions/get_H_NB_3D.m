function [ H_NB ] = get_H_NB_3D( raygain, rayAOA_az, rayAOA_el, rayAOD_az, rayAOD_el, cluster_num, ray_num, Nt_az, Nt_el, Nr_az, Nr_el )
%GET_H_FREQ Summary of this function goes here
%   This is the function used for most script
%   Channel model is H = sum
%   alpha*exp(-j*2*pi*n*tau_0/(Nfft*Ts))*ar(phi)*at(theta)
%   H_freq = get_H_freq2( raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr )
%   IP: raygain (cluster_num by ray_num) matrix with complex gain for each path
%   IP: raydelay (cluster_num by ray_num) matrix with prop delay for each path
%   IP: rayAOA (cluster_num by ray_num) matrix with AOA for each path
%   IP: rayAOD_az (cluster_num by ray_num) matrix with AOD for each path (azimuth)
%   IP: rayAOD_el (cluster_num by ray_num) matrix with AOD for each path (elevation)
%   IP: cluster_num is scaler for number of multipath cluster
%   IP: ray_num is scaler for number of rays in each cluster
%   IP: Nt_az is number of antenna in az of (UPA) in transmitter
%   IP: Nt_el is number of antenna in el of (UPA) in transmitter
%   IP: Nr is number of antenna (ULA) in receiver
    
    Nt = Nt_az*Nt_el;
    Nr = Nr_az*Nr_el;
    
    H_NB = zeros(Nr,Nt);
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            
            phi_az = rayAOA_az(cluster_index, ray_index);
            phi_el = rayAOA_el(cluster_index, ray_index);
            
            arx_az = exp(1j*(0:Nr_az-1)'*pi*sin(phi_az))/sqrt(Nr_az);
            arx_el = exp(1j*(0:Nr_el-1)'*pi*sin(phi_el))/sqrt(Nr_el);
            
            arx = kron( arx_el, arx_az );
            
            theta_az = rayAOD_az(cluster_index, ray_index);
            theta_el = rayAOD_el(cluster_index, ray_index);
            
            atx_az = exp(1j*(0:Nt_az-1)'*pi*sin(theta_az))/sqrt(Nt_az);   
            atx_el = exp(1j*(0:Nt_el-1)'*pi*sin(theta_el))/sqrt(Nt_el); 
            
            atx = kron( atx_el, atx_az );
            
            H_NB = H_NB + raygain(cluster_index, ray_index)*arx*atx';

        end
    end

end