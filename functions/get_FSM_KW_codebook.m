function [ vec ] = get_FSM_KW_codebook( phi_c, beamwidth, N, A)
%GET_FSM_KW_CODEBOOK Summary of this function goes here
%   Design steering vector based on Fourier Series Method with Kaiser
%   windowing approach
%   [ vec ] = get_FSM_KW_codebook( phi_c, beamwidth, N, A)
%   phi_c: center of steered beam (rad)
%   beamwidth: desired beamwidth (rad)
%   N: number of elements in array
%   A: desired attenuation outside mainlobe
    
    % ------- Initialization ------------
    if mod(N,2)==1
        odd_element = 1;
        M = (N-1)/2;
        K_win = zeros(M+1,1);
        steer_vec = zeros(M+1,1);
    else
        M = N/2;
        odd_element = 0;
        K_win = zeros(M,1);
        steer_vec_temp = zeros(M,1);
    end
    vec = zeros(N,1);
    phi_b = beamwidth/2;
    
    % ---------- Parameters in Design (D) ---------
    if A>21
        D = (A-7.95)/14.36;
    else
        D = 0.992;
    end
    
    % ---------- Parameters in Design (gamma) ---------
    if A>= 50
        gamma = 0.11*(A-8.7);
    elseif A<=21
        gamma = 0;
    else
        gamma = 0.58*(A-21)^0.4+0.079*(A-21);
    end
    
    % ------mapping from phi domain to Psi domain --------
    d_Psi = (2*pi*D)/(N-1);
    Psi_0 = -pi*sin(phi_c)*cos(phi_b/2);
    Psi_p = pi*cos(phi_c)*sin(phi_b/2);
    Psi_b = Psi_p+d_Psi;
    
    % ---------- Kaiser Window ---------
    if odd_element
        for mm=1:M+1
            m = mm-1;
            K_win(mm) = besseli(0,gamma*sqrt(1-m^2/M^2))/besseli(0,gamma);

        end
        K_win_symm = [flipud(K_win);K_win(2:end)];
    else
        for mm=1:M
            m = mm;
            K_win(mm) = besseli(0,gamma*sqrt(1-m^2/M^2))/besseli(0,gamma);

        end
        K_win_symm = [flipud(K_win);K_win];
    end

    % ---------- Steering Vector ---------
    if odd_element
        for mm=1:N
            m = mm-M-1;
            if m==0
                m=1e-5;
            end
            steer_vec(mm) = exp(-1j*m*Psi_0)*sin(Psi_b*m)/(pi*m);
        end
    else
        for mm=1:M
            m = mm;
            steer_vec_temp(mm) = exp(-1j*(m-0.5)*Psi_0)*sin(Psi_b*(m-0.5))/(pi*(m-0.5));
        end
        steer_vec = [flipud(conj(steer_vec_temp));steer_vec_temp];
    end
    
    % ---------- Apply Window on Steer Vector ---------
    vec = K_win_symm.*steer_vec;
    

end

