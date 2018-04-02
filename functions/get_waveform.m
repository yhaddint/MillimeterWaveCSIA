function [ sig_unitpow ] = get_waveform( QAM_level, sig_length, upsam, M )
%GET_WAVEFORM Summary of this function goes here
%   Detailed explanation goes here
    
    sig_unitpow = zeros(sig_length, M);
    L = sig_length;
    for mm=1:M
        symbols = fix(L*2/upsam);   
        hmod = modem.qammod('M', QAM_level, 'InputType', 'integer');
        hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
        hpulse = design(hdesign);
        data = randi(QAM_level,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,upsam);
        temp_data = conv(data,hpulse.Numerator);
        sig = temp_data(end-L+1-1e3:end-1e3)./sqrt(temp_data(end-L+1-1e3:end-1e3)'...
            *temp_data(end-L+1-1e3:end-1e3)/L);
        sig_pow_dBm = -5;
        R = 50;
        sig_pow_V = sqrt(10^(sig_pow_dBm/10)*1e-3*R);
        sig_unitpow(:,mm) = sig * sig_pow_V;
    end
end

