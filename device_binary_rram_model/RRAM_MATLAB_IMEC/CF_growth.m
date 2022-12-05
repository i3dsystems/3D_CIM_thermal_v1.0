function gap= CF_growth( gap, voltage, tstep, T, v0, sigma_temporal, gap_min, gap_max )

    q=  1.6E-19;
    kb= 1.3806503E-23;
    kT= kb*T;
    a0= 0.25e-9;
    L=  5E-9; % Oxide thickness (m)
    %Ea=1.5; % Acitivated energy of Vo in HfOx
    Ea_set= 1.501;  % eV
    Ea_reset= 1.5;  % eV
    %F_min=1.4E9;

    if gap>=gap_min && gap<=gap_max
        
        enh = 16.5-1.25*(gap/1E-9)^3; % Field enhancment factor
        
        %if enh*abs(voltage)/L<F_min;
        %       enh=0;
        %end
        
        gamma_V=    enh*a0/L;
        delta_gap=  v0*tstep*( exp(-q*Ea_set/kT)*exp(gamma_V*q*voltage/kT) - exp(-q*Ea_reset/kT)*exp(-gamma_V*q*voltage/kT) )/2;
        gap=        gap - delta_gap*(1+randn*sigma_temporal);

        if gap<gap_min
            gap= gap_min;
        elseif gap>gap_max
            gap= gap_max;
        end

    end
end