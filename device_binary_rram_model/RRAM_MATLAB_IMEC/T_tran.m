function T = T_tran( T, T0, V, I, dt )

% rho =   9.8e3;    % kg/m^3
% C_rho = 120;    % J/kg-K
% volume = (60e-9)^3; % m^3
% k =     1.1;    % W/m-K
% % Kth = k*((20e-9)^2*2)/7.5e-9;
% Kth =   k*1000e-9;
% Cth =   C_rho*rho*volume; % J/K
% Tau_th = Cth/Kth;   % sec

delay_scale=1;
T_scale= 1;
Cth =   3.1825e-16*delay_scale/T_scale; % J/K
Tau_th = 2.3e-10*delay_scale;           % sec

%T = T + dt*(abs(V*I)/Cth - (T-T0)/Tau_th); % explicit
T = (T + dt*(abs(V*I)/Cth + T0/Tau_th))/(1+dt/Tau_th);  % implicit

end