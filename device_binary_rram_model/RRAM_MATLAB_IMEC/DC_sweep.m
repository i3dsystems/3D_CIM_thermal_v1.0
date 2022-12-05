clear all;
close all;
clc;

%% Operation parameters
select= 2;
if select==1 % RESET
    gap0=   1e-9;   % Intial gap, 1.1589e-9
    V_stop= -2;     % DC sweep stop voltage, use postive voltage for set, negative voltage for reset
elseif select==2 % SET 
    gap0=   1.7e-9; % Inital gap, 17e-10 m
    V_stop= 1.5;    % DC sweep stop voltage, use postive voltage for set, negative voltage for reset
end

%% Other parameters/Initialization
T0= 273+25;
Rth=715e3; % Thermal resistance (K/W)

gap_min=    0.1e-9;
gap_max=    2e-9;

cycle=1;
ramp_rate=      0.02; % V/sec
voltage_step=   0.1;
time_step=      abs(voltage_step)/ramp_rate; 

if V_stop<0
    voltage_step=   -voltage_step;
end
num_voltage_step=   V_stop/voltage_step; % Number of voltage steps

internal_time_step= 1e-6;
internal_step_num=  time_step/internal_time_step; % More resolution in time step

sigma_spatial=  0; 
sigma_temporal= 0; 
v0= 300 + sigma_spatial*randn; 

voltage=    zeros(num_voltage_step,1);
gap_final=  zeros(num_voltage_step,1);
current_final=zeros(num_voltage_step,1);
gap=zeros(internal_step_num,1);

%%
for jj=1:cycle
    for ii=1:num_voltage_step % Voltage sweep

        voltage(ii)= voltage_step*ii;  % Voltage steps

         for gg= 2:internal_step_num % Internal loop at each voltage step
              if ii==1 % This way we avoid ii=1
                  gap(1)=   gap0;
              else
                  gap(1)=   gap_final(ii-1);
              end;
              I=   state_current( gap(gg-1), voltage(ii) );
              T=      T0 + abs(voltage(ii)*I*Rth);
              gap(gg)=    CF_growth( gap(gg-1), voltage(ii), internal_time_step, T, v0, sigma_temporal, gap_min, gap_max );
         end
         
         %  Record at end of each voltage step
         gap_final(ii)=         gap(gg);
         current_final(ii)=     state_current( gap_final(ii,1), voltage(ii) ); 
         %[current_final(ii),current_tunneling(ii), current_ohmic(ii)]=state_current(gap_final(ii),voltage(ii)); 

    end

    figure(1);
    semilogy(voltage, abs(current_final));
    hold on;
    % plot(voltage,abs(current_tunneling),'g');
    % plot(voltage,abs(current_ohmic),'r');
    xlabel('Voltage (V)');
    ylabel('Current (A)');

    figure(2);
    gap_final(:,1)= gap_final(:,1)*1e9;
    plot(voltage, gap_final);
    xlabel('voltage (V)');
    ylabel('gap distance (nm)');
    hold on;

end

%% IMEC data
figure(1);
if select==1 % RESET
    dataV=[-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1 -1.1 -1.2 -1.3 -1.4 -1.5 -1.6 -1.7 -1.8 -1.9 -2 -2 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1 -1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0];
    dataI=[2.13e-07 4.86e-07 8.65e-07 1.33e-06 1.92e-06 2.95e-06 4.01e-06 4.36e-06 5.90e-06 2.01e-06 1.10e-06 1.23e-06 1.26e-06 1.59e-06 2.05e-06 2.86e-06 3.99e-06 4.59e-06 4.92e-06 4.86e-06 5.45e-06 4.59e-06 4.01e-06 3.49e-06 3.02e-06 2.58e-06 2.18e-06 1.80e-06 1.47e-06 1.19e-06 9.44e-07 7.46e-07 5.86e-07 4.54e-07 3.42e-07 2.48e-07 1.62e-07 6.52e-08 2.99e-08 1.13e-08 5.23e-11];
    hold on;
    semilogy(dataV,dataI,'r');
    legend('Model','IMEC data','location','southwest');
end
if select==2 % SET
    dataV=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
    dataI=[2.61e-08 6.75e-08 1.13e-07 1.77e-07 2.60e-07 3.56e-07 4.57e-07 4.44e-07 6.91e-07 8.28e-07 9.35e-07 1.10e-06 1.24e-06 5.66e-06 9.13e-06 9.32e-06 9.39e-06 9.45e-06 9.50e-06 9.54e-06 9.56e-06 9.51e-06 9.47e-06 9.41e-06 9.34e-06 9.26e-06 9.17e-06 9.05e-06 8.89e-06 8.63e-06 8.23e-06 7.69e-06 7.03e-06 6.26e-06 5.37e-06 4.49e-06 3.56e-06 2.60e-06 1.70e-06 8.28e-07 2.99e-10];
    hold on;
    semilogy(dataV,dataI,'r');
    legend('Model','IMEC data','location','southeast');
end

%%
set(findall(findall(0,'Type','figure'),'type','axes'),'fontsize',14);
set(findall(findall(0,'Type','figure'),'type','text'),'fontsize',14);
set(findobj(findall(0,'Type','figure'),'type','line'),'LineWidth',2);

