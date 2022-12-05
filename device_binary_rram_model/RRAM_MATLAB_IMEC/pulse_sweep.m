%clear all;
close all;
clc;

% Input excel details
ip_file_name = 'Test/device_binary_rram_model/input_data/memory_max_temperatures.csv';
ip_sheet = 'temperature';
% m3d_liquid_sheet = 'm3d_liquid';
% tsv3d_air_sheet = 'tsv3d_air';
% tsv3d_liquid_sheet = 'tsv3d_liquid';
%emib_air_sheet = 'emib_air';
%soc_sheet = '3d_soc+';
    
%% Inputs
SELECT= 0;
    % 1: RESET
    % 2: SET
    % else: retention
FIX_GAP_MODE= 0;
    % 0: gap changeable
    % 1: fixed gap
DC_T_MODE= 0;
    % 0: transient temperature
    % 1: DC temperature
NO_PLOT_T= 0;
    % 0: plot temperature
    % 1: don't plot temperature


%% Operation parameters
if SELECT==1 
    %% RESET
    gap0= 0.1e-9;
    num_voltage= 5;
    V_amplitude= linspace(-1.4, -1.8, num_voltage);
    t_final=        1e+2;
    num_tsteps=     1e+5;
    min_unit_time=  1e-25;
    T0= (273.15+25)*ones(num_voltage,1);
elseif SELECT==2 
    %% SET
    gap0= 0.74e-9; 
    num_voltage= 5;
    V_amplitude= linspace(1.6, 2, num_voltage);
    t_final=        1e+2;
    num_tsteps=     1e+5;
    min_unit_time=  1e-25;
    T0= (273.15+25)*ones(num_voltage,1);
else 
    %% Retention
    gap0= 0.3e-9;

    T0= readmatrix(ip_file_name,'Range','A1:A1');

    num_voltage= numel(T0);
    V_amplitude= zeros(numel(T0),1);
    t_final=        1e+9;
    num_tsteps=     1e+5;
    min_unit_time=  t_final/num_tsteps;
end
ini_unit_time=  t_final/num_tsteps;


%% Other parameters/Initialization
Rth=715e3;  % Thermal resistance (K/W)
gap_min= 0.1e-9;
gap_max= 2e-9;
cycle=  1;
Vread=  0.1;
sigma_spatial=  0;
sigma_temporal= 0; 
v0= 300 + sigma_spatial*randn;
unit_time=      ini_unit_time*ones(num_voltage, 1);
[gap, T]=       deal( zeros(num_tsteps+1, 1) );
[pulse_time, R]=deal( zeros(num_tsteps+1, num_voltage) );
legend_labels=  cell(1,num_voltage);
VAR1=   0.99;   % (0,1),    Timestep control parameter in SET
VAR2=   0.1;    % (0,1),    Timestep control parameter in RESET
VAR3=   2;      % (1,inf),  Division for timestep reduction

%%
for jj=1:cycle
    for ii=1:num_voltage
		gap(1)=	gap0;
		I=    state_current(gap0, Vread);
		R(1,ii)=    Vread/I;
        T(1,ii)=    T0(ii);

        if (FIX_GAP_MODE)
            pulse_time(1,ii)=  unit_time(ii);
            for kk=2:num_tsteps+1
                pulse_time(kk,ii)= pulse_time(kk-1,ii) + unit_time(ii);
                I=     state_current( gap(kk-1), V_amplitude(ii) );
                if (DC_T_MODE)
                    T(kk,ii)=   T0(ii) + abs(V_amplitude(ii)*I*Rth);
                else
                    T(kk,ii)=   T_tran( T(kk-1,ii), T0(ii), V_amplitude(ii), I, unit_time(ii) );
                end
                gap(kk)=   gap(kk-1);
                I=     state_current( gap(kk), Vread );
                R(kk,ii)=  Vread/I;
            end
        else
            %% Auto timestepping
            while (1)
                pulse_time(1,ii)=  unit_time(ii);
                pulse_time(2,ii)= pulse_time(1,ii) + unit_time(ii);
                I=     state_current( gap(1), V_amplitude(ii) );
                if (DC_T_MODE)
                    T(2,ii)=   T0(ii) + abs(V_amplitude(ii)*I*Rth);
                else
                    T(2,ii)=   T_tran( T(1,ii), T0(ii), V_amplitude(ii), I, unit_time(ii) );
                end
                gap(2)=   CF_growth( gap(1), V_amplitude(ii), unit_time(ii), T(2,ii), v0, sigma_temporal, gap_min, gap_max );
                I=     state_current( gap(2), Vread );
                R(2,ii)=  Vread/I;
                %% Targets
                if (sign(V_amplitude(ii))==1 && gap(2)-gap_min>VAR1*(gap(1)-gap_min))...
                    || (sign(V_amplitude(ii))==-1 && gap(2)-gap(1)<VAR2*(gap_max-gap(1)))...
                    || unit_time(ii)<min_unit_time
                    break;
                end
                unit_time(ii)= unit_time(ii)/VAR3;
            end
            %%
            for kk=3:num_tsteps+1
                pulse_time(kk,ii)= pulse_time(kk-1,ii) + unit_time(ii);
                I=     state_current( gap(kk-1), V_amplitude(ii) );
                if (DC_T_MODE)
                    T(kk,ii)=   T0(ii) + abs(V_amplitude(ii)*I*Rth);
                else
                    T(kk,ii)=   T_tran( T(kk-1,ii), T0(ii), V_amplitude(ii), I, unit_time(ii) );
                end
                gap(kk)=   CF_growth( gap(kk-1), V_amplitude(ii), unit_time(ii), T(kk,ii), v0, sigma_temporal, gap_min, gap_max );
                I=     state_current( gap(kk), Vread );
                R(kk,ii)=  Vread/I;
            end
        end
        len = length(R);
        drift_coefficient = (R(len)-R(1))/R(1);
        fprintf("Drift Coefficient: %12.3e\n\n", drift_coefficient);
        
        %% always output min and max temperature
        switch ii
            case 1
                file_name_ret = ['../Test/device_binary_rram_model/output_results/', 'retention', '.txt'];
                file_name_temp_max = ['../Test/device_binary_rram_model/output_results/', 'max_temps_mem.txt'];
                file_name_time = ['../Test/device_binary_rram_model/output_results/', 'time.txt'];
                file_name_drift = ['../Test/device_binary_rram_model/output_results/', 'drift_coefficient', '.txt'];
            otherwise
                disp('No retention cases to analyze.');
        end

        fid=fopen(file_name_ret,'w+');
        if fid == -1
            error('Author:Function:OpenFile', 'Cannot open file: %s', file_name_ret);
        end
        fprintf(fid,'%g\r\n',R(:,ii));
        fclose(fid);
        
        fid=fopen(file_name_temp_max,'w+');
        if fid == -1
            error('Author:Function:OpenFile', 'Cannot open file: %s', file_name_temp_max);
        end
        fprintf(fid,'%g\r\n',(T0(ii)-273.15));
        fclose(fid);
        
        fid=fopen(file_name_time,'w+');
        if fid == -1
            error('Author:Function:OpenFile', 'Cannot open file: %s', file_name_time);
        end
        fprintf(fid,'%g\r\n',pulse_time(:,ii));
        fclose(fid);
        
        fid=fopen(file_name_drift,'w+');
        if fid == -1
            error('Author:Function:OpenFile', 'Cannot open file: %s', file_name_drift);
        end
        fprintf(fid,'%g\r\n',drift_coefficient);
        fclose(fid);
            
        %% Plots
        CM= prism(num_voltage);
        
        figure(1);
        loglog(pulse_time(:,ii), R(:,ii), 'o-', 'color', CM(ii,:));
        xlabel('Time (sec)');
        ylabel('Resistance (Ohm)');
        hold on;
        
        if (~NO_PLOT_T)
            figure(2);
            semilogx(pulse_time(:,ii), T(:,ii)-273.15, 'color', CM(ii,:));
            xlabel('Time (sec)');
            ylabel('Temperature (^{o}C)');
            hold on;
        end
        
        % Legends
        if SELECT==1 || SELECT==2
            legend_labels{ii}= [num2str(V_amplitude(ii)) 'V'];
        else
            legend_labels{ii}= [num2str(T0(ii)-273.15) '^{o}C'];
        end
        
    end
end


%% Formating
for ff=1:numel(findobj('Type','figure'))
    figure(ff);
    legend(legend_labels, 'location', 'northwest');
    grid on;
end
set(findall(findall(0,'Type','figure'),'type','axes'),'fontsize',14);
set(findall(findall(0,'Type','figure'),'type','text'),'fontsize',14);
set(findobj(findall(0,'Type','figure'),'type','line'),'LineWidth',2);

