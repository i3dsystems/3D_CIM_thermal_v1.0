function startup(var, val)
    % Add files to path
    if isempty(strfind(path, './device_binary_rram_model);'))
        addpath(genpath('./device_binary_rram_model'));
    end
    
    % Add files to path
    if isempty(strfind(path, './TSV_stacked_3D_v4/air_cooling;'))
        addpath(genpath('./TSV_stacked_3D_v4/air_cooling'));
    end
    
    % Calling testbench
    testbench_rram_mem_on_logic_air();
    
    % Save files in results
    copyfile('device_binary_rram_model/input_data/memory_max_temperatures.csv','TSV_stacked_3D_v4/air_cooling/result')
    % Save thermal results
    if var == 0
        filename = 'TSV_stacked_3D_v4/air_cooling/result/0.csv';
    else
        filename = 'TSV_stacked_3D_v4/air_cooling/result/' + string(var) + '_' + string(val) + '.csv';
    end
    movefile('TSV_stacked_3D_v4/air_cooling/result/memory_max_temperatures.csv', filename)
    
    
    % Call pulse sweep
    fprintf("\n\n-----Calling Binary RRAM Model...----\n\n");
    pulse_sweep()
    
    % Save files in results
    copyfile('device_binary_rram_model/output_results/max_temps_mem.txt','device_binary_rram_model/output_results/results');
    copyfile('device_binary_rram_model/output_results/retention.txt','device_binary_rram_model/output_results/results');
    copyfile('device_binary_rram_model/output_results/time.txt','device_binary_rram_model/output_results/results');
    copyfile('device_binary_rram_model/output_results/drift_coefficient.txt','device_binary_rram_model/output_results/results');
    
    if var == 0
        filename_ret = 'device_binary_rram_model/output_results/results/retention_0.txt';
        filename_max = 'device_binary_rram_model/output_results/results/max_temps_mem_0.txt';
        filename_t = 'device_binary_rram_model/output_results/results/time_0.txt';
        filename_dc = 'device_binary_rram_model/output_results/results/drift_coefficient_0.txt';
    else
        filename_ret = 'device_binary_rram_model/output_results/results/retention_' + string(var) + '_' + string(val) + '.txt';
        filename_max = 'device_binary_rram_model/output_results/results/max_temps_mem_' + string(var) + '_' + string(val)+'.txt';
        filename_t = 'device_binary_rram_model/output_results/results/time_' + string(var) + '_' + string(val)+'.txt';
        filename_dc = 'device_binary_rram_model/output_results/results/drift_coefficient_' + string(var) + '_' + string(val)+'.txt';
        
    end
    movefile('device_binary_rram_model/output_results/results/max_temps_mem.txt', filename_max);
    movefile('device_binary_rram_model/output_results/results/retention.txt', filename_ret);
    movefile('device_binary_rram_model/output_results/results/time.txt', filename_t);
    movefile('device_binary_rram_model/output_results/results/drift_coefficient.txt', filename_dc);
    
    fprintf("-----FIN-----\n\n");
    
    end