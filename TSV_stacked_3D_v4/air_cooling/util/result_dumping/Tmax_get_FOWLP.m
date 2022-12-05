function T = Tmax_get_FOWLP(x, system, chip, Layer, draw, t, cell_count, sim_case)
     
    T = zeros(system.chip.N*2,1);
    pointer = 1;
    temp_matrix = zeros(1,(system.chip.N*3));
    max_temp = 0;
    
    for i = 1:1:system.chip.N
        const = chip(i).Nx*chip(i).Ny;
%         index_start = system.Nhsp + chip(i).Nx*chip(i).Ny * sum(Layer.N(1:chip(i).layer.id-1));            
        index_start = system.Nhsp;
        for l= 1 : 1 : i-1
            index_start = index_start + chip(l).Nx*chip(l).Ny.*(Layer.N(l)-1);
        end
        id = const +index_start;
        left = 1+id;
        right = const+id;
        
        T(pointer) = max(x(left:right))-273;
        
        T(pointer+system.chip.N) = min(x(left:right))-273;
        if draw.displayT == 1
            fprintf('chip%i Min ~ Max: %.2f ~ %.2f\n', ...
                          i,  T(pointer+system.chip.N), T(pointer));
        end
        
        if strcmp('Memory Tier', chip(i).name)
            max_temp = T(pointer) + 273.15;
        end
        %% always output min and max temperature
%         file_name_max = ['./result/', num2str(chip(1).Xsize), '_', num2str(chip(1).power), '_', num2str(chip(1).layer.die(1)), '_', 'chip', num2str(i), '_max.txt'];
%         file_name_min = ['./result/', num2str(chip(1).Xsize), '_', num2str(chip(1).power), '_', num2str(chip(1).layer.die(1)), '_', 'chip', num2str(i), '_min.txt'];
%         if t == 0
%             fid=fopen(file_name_max,'w+');
%             fprintf(fid,'%g\r\n',T(pointer));
%             fclose(fid);
%             fid=fopen(file_name_min,'w+');
%             fprintf(fid,'%g\r\n',T(pointer+system.chip.N));
%             fclose(fid);
%         else
%             fid=fopen(file_name_max,'a+');
%             fprintf(fid,'%g\r\n',T(pointer));
%             fclose(fid);
%             fid=fopen(file_name_min,'a+');
%             fprintf(fid,'%g\r\n',T(pointer+system.chip.N));
%             fclose(fid);
%         end
        if draw.write == 1
%             file_name = ['./result/', num2str(chip(1).Xsize), '_', num2str(chip(1).power), '_', num2str(chip(1).layer.die(1)), '_', 'chip', num2str(i), '.txt'];
%             file_name_avg = '../results_air.xlsx'; % file to write average temperatures to excel sheet
%             file_name_avg = '../results_simple_case_air.xlsx'; % file to write average temperatures to excel sheet
            file_name_avg = ['../results_worst_case_air_case_', num2str(sim_case), '.xlsx']; % file to write average temperatures to excel sheet
            if t == 0
                fid=fopen(file_name,'w+');
                x_temp = x(left : right)-273;
                fprintf(fid,'%g\r\n',x_temp);
                fclose(fid);                 
            else % writing average temperatures to excel sheet
%                 fid=fopen(file_name,'a+');
                
                max_temperature = T(pointer);
                min_temperature = T(pointer+system.chip.N);
                x_temp = x(left : right)-273;
                avg_temperature = sum(x_temp, 'all')/length(x_temp);
                temp_matrix((i-1)*3+1) = max_temperature;
                temp_matrix((i-1)*3+2) = min_temperature;
                temp_matrix(i*3) = avg_temperature;

            end
        end
        pointer = pointer+1;
    end
    
    %Chip area, pwoer density, tier substrate thickness, TSV diameter
    chip_area = chip(system.chip.N).Xsize*chip(system.chip.N).Ysize;
    chip_power_density = (chip(system.chip.N).power/chip_area)*1e-6;
    chip_matrix = [[chip_area  chip_power_density  chip(system.chip.N).layer.die(1)  chip(system.chip.N).TSV.d] temp_matrix];
    chip_cell = ['B', num2str(cell_count), ':', 'T', num2str(cell_count)];
    writematrix(chip_matrix,file_name_avg,'Sheet','temperatures_air_cooling','Range',chip_cell);
    writematrix(max_temp,'device_binary_rram_model/input_data/memory_max_temperatures.csv')
end

