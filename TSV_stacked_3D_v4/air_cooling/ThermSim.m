function T = ThermSim(system, chip, Material, draw, h, cell_count, sim_case)
    start = tic;
    %Layer.material, Layer.thick, Layer.N
%     [Layer, Layer_add] = stack_build(system, chip);
    [Layer] = stack_build(system, chip);
    
    %mesh the chip
    [ chip, system] = mesh(chip, system);
%     [Ntime, Tmesh] = mesh_T(chip, system);
    
    %array initialization
    [Cdt_L, Cdt_V, Cap, Den, row, column, value, var, system] ...
    = initial(chip, system, Layer, Material);

    %non-uniform layer handling
%     [Cdt_L, Cdt_V, Cap, Den, Layer] = Material_fill(Cdt_L, Cdt_V, Cap, Den, system, ...
%                                   chip, Material, Layer, draw.C);    
    [Cdt_L, Cdt_V, Cap, Den] = Via_insetion(Cdt_L, Cdt_V, Cap, Den, system, ...
                                  chip, Material, Layer, draw.C);
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%steady state
    if system.tran ~= 1
        %#1 build the conductance
        [Y, ~, H] = Matrix_build_FOWLP(chip, system, h, Cdt_L, Cdt_V, Cap, Den, ...
                                     Layer, row, column, value, var);

        %#2 add power excitation
        heat = add_heat_FOWLP(chip, system, chip(1), var, draw.P, Layer, draw);
        
        %#3 steady state thermal analysis
        x = thermal_solve_ss(Y, H, heat, h.Ta, var);
        
        %#4 draw map        
        draw_map( chip, system, Layer, h, x, draw, 0);
        
        %#5 return the maximum temperature
        T = Tmax_get_FOWLP(x, system, chip, Layer, draw, 1, cell_count, sim_case);
    
    %%transient state
    else      
        %matrix build, including C matrix
        [Y, C, H] = Matrix_build_FOWLP(system, h, Cdt_L, Cdt_V, Cap, Den, ...
                                     Layer, row, column, value, var);
        %find initial trace
        chip1 = load('chip1.trace');
        chip(1).map(1:chip(1).blk_num(1),5) = chip1(3+1:3+chip1(2));
        chip(1).power(1) = chip1(3);
        
        trace_pointer = 1;                
        
        %add power excitation
        heat = add_heat_FOWLP(chip, system, system.pack, var, draw.P, Layer);
        
        %solver steady state problem firstly
        x = thermal_solve_ss(Y, H, heat, h.Ta, var);     
%         chip1 = load('chip1new.trace');
        %draw map
        draw_map( chip, system, Layer, h, x, draw, 0);
        
        %#5 return the maximum temperature
        fprintf('Time: %.1f ms\n', 0);        
        T = Tmax_get_FOWLP(x, system, chip, Layer, draw, 0);        
        [R, R_tran, S, S_tran] = thermal_fact(Y, C, H, system.dt);
        i = 2;   
        %transient iterations
        while i<= system.Ntime
%             if ( i == 20)
%                 fprintf('STARTING %.1f th SIMULATION \n', i);
%             end
            
            fprintf('Time: %.1f ms\n', Tmesh(i)*1e3);
            prev_heat = heat;
            prev_x = x;
            if Tmesh(i) >= chip1(trace_pointer*(chip1(2)+3)+1)
                const = trace_pointer*(chip1(2)+3);
                trace_pointer = trace_pointer+1;
                chip(1).map(1:chip(1).blk_num(1),5) = chip1(3+1+const:3+const+chip1(2));
                chip(1).power(1) = chip1(3+const);                                                                              
                fprintf('Time: %.1f ms: assigning power map\n', Tmesh(i)*1e3);
                heat = add_heat_FOWLP(chip, system, system.pack, var, draw.P, Layer);
            end
            
%             x = thermal_solve(Y, C, H, heat, prev_x, prev_heat, h.Ta, var, Tmesh(i)-Tmesh(i-1)); 
            x = thermal_solve_fact(R, R_tran, S, S_tran, Y, C, H, heat, prev_x, prev_heat, h.Ta, var, system.dt); 
            draw_map( chip, system, Layer, h, x, draw, Tmesh(i));
            
            %return the maximum temperature
            T = Tmax_get_FOWLP(x, system, chip, Layer, draw, Tmesh(i));
            
            i = i+1;
        end
    end
    fprintf('total run time: %.2f \n', toc(start));
end