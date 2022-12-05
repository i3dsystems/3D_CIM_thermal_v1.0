function A = mem_power_map_gen(x_margin_MEM, y_margin_MEM, nMEM, nyMEM, nxMEM, block_dimension, y_pitch, x_pitch, per_MEM_power)

x_ctr = 0;
y_ctr = 1;
i = 1;
while (i<=nMEM)
    while(x_ctr < nxMEM)
        while(y_ctr <= nyMEM)
            
            if(x_ctr == 0)
                if(y_ctr == 1)
                    A = [x_margin_MEM y_margin_MEM block_dimension block_dimension per_MEM_power 0];
                    i = i+2;
                else
                    B = [x_margin_MEM (A((x_ctr*nyMEM*2 + y_ctr)-1, 2)+((mem_y_dimension+y_pitch))) block_dimension block_dimension per_MEM_power 0];
                    A = [A; B];
                    i = i+2;
                end
            else
                if(y_ctr == 1)
                    B = [(A((x_ctr*nyMEM*2 + 1)-1, 1)+((block_dimension+x_pitch))) y_margin_MEM block_dimension block_dimension per_MEM_power 0];
                    A = [A; B];
                    i = i+2;
                else
                    B = [(A((x_ctr*nyMEM*2 + 1)-1, 1)+((block_dimension+x_pitch))) (A((x_ctr*nyMEM*2 + y_ctr)-1, 2)+((mem_y_dimension+y_pitch))) block_dimension block_dimension per_MEM_power 0];
                    A = [A; B];
                    i = i+2;
                end
            end
            y_ctr = y_ctr + 2;
        end
        x_ctr = x_ctr + 1;
        y_ctr = 1;
    end
end
% whos A;
end