function A = logic_power_map_gen(x_margin_adc, y_margin_adc, nyADC, nxADC, block_dimension, y_pitch, x_pitch, per_ADC_power)

x_ctr = 0;
y_ctr = 1;
% i = 1;
% while (i<=nADC)
    while(x_ctr < nxADC)
        while(y_ctr <= nyADC)
            
            if(x_ctr == 0)
                if(y_ctr == 1)
                    B = [x_margin_adc y_margin_adc block_dimension block_dimension per_ADC_power 0];
                    A = B;
%                     i = i+1;
                else
                    B = [x_margin_adc (A((x_ctr*nyADC + y_ctr)-1, 2)+((block_dimension+y_pitch))) block_dimension block_dimension per_ADC_power 0];
                    A = [A; B];
%                     i = i+1;
                end
            else
                if(y_ctr == 1)
                    B = [(A((x_ctr*nyADC + 1)-1, 1)+((block_dimension+x_pitch))) y_margin_adc block_dimension block_dimension per_ADC_power 0];
                    A = [A; B];
%                     i = i+1;
                else
                    B = [(A((x_ctr*nyADC + 1)-1, 1)+((block_dimension+x_pitch))) (A((x_ctr*nyADC + y_ctr)-1, 2)+((block_dimension+y_pitch))) block_dimension block_dimension per_ADC_power 0];
                    A = [A; B];
%                     i = i+1;
                end
            end
            y_ctr = y_ctr + 1;
        end
        x_ctr = x_ctr + 1;
        y_ctr = 1;
    end
% end
% whos A;
end