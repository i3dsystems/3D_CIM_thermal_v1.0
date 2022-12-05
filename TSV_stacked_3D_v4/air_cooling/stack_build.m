function [Layer] = stack_build(system, chip)
%coding log: Dec 12 2015
%add TIM2 layer
    %for each die, it has three layers: die & ILD & Under
    % TIM + interposer;
    fprintf('Building stack layers\n');
    Layer.hsp = system.num.hsp + system.num.tim2;
    Layer.pack = system.num.pack;
    Layer.Nt = system.num.hsp + system.num.tim2 + system.num.pack; 
    
%     Layer_add.Nt = 0;
    %for interposer + hsp
    for i = 1 : 1 : system.layer.N %system tier or layer
        for j = 1 : 1 : system.chip.N %system chip
            if (chip(j).layer.id == i) %if first chip in that layer
                Layer.N(i) = chip(j).model;
                Layer.Nt = Layer.Nt + Layer.N(i);            
                break;
            end
            % die & ILD & RDL & under ||||| each die's TIM
        end
    end
    %%%add one layer for final ILD layer%%%%
%     Layer.Nt = Layer.Nt + 1;
%     %%%add one layer for package microbumps%%%%
%     Layer.Nt = Layer.Nt + 1;
    
    fprintf('there are %i physical layers\n', Layer.Nt);
    fprintf('there are %i nodal layers to be solved\n', Layer.Nt+1);
    fprintf('\n');
    Layer.material = zeros(Layer.Nt,1);
    Layer.thick = zeros(Layer.Nt,1);
    Layer.cal_type = zeros(Layer.Nt, 1);
    %there are three types; 
    %1: pure material, vertical = horizontal
    %2: vertically inserted material, vertical: parallel; horizontal: series
    %3: horizontally inserted material, vertical: series; horizontal: parralel
    
    Layer.other = zeros(Layer.Nt, 1); % -1 for TSV, -2 for microbumps; %other variable is used to specify what is the other type of material apart from the primary material
    Layer.portion = zeros(Layer.Nt, 1);
    %0 for no additional materials, corresponding to cal_type 1
    
    for i=1:1:Layer.hsp
        if i<=system.num.tim2
            Layer.material(i) = system.layer.tim2(2); %tim2
            Layer.thick(i) = system.layer.tim2(1)/system.num.tim2;
        else
            Layer.material(i) = system.layer.hsp(2); %heat spreader
            Layer.thick(i) = system.layer.hsp(1)/system.num.hsp;
        end
        Layer.cal_type(i) = 1;
    end
    
    for i = Layer.Nt : -1 : Layer.Nt - Layer.pack + 1
        Layer.material(i) = system.layer.pack(2); %package
        Layer.thick(i) = (system.layer.pack(1))/Layer.pack; %package;
        Layer.cal_type(i) = 3;
        Layer.other(i) =  4; %copper, needs to find an un-hardcoded way
        Layer.portion(i) = system.layer.metal_portion;
    end
    
    current_id = Layer.hsp+1;
    for i = 1 : 1 : system.layer.N
        %% 2. Die layers in every package
        for k = 1 : 1 : system.chip.N
            if (chip(k).layer.id == i) %%Check for chip 1, if chip 1, use the below steps, else a different step mentioned below
                for j=1:1:chip(k).model
                    if k ~= 1
                        j = j + 1;
                    end
                    if j==1
                        %% TIM1 layer between hsp and top chip first layer
                        if (current_id == 3)
                            Layer.material(current_id) = chip(k).layer.tim(2);
                            Layer.thick(current_id) = chip(k).layer.tim(1);
                            Layer.cal_type(current_id) = 1;
                            Layer.other(current_id) = 0;
                            current_id = current_id + 1;
                        end
                    elseif j == 2
                        %% die layer
                        Layer.material(current_id) = chip(k).layer.die(2);
                        Layer.thick(current_id) = chip(k).layer.die(1);
                        if i > 1
                            %top die is uniform; others contain TSV
                            Layer.cal_type(current_id) = 2;
                            Layer.other(current_id) = -1; 
                        else
                            Layer.cal_type(current_id) = 1;
                        end
                        current_id = current_id + 1;
                    elseif j == 3
                        %% metal layer (ILD)
                        Layer.material(current_id) = chip(k).layer.ild(2);
                        Layer.thick(current_id) = chip(k).layer.ild(1);
                        Layer.cal_type(current_id) = 3;
                        Layer.other(current_id) = 4; %copper
                        Layer.portion(current_id) = chip(k).metal_portion;
                        current_id = current_id + 1;
                    else
                        %% micro bump layer (for e.g. cu pillars)
                        Layer.material(current_id) = chip(k).layer.under(2);
                        Layer.thick(current_id) = chip(k).layer.under(1);
                        Layer.cal_type(current_id) = 2;
                        Layer.other(current_id) = -2;
                        current_id = current_id + 1;
                    end
                end
                break;
            end
        end
        
        
    end
end

