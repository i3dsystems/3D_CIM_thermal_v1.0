function [Cdt_L, Cdt_V, Cap, Den, row, column, value, var, system_] ...
          = initial(chip, system, Layer, Material)
%initialization
    start = 0;
    for l= 1 : 1 : system.layer.N %%To find out the total number of dies in the system (just skip to line 12)
        for i = 1:system.chip.N
            if (chip(i).layer.id == l)
                start = start + 1; %chip(l).N;
            end
        end
    end
    system.Ndie = start;
    var = 0;
    var = var + system.hsp.Nx * system.hsp.Ny * (Layer.hsp+1) + system.pack.Nx * system.pack.Ny * (Layer.pack + 1); % In 'Layer.hsp+1', the +1 is for TIM2 layer between heat spreader and first silicon layer.
    for i = 1 : 1 : system.layer.N
        if i==1
            var = chip(i).Nx * chip(i).Ny * (Layer.N(i)-1) + var;
        else
            var = chip(i).Nx * chip(i).Ny * (Layer.N(i)) + var;
        end
    end
%     var = var + (system.pack.Nx * system.pack.Ny * (1)); % The '1' here is for the final ILD layer after the embedded tiers.

    row = zeros(var*7,1);
    column = zeros(var*7,1);
    value = zeros(var*7,3);
    
    var_hsp = system.hsp.Nx * system.hsp.Ny;
    var_chip = chip(system.chip.N).Nx * chip(system.chip.N).Ny;
    var_pack = system.pack.Nx * system.pack.Ny;
    
%conductivity arrays; the number of conductivity arrays elements are less than the
%number of var elements.
    Cdt_L = zeros(var - var_chip - var_pack, 2);
    Cdt_V = zeros(var - var_chip - var_pack, 1);
    Cap = zeros(var - var_chip - var_pack, 1);
    Den = zeros(var - var_chip - var_pack, 1);             
    Base_material = Layer.material;
    
    right_id = 0;
    for layer_id = 1:1:Layer.hsp+1
        left_id = right_id + 1;
        right_id = system.hsp.Nx*system.hsp.Ny+left_id - 1;
        if Layer.cal_type(layer_id) == 3
            k1 = Material.K(Base_material(layer_id) );
            k2 = Material.K(Layer.other(layer_id) );
            cp1 = Material.C(Base_material(layer_id) );
            cp2 = Material.C(Layer.other(layer_id) );
            den1 = Material.D(Base_material(layer_id) ); 
            den2 = Material.D(Layer.other(layer_id) );
            
            k = [k1, k2];
            portion = Layer.portion(layer_id);
            p = [1-portion, portion];
            [cv, cl] = Cal_Con(k, p);
            
            cap = cp1*portion + cp2*(1-portion);
            den = den1*portion + den2*(1-portion);
        else
            cv = Material.K(Base_material(layer_id) );
            cl = Material.K(Base_material(layer_id) );
            cap = Material.C(Base_material(layer_id) );
            den = Material.D(Base_material(layer_id) );
        end
            
        Cdt_L(left_id:right_id, 1:2) = cl;
        Cdt_V(left_id:right_id, 1) = cv;
        Cap(left_id:right_id, 1) = cap;
        Den(left_id:right_id, 1) = den;
    end    
    system.Nhsp = right_id;
%     layer_id = layer_id + 1;

    %initial the conductivity arrays in chip domain
    for i = 1 : 1 : system.layer.N
        for j = 1 : 1 : Layer.N(i)
            %the TIM is pre-defined with a big layer, and is skipped inside
            %this loop
            if j == 1 && i == 1
                layer_id = layer_id + 1;
                continue;
            end
%             layer_id = layer_id + 1;
            left_id = right_id + 1;
            right_id = chip(i).Nx*chip(i).Ny + left_id - 1;
            if Layer.cal_type(layer_id) == 3
                k1 = Material.K(Base_material(layer_id) );
                k2 = Material.K(Layer.other(layer_id) );
                cp1 = Material.C(Base_material(layer_id) );
                cp2 = Material.C(Layer.other(layer_id) );
                den1 = Material.D(Base_material(layer_id) );
                den2 = Material.D(Layer.other(layer_id) );
                
                k = [k1, k2];
                portion = Layer.portion(layer_id);
                p = [1-portion, portion];
                [cv, cl] = Cal_Con(k, p);
                
                cap = cp1*portion + cp2*(1-portion);
                den = den1*portion + den2*(1-portion);
            else
                cv = Material.K(Base_material(layer_id) );
                cl = Material.K(Base_material(layer_id) );
                cap = Material.C(Base_material(layer_id) );
                den = Material.D(Base_material(layer_id) );
            end
            
%             if (left_id == 529)
%                 a = 1;
%             end
            
            Cdt_L(left_id:right_id, 1:2) = cl;
            Cdt_V(left_id:right_id, 1) = cv;
            Cap(left_id:right_id, 1) = cap;
            Den(left_id:right_id, 1) = den;
            
            layer_id = layer_id + 1;
            
        end
        
    end
    system.Nchip = right_id;
%initial the conductivity arrays in package arrays
    for layer_id = Layer.Nt - Layer.pack + 1 : 1 : Layer.Nt
        left_id = right_id + 1;
        right_id = system.pack.Nx*system.pack.Ny+left_id - 1;
        if Layer.cal_type(layer_id) == 3
            k1 = Material.K(Base_material(layer_id) );
            k2 = Material.K(Layer.other(layer_id) );
            cp1 = Material.C(Base_material(layer_id) );
            cp2 = Material.C(Layer.other(layer_id) );
            den1 = Material.D(Base_material(layer_id) ); 
            den2 = Material.D(Layer.other(layer_id) );
            
            k = [k1, k2];
            portion = Layer.portion(layer_id);
            p = [1-portion, portion];
            [cv, cl] = Cal_Con(k, p);
            
            cap = cp1*portion + cp2*(1-portion);
            den = den1*portion + den2*(1-portion);
        else
            cv = Material.K(Base_material(layer_id) );
            cl = Material.K(Base_material(layer_id) );
            cap = Material.C(Base_material(layer_id) );
            den = Material.D(Base_material(layer_id) );
        end
            
        Cdt_L(left_id:right_id, 1:2) = cl;
        Cdt_V(left_id:right_id, 1) = cv;
        Cap(left_id:right_id, 1) = cap;
        Den(left_id:right_id, 1) = den;     
    end
    
                
    system_ = system;
end


