function [Cl, Cv, Cp, De] = Via_insetion(Cdt_L, Cdt_V, Cap, Den, system, ...
                                 chip, Material, Layer, draw)
% change the conductivity of meshes which contain TSVs
    fprintf('Vias insertion \n');
    layer_id = Layer.hsp + 1;
    right_id = system.Nhsp;
%     chipid = 0;
    for i = 1:1:system.layer.N        
        %% via conductivity modeling is based on the following paper:
        %doi: 10.1109/ASPDAC.2013.6509643
        if draw == 1
            drawC_TSV = NaN * ones(system.pack.Ny, system.pack.Nx, 2);
            drawC_Bump = NaN * ones(system.pack.Ny, system.pack.Nx, 2);
        end        
        for j = 1 : system.chip.N
            if (chip(j).layer.id == i)
                break;
            end
        end
%         if(chipid <= system.chip.N)
%             chipid = chipid + 1;
%         end
%             
        chipid = j;
        
        array = ChipInPlane(chip(chipid), system.pack);
        xl = array(1); yb = array(3);
        
        %% bump
        bump.kl = Material.K( chip(chipid).bump.material );
        bump.kv = Material.K( chip(chipid).bump.material );
        bump.cap = Material.C( chip(chipid).bump.material );
        bump.den = Material.D( chip(chipid).bump.material );        
        fprintf('layer %d bump informtion: ', i);
        bump_map = zeros(chip(chipid).Ny, chip(chipid).Nx, 2);
        TSV_map = zeros(chip(chipid).Ny, chip(chipid).Nx, 2);
        for jj = 1:size(chip(chipid).bump.map)
            chip(chipid).bump.xl = chip(chipid).bump.map(jj,1);
            chip(chipid).bump.yb = chip(chipid).bump.map(jj,2);
            chip(chipid).bump.xr = chip(chipid).bump.xl + chip(chipid).bump.map(jj,3);
            chip(chipid).bump.yt = chip(chipid).bump.yb + chip(chipid).bump.map(jj,4);
%             tmp_map = via_count_assign(system.pack, chip(chipid).bump);
            tmp_map = via_count_assign(chip(chipid), chip(chipid).bump);
            bump_map(tmp_map ~= 0) = tmp_map(tmp_map ~= 0);
        end      
        if (i > 1)
            liner = Material.K( chip(chipid).TSV.material(2) );
            core = Material.K( chip(chipid).TSV.material(1) );        
            dm = chip(chipid).TSV.d - 2*chip(chipid).TSV.liner;
            dl = chip(chipid).TSV.liner;
            r = core/liner;

            tsv.kl = (dm^2*(r-1)/((dm+2*dl)*(dm+2*dl*r))+1)*liner;
            tsv.kv = (dm^2*(r-1)/(dm+2*dl)^2+1)*liner;     
            tsv.cap = Material.C( chip(chipid).TSV.material(1) );
            tsv.den = Material.D( chip(chipid).TSV.material(1) );
            fprintf('layer% i TSV informtion: ', i);
            for jj = 1:size(chip(chipid).TSV.map)
                chip(chipid).TSV.xl = chip(chipid).TSV.map(jj,1);
                chip(chipid).TSV.yb = chip(chipid).TSV.map(jj,2);
                chip(chipid).TSV.xr = chip(chipid).TSV.xl + chip(chipid).TSV.map(jj,3);
                chip(chipid).TSV.yt = chip(chipid).TSV.yb + chip(chipid).TSV.map(jj,4);
%                 tmp_map = via_count_assign(system.pack, chip(chipid).TSV);
                tmp_map = via_count_assign(chip(chipid), chip(chipid).TSV);
                TSV_map(tmp_map ~= 0) = tmp_map(tmp_map ~= 0);
            end
        end
                
        for j = 1 : 1 : Layer.N(i)
            %% the tim layer is pre-defined with a big layer
            if j == 1 && i == 1
                layer_id = layer_id + 1;
                continue;
            end
%             if ( layer_id == 7)
%                 a = 1;
%             end
            %% TSV assignment
            if ((Layer.cal_type(layer_id) == 2 && Layer.other(layer_id) == -1) ...
                    || (i > 1 && j == 2 ) ) %second term is for metal and RDL
%                 for ii = 1 : system.pack.Nx-1
                for ii = 1 : chip(chipid).Nx-1
%                     for jj = 1 : system.pack.Ny-1
                    for jj = 1:chip(chipid).Ny-1
%                         id = right_id + ii + (jj-1) * system.pack.Nx;
                        id = right_id + ii + (jj-1) * chip(chipid).Nx;
                        if ( id == 5634364)
                            a = 1;
                        end
                        kcell_L = Cdt_L(id,1);
                        kcell_V = Cdt_V(id);
                        
                        %Ankit: Change in grid conductivity due to TSV
                        %insertion needs to account for both vertical and
                        %lateral conductivities. Thus, in the formulas
                        %below we account for both Cdt_L and Cdt_V by using
                        %kcell_L and kcell_V. The paper referred (above)
                        %for via conductivity modeling does not define 
                        %whether the variable 'kcell' is lateral or vertical
                        %conductivity. Ideally, the value of kcell should
                        %be the conductivity in the direction of heat
                        %conduction considered (x, y, or z, i.e. lateral or
                        %vertical).
                        
%                         W = system.pack.Xmesh( ii + 1 ) - system.pack.Xmesh( ii );
%                         H = system.pack.Ymesh( jj + 1 ) - system.pack.Ymesh( jj );
                        W = chip(chipid).Xmesh( ii + 1 ) - chip(chipid).Xmesh( ii );
                        H = chip(chipid).Ymesh( jj + 1 ) - chip(chipid).Ymesh( jj );
                        h = chip(chipid).TSV.d;
                        w = chip(chipid).TSV.d;
                        betaL = tsv.kl/kcell_L;
                        betaV = tsv.kv/kcell_V;  
                        Nx = TSV_map(jj, ii, 1);
                        Ny = TSV_map(jj, ii, 2);
                        Cdt_L(id, 1) = (h*w*(betaL-1)*Nx*Ny/(W*H*(1+(1-Nx*w/W)*(betaL-1)))+1)*kcell_L;
                        Cdt_L(id, 2) = (h*w*(betaL-1)*Nx*Ny/(W*H*(1+(1-Ny*h/H)*(betaL-1)))+1)*kcell_L;
                        Cdt_V(id) = (h*w*(betaV-1)*Nx*Ny/(W*H)+1)*kcell_V;
                        Cap(id) = (h*w*(tsv.cap/Cap(id)-1)*Nx*Ny/(W*H)+1)*Cap(id);
                        Den(id) = (h*w*(tsv.den/Den(id)-1)*Nx*Ny/(W*H)+1)*Den(id);
                        if draw == 1
                            drawC_TSV(yb+jj-1, xl+ii-1, 1) = Cdt_L(id, 1);
                            drawC_TSV(yb+jj-1, xl+ii-1, 2) = Cdt_V(id);
                        end
                    end
                end
            end
            
            %% ubump assignment
            if Layer.cal_type(layer_id) == 2 && Layer.other(layer_id) == -2
                for ii = 1:chip(i).Nx-1
                    for jj = 1:chip(i).Ny-1
                        id = right_id + ii + (jj-1) * chip(chipid).Nx;
                        if ( id == 5634364)
                            a = 1;
                        end
                        kcell_L = Cdt_L(id,1);
                        kcell_V = Cdt_V(id);
                        kcell = Cdt_V(id);
                        W = chip(chipid).Xmesh( ii + 1 ) - chip(chipid).Xmesh( ii );
                        H = chip(chipid).Ymesh( jj + 1 ) - chip(chipid).Ymesh( jj );
                        Nx = bump_map(jj, ii, 1);
                        Ny = bump_map(jj, ii, 2);
                        h = chip(chipid).bump.d;
                        w = chip(chipid).bump.d;
                        betaL = bump.kl/kcell;
                        betaV = bump.kv/kcell;    
                        viaCap = bump.cap;
                        viaDen = bump.den;                                               

                        Cdt_L(id, 1) = (h*w*(betaL-1)*Nx*Ny/(W*H*(1+(1-Nx*w/W)*(betaL-1)))+1)*kcell;
                        Cdt_L(id, 2) = (h*w*(betaL-1)*Nx*Ny/(W*H*(1+(1-Ny*h/H)*(betaL-1)))+1)*kcell;
                        Cdt_V(id) = (h*w*(betaV-1)*Nx*Ny/(W*H)+1)*kcell;
                        Cap(id) = (h*w*(viaCap/Cap(id)-1)*Nx*Ny/(W*H)+1)*Cap(id);
                        Den(id) = (h*w*(viaDen/Den(id)-1)*Nx*Ny/(W*H)+1)*Den(id);
                        if draw == 1
                            drawC_Bump(yb+jj-1, xl+ii-1, 1) = Cdt_L(id, 1);
                            drawC_Bump(yb+jj-1, xl+ii-1, 2) = Cdt_V(id);
                        end                        
                    end
                end
            end            
%             right_id = system.pack.Nx*system.pack.Ny + right_id;
            right_id = chip(chipid).Nx*chip(chipid).Ny + right_id;
            layer_id = layer_id + 1;
        end          
        if draw == 1
            if (i>1)
                figure()
                pcolor(system.pack.Xmesh*100, system.pack.Ymesh*100, drawC_TSV(:,:,1));
                h=colorbar;
                set(get(h,'Title'),'string','k*m/W','FontSize',16)
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
                str = ['layer ', num2str(i), ' TSV layer X/Y conductivity'];
                title(str);

                figure()
                pcolor(system.pack.Xmesh*100, system.pack.Ymesh*100, drawC_TSV(:,:,2));
                h=colorbar;
                set(get(h,'Title'),'string','k*m/W','FontSize',16)
                set(gca,'FontSize',16);
                xlabel('x(cm)');
                ylabel('y(cm)');set(gca,'FontSize',16);
                str = ['layer ', num2str(i), ' TSV layer Z conductivity'];
                title(str);
            end

            figure()
            pcolor(system.pack.Xmesh*100, system.pack.Ymesh*100, drawC_Bump(:,:,1));
            h=colorbar;
            set(get(h,'Title'),'string','k*m/W','FontSize',16)
            set(gca,'FontSize',16);
            xlabel('x(cm)');
            ylabel('y(cm)');set(gca,'FontSize',16);
            str = ['layer ', num2str(i), ' Bump layer X/Y conductivity'];
            title(str);            

            figure()
            pcolor(system.pack.Xmesh*100, system.pack.Ymesh*100, drawC_Bump(:,:,2));
            h=colorbar;
            set(get(h,'Title'),'string','k*m/W','FontSize',16)
            set(gca,'FontSize',16);
            xlabel('x(cm)');
            ylabel('y(cm)');set(gca,'FontSize',16);
            str = ['layer ', num2str(i), ' Bump layer Z conductivity'];
            title(str);            
        end        
    end
        
    fprintf('\n');
    Cl = Cdt_L;
    Cv = Cdt_V;
    Cp = Cap;
    De = Den;
end

