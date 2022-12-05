function [Cl, Cv, Cp, De, Layer_out] = Material_fill(Cdt_L, Cdt_V, Cap, Den, system, ...
                                 chip, Material, Layer_in, draw)
% change the conductivity of meshes which contain TSVs
    fprintf('Molding material filling \n');
    if draw == 1
        drawC = NaN * ones(system.pack.Ny, system.pack.Nx, 2*system.layer.N);
    end
    layer_id = 1;
    
    Layer_in.mask = zeros(system.pack.Ny, system.pack.Nx,system.chip.N);
    offset = system.Nhsp;
%     for layer = 1 : system.layer.N % loop till all layers
%         for i = 1 : layer % loop till current layer
%             if (layer == 1)
%                 offset = system.Nhsp + (Layer_in.N(layer) - 1)*system.pack.Nx * system.pack.Ny;
%             else
% %                 for i = 1 : 1 : system.chip.N % loop through all chips
% %                     if (chip(i).layer.id == j) % if current chip is first chip, use its modeled number of layers to calculate offset
%                 offset = (Layer_in.N(layer) - 1)*system.pack.Nx * system.pack.Ny + offset;
%                         %             offset = Layer_in.hsp * system.hsp.Nx * system.hsp.Ny + Layer_in.N(j-1)*system.pack.Nx * system.pack.Ny + offset;
%             end

    for layer = 1 : system.layer.N % loop till all layers
        switch layer
            case {1}
                
            case {2}
                offset = offset + (Layer_in.N(layer-1)-1)*system.pack.Nx * system.pack.Ny;
            otherwise
                offset = offset + (Layer_in.N(layer-1))*system.pack.Nx * system.pack.Ny;
        end
        
        mask = zeros(system.pack.Ny, system.pack.Nx);
        for i = 1:1:system.chip.N
            temp_mask = zeros(system.pack.Ny, system.pack.Nx);
            if (chip(i).layer.id == layer)
                temp_mask = mask;
                array = ChipInPlane(chip(i), system.pack);
                xl = array(1); yb = array(3);
                xr = array(2); yt = array(4);                
                mask(yb:yt, xl:xr) = 1;
                Layer_in.mask(:,:,i) = mask - temp_mask;
            
            else
                continue;
            end
            
        end
        
%         Layer_in(layer).mask = mask;
        
        for k = 1 : 1 : system.chip.N
            if (chip(k).layer.id == layer-1)
                for j = 1 : 1 : chip(k).model
                    const = system.pack.Nx * system.pack.Ny;
                    drawC(:,:,layer_id) = reshape(Cdt_V(offset+1+const*(j-1):offset+j*const), system.pack.Nx, system.pack.Ny)';
                    for ii = 1 : system.pack.Nx
                        for jj = 1 : system.pack.Ny
                            id = offset + (j-1)* system.pack.Nx * system.pack.Ny + ii + (jj-1)*system.pack.Nx;
                            if (mask(jj, ii)==1)
                                continue;
                            end
                            %cv = Material.K(chip(k).layer.ild(2));
                            %cl = Material.K(chip(k).layer.ild(2));
                            %cap = Material.C(chip(k).layer.ild(2));
                            %den = Material.D(chip(k).layer.ild(2));
                            cv = 2.7505;
                            cl = 200.69;
                            cap = 545;
                            den = 5669;
                            %end
                            Cdt_V(id) = cv;
                            Cdt_L(id, 1:2) = cl;
                            Cap(id) = cap;
                            Den(id) = den;
                            if draw == 1
                                drawC(jj, ii, layer_id) = Cdt_L(id, 1);
                            end
                        end
                    end
                    
                    
                    
                    %                     cv = Material.K(chip(i).layer.mold_b(2));
                    %                     cl = Material.K(chip(i).layer.mold_b(2));
                    %                     cap = Material.C(chip(i).layer.mold_b(2));
                    %                     den = Material.D(chip(i).layer.mold_b(2));
                    if draw == 1
                        figure(layer_id)
                        pcolor(system.pack.Xmesh*100, system.pack.Ymesh*100, drawC(:,:,layer_id));
                        h=colorbar;
                        set(get(h,'Title'),'string','k*m/W','FontSize',16)
                        set(gca,'FontSize',16);
                        xlabel('x(cm)');
                        ylabel('y(cm)');set(gca,'FontSize',16);
                        if j == 1
                            str = ['Layer ', num2str(layer), ', Die'];
                            title(str);
                        else
                            str = ['Layer ', num2str(layer), ', Metal'];
                            title(str);
                        end
                    end
                    layer_id = layer_id + 1;
                end
            else
                continue;
            end
        end
    end
    

    fprintf('\n');
    Cl = Cdt_L;
    Cv = Cdt_V;
    Cp = Cap;
    De = Den;
    Layer_out = Layer_in;
end

