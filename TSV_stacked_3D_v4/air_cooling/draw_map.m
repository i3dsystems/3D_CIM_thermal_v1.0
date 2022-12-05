function draw_map( chip, system, Layer, h, x, draw, t)
    %%%%%%%%%%variable conversion%%%%%%%%%%%%%%%%%%
    if draw.T == 0
        return
    end
    Ta = h.Ta;
    hs = h.up;
%     drawT_pack= zeros(chip(1).Nx, chip(1).Ny, system.layer.N);
%     const = chip(1).Nx*chip(1).Ny;    
%     for m = 1:1:system.chip.N        
%         index_start = system.Nhsp + chip(m).Nx*chip(m).Ny * sum(Layer.N(1:chip(m).layer.id-1));    
%  
%         Chip_xl = find(abs(chip(m).Xmesh-chip(m).xl)<1e-12);
%         Chip_yb = find(abs(chip(m).Ymesh-chip(m).yb)<1e-12);
%         Chip_xr = find(abs(chip(m).Xmesh-chip(m).xr)<1e-12);
%         Chip_yt = find(abs(chip(m).Ymesh-chip(m).yt)<1e-12);
%         
%         id = chip(m).Nx*chip(m).Ny + index_start;
%         drawT_pack(:,:, chip(m).layer.id) = reshape(x(id+1:id+const), chip(m).Nx, chip(m).Ny)'; 
%         drawT_die = drawT_pack(Chip_yb:Chip_yt, Chip_xl:Chip_xr, chip(m).layer.id);
% %         die_temp_copy = drawT_die-273;
% %         die_temp_copy = reshape((drawT_die-273), (Chip_xr - Chip_xl + 1)*(Chip_yt-Chip_yb+1), 1);


    drawT_pack= NaN*ones(system.pack.Ny, system.pack.Nx);
    index = 20;
    for m = 1:1:system.chip.N
        drawT_die = zeros(chip(m).Ny, chip(m).Nx);
        drawT_metal = zeros(chip(m).Ny, chip(m).Nx);
        const = chip(m).Nx*chip(m).Ny;
        
        index_start = system.Nhsp;
        for l= 1 : 1 : m-1
            index_start = index_start + chip(l).Nx*chip(l).Ny.*(Layer.N(l)-1);
        end
        
%         for k=1:1:chip(m).N
            id = const+index_start;
            drawT_die(:,:) = reshape(x(id+1:id+const), chip(m).Nx, chip(m).Ny)';
%             drawT_die = drawT_pack(Chip_yb:Chip_yt, Chip_xl:Chip_xr, chip(m).layer.id);
            drawT_metal(:,:) = reshape(x(id+const+1:id+2*const), chip(m).Nx, chip(m).Ny)';
%         end


        %% die plotting            
        map = chip(m).map;
        if size(chip(m).map,1)>0
            blk_name = chip(m).blk_name;
        else
            blk_name = chip(m).name;
        end
        index = index + 1;   
        if draw.gif ~= 1
            [Tmax, Tmin] = DrawSteady(index, chip(m).Xmesh, chip(m).Ymesh, drawT_die, ...
                                      draw, Ta, map, blk_name);
            string = ['chip', num2str(m), ' Map'];
            title(string);
        else
            name = ['chip', num2str(m), '.gif'];                
            [Tmax, Tmin] = DrawSteady(index, chip(m).Xmesh, chip(m).Ymesh, drawT_die, ...
                                      draw, Ta, map, name, t);                
        end
        fprintf('chip%i Max: %f\n', m, Tmax);
        fprintf('chip%i Min: %f\n', m, Tmin);
        
        %% metal plotting
        index = index + 1;
        if draw.gif ~= 1
            [Tmax, Tmin] = DrawSteady(index, chip(m).Xmesh, chip(m).Ymesh, drawT_metal, ...
                                      draw, Ta, map, blk_name);
            string = ['chip', num2str(m), ' Metal Map'];
            title(string);
        else
            name = ['chip', num2str(m), '_Metal.gif'];         
            [Tmax, Tmin] = DrawSteady(index, chip(m).Xmesh, chip(m).Ymesh, drawT_metal(:,:,k), ...
                                      draw, Ta, map, name, t);
        end

        fprintf('chip%i Metal Max: %f\n', m, Tmax);
        fprintf('chip%i Metal Min: %f\n', m, Tmin);            

%         if k == chip(m).N
%             Chip_xl = find(abs(system.pack.Xmesh-chip(m).xl)<1e-12);
%             Chip_yb = find(abs(system.pack.Ymesh-chip(m).yb)<1e-12);
%             Chip_xr = find(abs(system.pack.Xmesh-chip(m).xr)<1e-12);
%             Chip_yt = find(abs(system.pack.Ymesh-chip(m).yt)<1e-12);    
%             drawT_pack(Chip_yb:Chip_yt,Chip_xl:Chip_xr) = drawT_die(:,:,k);
%         end
    end

    %% plot the package figures
%     for j = 1:system.layer.N
%         l = 0;
%         tmp = chip;
%         for ii = 1:system.chip.N
%             if chip(ii).layer.id == j
%                 l = l + 1;                    
%                 tmp(l) = chip(ii);
%             end
%         end
%         [map, chip_name] = chip_array(tmp, l);
%         index = 20+system.Ndie+j;
%         if draw.gif ~= 1
%             [Tmax, Tmin] = DrawSteady(index, system.pack.Xmesh, system.pack.Ymesh, drawT_pack(:,:,j), ...
%                                       draw, Ta, map, chip_name);
%             %title('Package');
%         else
%             name = 'package.gif';
%             [Tmax, Tmin] = DrawSteady(index, system.pack.Xmesh, system.pack.Ymesh, drawT_pack(:,:,j), ...
%                                       draw, Ta, map, name, t); 
%         end    
%         fprintf('Active Layer #%d Max: %f\n', j, Tmax);
%         fprintf('Active Layer #%d Min: %f\n', j, Tmin);
%     end

    [map, chip_name] = chip_array(chip, system.chip.N);
    const = system.pack.Nx * system.pack.Ny;
    drawT_pack = reshape(x(system.Nchip+1:system.Nchip+const), system.pack.Nx, system.pack.Ny)';
    index = index + 1;
    if draw.gif ~= 1
        [Tmax, Tmin] = DrawSteady(index, system.pack.Xmesh, system.pack.Ymesh, drawT_pack, ...
                                  draw, Ta, map, chip_name);
        title('Package');
    else
        name = 'package.gif';
        [Tmax, Tmin] = DrawSteady(index, system.pack.Xmesh, system.pack.Ymesh, drawT_pack, ...
                                  draw, Ta, map, name, t); 
    end    
    fprintf('Package Max: %f\n', Tmax);
    fprintf('Package Min: %f\n', Tmin);
        
    %% plot the heat spreader
    const = system.hsp.Nx*system.hsp.Ny;%system.hsp.Ny * system.hsp.Nx;
    drawT_pack = reshape(x(1:const), system.hsp.Nx, system.hsp.Ny)';
    escape_power(system.hsp, drawT_pack, Ta, hs)
    
    index = index + 1;
    if draw.gif ~= 1
        [Tmax, Tmin] = DrawSteady(index, system.hsp.Xmesh, system.hsp.Ymesh, drawT_pack, ...
                                  draw, Ta, map, chip_name);
        title('Heat Spreader');
    else
        name = 'hsp.gif';
        [Tmax, Tmin] = DrawSteady(index, system.hsp.Xmesh, system.hsp.Ymesh, drawT_pack, ...
                                  draw, Ta, map, name, t); 
    end
    fprintf('Heat spreader Max: %f\n', Tmax);
    fprintf('Heat spreader Min: %f\n', Tmin);    
end
