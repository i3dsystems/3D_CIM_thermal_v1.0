function [chip_, system_] = mesh_chip( chip, system )
%mesh the chip area
%the basic idea is mesh the system.box.x between adjacent points
    %% mesh the X arrays
    Xarray = [];
    Yarray = [];
    lenX = length(system.box.x);
    lenY = length(system.box.y);
    for i = 1 : 1 : lenX - 1
        if system.box.x(i) + system.chip.Xgrid > system.box.x(i+1)
            meshX = [system.box.x(i), system.box.x(i+1)];
        else
            meshX = [system.box.x(i) : system.chip.Xgrid : system.box.x(i+1), system.box.x(i+1)];
        end
        Xarray = [Xarray, meshX];
    end
    Xarray = my_unique(Xarray, 1e-10);
    system.box.Xmesh = Xarray; % 'system.box.Xmesh' contains the mesh element center positions
    system.box.Nx = length(Xarray);
    
    for i = 1 : 1 : lenY-1
        if system.box.y(i) + system.chip.Ygrid > system.box.y(i+1)
            meshY = [system.box.y(i), system.box.y(i+1)];
        else
            meshY = [system.box.y(i) : system.chip.Ygrid : system.box.y(i+1), system.box.y(i+1)];
        end
        Yarray = [Yarray, meshY];                
    end
    
    Yarray = my_unique(Yarray, 1e-10);
    system.box.Ymesh = Yarray; % 'system.box.Ymesh' contains the mesh element center positions
    system.box.Ny = length(Yarray);    
    
    for i = 1 : 1 : system.chip.N
        Nxl = find(abs(Xarray - chip(i).xl)<1e-12);
        Nxr = find(abs(Xarray - chip(i).xr)<1e-12);
        chip(i).Xmesh = Xarray(Nxl:Nxr);
        chip(i).Nx = Nxr - Nxl + 1;
        Nyb = find(abs(Yarray - chip(i).yb)<1e-12);
        Nyt = find(abs(Yarray - chip(i).yt)<1e-12);
        chip(i).Ymesh = Yarray(Nyb:Nyt);
        chip(i).Ny = Nyt - Nyb + 1;
    end
    
    chip_ = chip;
    system_ = system;
    fprintf('chip domain meshing done\n')
end

