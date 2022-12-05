function [Nx, Xmesh, Ny, Ymesh] = mesh_box(inner, outer)
%mesh the heat spreader using the results of system.pack.Xmesh/Ymesh
%useful references: system.pack.Xmesh, system.pack.Ymesh
    mesh_size_min = inner.Xmesh(2) - inner.Xmesh(1);
    %% mesh the left part of the system.pack    
    Num= floor( ( inner.xl - outer.xl )/mesh_size_min ) + 1;
    %%%this is for ticking out the unused element in the array
    mesh = -ones(1, Num);
    
    mesh(1) = inner.xl;
    mesh_size = mesh_size_min;
    i = 1;
    while( mesh_size < outer.Xgrid )
        if mesh(i) - mesh_size <= outer.xl %If the distance between leftmost inner element and leftmost outer element is smaller than the size of outer grid, the outer xl is used as next mesh element location.
            mesh(i+1) = outer.xl;
            i = i+1;
            break;
        end
        mesh(i+1) = mesh(i) - mesh_size;
        mesh_size = mesh_size + mesh_size_min;
        i = i + 1;
    end
    Num = floor( (mesh(i) - outer.xl)/outer.Xgrid ) + 1;
    grid = (mesh(i) - outer.xl)/Num;
    rest = outer.xl : grid : mesh(i); % Rest of mesh element positions are built starting from package (or hsp) xl to last element of chip (or package) (going in -ve x direction from center)
    
    mesh = [mesh, rest];
    id = mesh >= 0;
    mesh = mesh(id);
    
    mesh = sort(mesh);
    mesh1 = my_unique(mesh, 1e-10);
    
    %% mesh the right part of the system.pack
    Num= floor(( outer.xr - inner.xr )/mesh_size_min) + 1;
    %%%this is for ticking out the unused element in the array
    mesh = -ones(1, Num);
    
    mesh(1) = inner.xr;
    mesh_size = mesh_size_min;
    i = 1;
    while( mesh_size < outer.Xgrid )
        if mesh(i) + mesh_size >= outer.xr
            mesh(i+1) = outer.xr;
            i = i+1;
            break;
        end
        mesh(i+1) = mesh(i) + mesh_size;
        mesh_size = mesh_size + mesh_size_min;
        i = i + 1;
    end
    Num = floor( (outer.xr - mesh(i))/outer.Xgrid ) + 1;
    grid = (outer.xr - mesh(i))/Num;
    rest = [mesh(i) : grid : outer.xr, outer.xr];
    
    mesh = [mesh, rest];
    id = mesh >= 0;
    mesh = mesh(id);
    
    mesh = sort(mesh);
    mesh2 = my_unique(mesh, 1e-10);
    
    %% merge the three mesh parts as Xmesh
    Xmesh = [mesh1, inner.Xmesh, mesh2];
    [Xmesh, Nx] = my_unique(Xmesh, 1e-10);
    
    %% mesh size definition
    mesh_size_min = inner.Ymesh(2) - inner.Ymesh(1);
    %% mesh the bottom part of the system.pack    
    Num= floor( ( inner.yb - outer.yb )/mesh_size_min ) + 1;
    %%%this is for ticking out the unused element in the array
    mesh = -ones(1, Num);
    
    mesh(1) = inner.yb;
    mesh_size = mesh_size_min;
    i = 1;
    while( mesh_size < outer.Ygrid )
        if mesh(i) - mesh_size <= outer.yb
            mesh(i+1) = outer.yb;
            i = i+1;
            break;
        end
        mesh(i+1) = mesh(i) - mesh_size;
        mesh_size = mesh_size + mesh_size_min;
        i = i + 1;
    end
    Num = floor( (mesh(i) - outer.yb)/outer.Ygrid ) + 1;
    grid = (mesh(i) - outer.yb)/Num;
    rest = outer.yb : grid : mesh(i);
    
    mesh = [mesh, rest];
    id = mesh >= 0;
    mesh = mesh(id);
    
    mesh = sort(mesh);
    mesh1 = my_unique(mesh, 1e-10);
    
    %% mesh the top part of the system.pack
    Num= floor (( outer.yt - inner.yt )/mesh_size_min ) + 1;
    %%%this is for ticking out the unused element in the array
    mesh = -ones(1, Num);
    
    mesh(1) = inner.yt;
    mesh_size = mesh_size_min;
    i = 1;
    while( mesh_size < outer.Ygrid )
        if mesh(i) + mesh_size >= outer.yt
            mesh(i+1) = outer.yt;
            i = i+1;
            break;
        end
        mesh(i+1) = mesh(i) + mesh_size;
        mesh_size = mesh_size + mesh_size_min;
        i = i + 1;
    end
    Num = floor( (outer.yt - mesh(i))/outer.Ygrid ) + 1;
    grid = (outer.yt - mesh(i))/Num;
    rest = [mesh(i) : grid : outer.yt, outer.yt];
    
    mesh = [mesh, rest];
    id = mesh >= 0;
    mesh = mesh(id);
    
    mesh = sort(mesh);
    mesh2 = my_unique(mesh, 1e-10);
    
    %% merge the three mesh parts as Xmesh
    Ymesh = [mesh1, inner.Ymesh, mesh2];
    [Ymesh, Ny] = my_unique(Ymesh, 1e-10);
end

