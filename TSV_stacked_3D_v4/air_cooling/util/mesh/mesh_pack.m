function system_ = mesh_pack(system)
%mesh the package using the results of interposer meshing
    [Nx, Xmesh, Ny, Ymesh] = mesh_box(system.inter, system.pack);
    system.pack.Nx = Nx;
    system.pack.Ny = Ny;
    system.pack.Xmesh = Xmesh;
    system.pack.Ymesh = Ymesh;
    system_ = system;    
    fprintf('package domain meshing done\n')
end

