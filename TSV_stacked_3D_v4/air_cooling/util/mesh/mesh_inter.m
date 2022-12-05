function system_ = mesh_inter(system)
%mesh the interposer using the results of system.box.Xmesh/Ymesh
%useful references: system.box.Xmesh, system.box.Ymesh
    [Nx, Xmesh, Ny, Ymesh] = mesh_box(system.box, system.inter);
    system.inter.Nx = Nx;
    system.inter.Ny = Ny;
    system.inter.Xmesh = Xmesh;
    system.inter.Ymesh = Ymesh;
    system_ = system;
    fprintf('interposer domain meshing done\n')    
end

