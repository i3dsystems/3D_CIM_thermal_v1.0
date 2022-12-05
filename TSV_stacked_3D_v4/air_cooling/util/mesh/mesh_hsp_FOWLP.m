function system_ = mesh_hsp_FOWLP(system)
%mesh the heat spreader using the results of system.pack.Xmesh/Ymesh
    [Nx, Xmesh, Ny, Ymesh] = mesh_box(system.pack, system.hsp);
    system.hsp.Nx = Nx;
    system.hsp.Ny = Ny;
    system.hsp.Xmesh = Xmesh;
    system.hsp.Ymesh = Ymesh;
    system_ = system;    
    fprintf('heat spreader domain meshing done\n')
end

