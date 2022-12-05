function [chip_, system_] = mesh(chip, system)
%%%%all coordinates are refering to the heat spreader
    fprintf('meshing the stack\n');
    tic;
    %% build whole coordinates system for easy reference    
    % meshing the chip
    [chip, system] = block_define(chip, system);
    [chip, system] = mesh_chip(chip, system);
    %% package meshing
    system = mesh_pack_no_inter(system);
    %the results are written into system.pack.Xmesh/Ymesh
    %% heat spreader meshing
    system = mesh_hsp_FOWLP(system);
    
    system_ = system;
    chip_ = chip;
    toc; 
    fprintf('meshing done\n');
    fprintf('\n')      
end

