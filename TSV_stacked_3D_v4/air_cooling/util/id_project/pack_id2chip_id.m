function [t, lz_u] = pack_id2chip_id(die, offset, chip, idx, idy, box)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    xl = box(1);
    yb = box(3);
%     lz_u = Layer.thick( Layer.hsp + sum(Layer.N(1:die)));
%     const = 0;
%     for i = 1:1:die
%         const = const + (Layer.N(i)-1)*chip(i).Nx*chip(i).Ny;
%     end
    t = (idx-xl+1)+(idy-yb)*chip(die).Nx + ...
        offset - chip(die).Nx*chip(die).Ny;
end