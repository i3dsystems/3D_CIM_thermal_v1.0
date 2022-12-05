function [b, lz_d] = hsp_id2chip_id(die, Layer, offset, chip, idx, idy, box)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    xl = box(1);
    yb = box(3);
    lz_d = Layer.thick( Layer.hsp + sum(Layer.N(1:die-1))+1);    
    const = 0;
    for i = 1:1:die-1
        const = const + (Layer.N(i)-1)*chip(i).Nx*chip(i).Ny;
    end
    b = ((idx-xl+1)+(idy-yb)*chip(die).Nx)+const+offset;
end