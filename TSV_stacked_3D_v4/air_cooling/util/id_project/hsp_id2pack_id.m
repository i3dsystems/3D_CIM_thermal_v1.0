function [b, lz_d] = hsp_id2pack_id(Layer, offset, pack, idx, idy, box)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    xl = box(1);
    yb = box(3);
    lz_d = Layer.thick( Layer.hsp + 1 );    
    b = ((idx-xl+1)+(idy-yb)*pack.Nx)+offset;
end