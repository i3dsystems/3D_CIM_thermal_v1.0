function array = ChipInPlane(chip, big_plane)
% This function find the corresponding position of 
    xl = find( abs(chip.xl - big_plane.Xmesh)<1e-12 );
    xr = find( abs(chip.xr - big_plane.Xmesh)<1e-12 );
    yb = find( abs(chip.yb - big_plane.Ymesh)<1e-12 );
    yt = find( abs(chip.yt - big_plane.Ymesh)<1e-12 );
    array = [xl, xr, yb, yt];
end