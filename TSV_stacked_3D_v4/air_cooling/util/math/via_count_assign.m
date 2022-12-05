function via_map = via_count_assign(chip, via)
%This function is to assign via to each meshes
    via_map = zeros(chip.Ny, chip.Nx, 2);
    %1 stores x direction via number
    %2 stores y direction via number
    via_xl = via.xl + via.px/2;
    via_yb = via.yb + via.py/2;
    via_x = via.px;
    via_y = via.py;
    if via_x == 0 || via_y == 0
        return;
    end
    Nx = mymod((via.xr - via.px/2 - via_xl), via_x)+1;
    Ny = mymod((via.yt - via.px/2 - via_yb), via_y)+1;

    x = via_xl + (Nx - 1)*via_x + via.d/2 + chip.xl;
    y = via_yb + (Ny - 1)*via_y + via.d/2 + chip.yb;
    xr = sum( x >= chip.Xmesh);
    yt = sum( y >= chip.Ymesh);

    x = via_xl + via.d/2 + chip.xl;
    y = via_yb + via.d/2 + chip.yb;
    xl = sum( x >= chip.Xmesh);
    yb = sum( y >= chip.Ymesh);        

    for i = 1:Nx
        x = via_xl + (i - 1)*via_x + via.d/2 + chip.xl;
        xid = sum( x>=chip.Xmesh);
        via_map(yb:yt, xid, 1) = 1;%(yb:yt, xid, 1) + 1;
    end

    for i = 1:Ny
        y = via_yb + (i - 1)*via_y + via.d/2 + chip.yb;
        yid = sum( y>=chip.Ymesh);
        via_map(yid, xl:xr, 2) = 1;%via_map(yid, xl:xr, 2) + 1;
    end    

    fprintf('there are %i vias\n', sum(sum(via_map(:,:,1).*via_map(:,:,2))));
end

