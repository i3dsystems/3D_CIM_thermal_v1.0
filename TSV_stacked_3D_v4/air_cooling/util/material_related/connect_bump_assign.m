function  connect_bump_map = connect_bump_assign(chip, connect_bump_array)
%this function is to find the MFI count and store it in a 2-d array
    N = size(connect_bump_array, 1);
    connect_bump_map = zeros(chip.Ny, chip.Nx, 2);
    connect_bump_x = zeros(chip.Ny, chip.Nx);
    connect_bump_y = zeros(chip.Ny, chip.Nx);
    for i=1:N
        % this is important since via_count_assign use relative coordinates
        bump.xl = connect_bump_array(i,1)-chip.xl;bump.xr = connect_bump_array(i,2)-chip.xl;
        bump.yb = connect_bump_array(i,3)-chip.yb;bump.yt = connect_bump_array(i,4)-chip.yb;
        bump.px = connect_bump_array(i,5);bump.py = connect_bump_array(i,6);
        bump.d = connect_bump_array(i,7);
        via_map = via_count_assign(chip, bump);
        via_map_x = via_map(:,:,1); via_map_y = via_map(:,:,2);
        connect_bump_x (via_map_x~=0) = via_map_x (via_map_x ~= 0);
        connect_bump_y (via_map_y~=0) = via_map_y (via_map_y ~= 0);
    end
    connect_bump_map(:,:,1) = connect_bump_x;
    connect_bump_map(:,:,2) = connect_bump_x;    
end

