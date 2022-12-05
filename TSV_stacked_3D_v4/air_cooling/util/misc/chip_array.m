function [map, name] = chip_array(chip, N)
    map = zeros(N, 6);
    name = cell(N,1);
    for i=1:N
        map(i, :) = [chip(i).xl, chip(i).yb, chip(i).Xsize, chip(i).Ysize, sum(chip(i).power), 0];
        %
        name(i) = chip(i).name;
    end
end

