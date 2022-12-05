function [accu_map,accu_name] = chip_map_stack(chip, N)
    N_block = 0;
    Ncol = 0;
    for i=1:N
        if size(chip(i).blk_num, 1) > 0
            N_block = chip(i).blk_num(end)+N_block;
            Ncol = max(size(chip(i).map, 2), Ncol);
        end
    end
    accu_map = zeros(N_block, Ncol);
    accu_name = cell(N_block, 1);
    l = 0;
    for i=1:N
        if size(chip(i).blk_num, 1) > 0
            End = sum(chip(i).blk_num);
            Stt = End - chip(i).blk_num(end) + 1;
            len = chip(i).blk_num(end);
            if len == 0
                len = len + 1;
                accu_map(l+1:l+len, :) = [chip(i).xl, chip(i).yb, chip(i).Xsize, chip(i).Ysize, 0, 0];
                accu_name(l+1:l+len) = chip(i).name;
            else
                accu_map(l+1:l+len, :) = chip(i).map(Stt:End, :);
                accu_name(l+1:l+len) = chip(i).blk_name(Stt:End);
            end
            l = l + len;
        end
    end
end

