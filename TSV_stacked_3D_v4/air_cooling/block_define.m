function [chip_, system_] = block_define(chip, system)
    %%this function is used for chip placement check
    for i = 1 : 1 : system.chip.N
        chip(i).xr = chip(i).xl + chip(i).Xsize;
        chip(i).yt = chip(i).yb + chip(i).Ysize;
        if size(chip(i).map, 1) > 0
            chip(i).map(:,1) = chip(i).map(:,1) + chip(i).xl;
            chip(i).map(:,2) = chip(i).map(:,2) + chip(i).yb;
        end
    end
    system.pack.xr = system.pack.xl + system.pack.Xsize;
    system.pack.yt = system.pack.yb + system.pack.Ysize;
    
    system.hsp.xl = 0; system.hsp.xr = system.hsp.Xsize;
    system.hsp.yb = 0; system.hsp.yt = system.hsp.Ysize;
    
    %% find the box to cover the chip domains
    xl_all = zeros(system.chip.N, 1);
    xr_all = zeros(system.chip.N, 1);
    yb_all = zeros(system.chip.N, 1);
    yt_all = zeros(system.chip.N, 1);
    for i = 1: 1: system.chip.N
        xl_all(i) = chip(i).xl;
        xr_all(i) = chip(i).xr;
        yb_all(i) = chip(i).yb;
        yt_all(i) = chip(i).yt;
    end
    xl = min(xl_all); xr = max(xr_all);
    system.box.x = sort([xl_all; xr_all]);    
    system.box.x = my_unique(system.box.x, 1e-10);
    yb = min(yb_all); yt = max(yt_all);
    system.box.y = sort([yb_all; yt_all]);
    system.box.y = my_unique(system.box.y, 1e-10);
    
    %% chip overlapping check
    for i = 1 : 1: system.chip.N
        for j = i+1:1:system.chip.N
            if chip(i).layer.id == chip(j).layer.id
                blk1 = [chip(i).xl, chip(i).xr, chip(i).yb, chip(i).yt];
                blk2 = [chip(j).xl, chip(j).xr, chip(j).yb, chip(j).yt];
                area = cal_overlap(blk1, blk2);
                if area > 0
                    msg = 'chip overlaps checks chip'+num2str(i)+' and chip'+num2str(j);
                    error(msg);
                end
            end
        end
    end
    fprintf('The locations of the chip are fine \n');
    %% interposer and hear spreader check
    if xl < system.pack.xl || xr > system.pack.xr || ...
       yb < system.pack.yb || yt > system.pack.yt
        error('interposer cannot hold all the chips, increase the interposer size');
    end
    fprintf('The locations between the chip and package are fine \n');
    
    if xl <= system.hsp.xl || xr >= system.hsp.xr || ...
       yb <= system.hsp.yb || yt >= system.hsp.yt
        error('heat spreader cannot hold all the chips, increase the interposer size \n');
    end
    fprintf('The locations between the chip and heat spreader are fine \n');
    
    system.box.xl = xl; system.box.xr = xr;
    system.box.yb = yb; system.box.yt = yt;
    
    chip_ = chip;
    system_ = system;
end

