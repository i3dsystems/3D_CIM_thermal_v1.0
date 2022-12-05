function heat = add_heat(chip, system, connect_plane, var, draw_P, Layer)
    fprintf('assign power map\n');
    tic;
    heat = zeros(var,1);
    
    if draw_P == 1
        drawP_pack = NaN*ones(connect_plane.Ny, connect_plane.Nx);
    end
    
    for m = 1 : 1 : system.chip.N
    %caculate power
        Chip_xl = find(abs(connect_plane.Xmesh-chip(m).xl)<1e-12);
        Chip_yb = find(abs(connect_plane.Ymesh-chip(m).yb)<1e-12);
        
        die_num = chip(m).N;
        if size(chip(m).map,1) > 0
            bgrd_block = zeros(size(chip(m).map,1), 1);
        end
        gridNx_chip = chip(m).Nx;
        gridNy_chip = chip(m).Ny;
        const = gridNx_chip * gridNy_chip;
        index_start = system.Nhsp;
        for l= 1 : 1 : m-1
            index_start = index_start + chip(l).Nx*chip(l).Ny.*(Layer.N(l)-1);
        end
        
        chip_xmesh = chip(m).Xmesh;
        chip_ymesh = chip(m).Ymesh;

        if draw_P == 1
            drawP_die = zeros(chip(m).Ny, chip(m).Nx, chip(m).N);
        end
        
        if size(chip(m).map, 1) > 0
            density = zeros(size(chip(m).map, 1), 1);
        end
        
        for DIE = 1:1:die_num
            metal = chip(m).p_metal(DIE);
            background = 0;
            if DIE == 1
                start = 1;
            else
                start = sum(chip(m).blk_num(1:DIE-1))+1;
            end
            if DIE == die_num
                End = sum(chip(m).blk_num(1:DIE));
            else
                End = start + chip(m).blk_num(DIE)-1;
            end

            %calculate whether need to assign background power        
            if End-start < 0
                background = chip(m).power(DIE)/(chip(m).Xsize*chip(m).Ysize);
            else
                if abs(chip(m).power(DIE) - sum(chip(m).map(start:End, 5))) > 1e-8
                    power_back = chip(m).power(DIE) - sum(chip(m).map(start:End, 5));
                    tmp = chip(m).Xsize*chip(m).Ysize;
                    for i = start:End
                        if chip(m).map(i, 6) == 0
                            tmp = tmp - chip(m).map(i, 3)*chip(m).map(i, 4);
                        end
                    end
                    if tmp <= 10e-12
                        area_back = chip(m).Xsize*chip(m).Ysize;
                    else
                        area_back = tmp;
                    end
                    background = power_back / area_back;
                end
            end
            %disp(background);
            if background > 1e-8 
                for i = 1:1:gridNx_chip
                    for j=1:1:gridNy_chip
                        id = i+(j-1)*gridNx_chip+const*(1+(DIE-1)*chip(m).model)+index_start;
                        if i == 1
                            boundary(1) = chip_xmesh(i);
                        else
                            boundary(1) = (chip_xmesh(i-1)+chip_xmesh(i))/2;
                        end
                        if i == gridNx_chip
                            boundary(2) = chip_xmesh(i);
                        else
                            boundary(2) = (chip_xmesh(i)+chip_xmesh(i+1))/2;
                        end

                        if j == 1
                            boundary(3) = chip_ymesh(j);
                        else
                            boundary(3) = (chip_ymesh(j-1)+chip_ymesh(j))/2;
                        end
                        if j == gridNy_chip
                            boundary(4) = chip_ymesh(j);
                        else
                            boundary(4) = (chip_ymesh(j)+chip_ymesh(j+1))/2;
                        end

                        gridx = boundary(2) - boundary(1);
                        gridy = boundary(4) - boundary(3);
                        grid_area = gridx*gridy;
                        heat(id) = heat(id) + background*grid_area*(1-metal);
                        heat(id+const) = heat(id+const) + background*grid_area*metal;
                        if draw_P == 1
                            drawP_die(j,i,DIE) = heat(id)/(gridx*gridy*10^4);
                            if DIE == die_num
                                drawP_pack(j-1+Chip_yb, i-1+Chip_xl) = heat(id)/(gridx*gridy*10^4);
                            end
                        end
                    end
                end  
            end

            %assign the block power maps
            if End < start
                continue
            end
            for k1 = start:End
                if chip(m).map(k1, 6) ~= 0
                    index = chip(m).map(k1, 6);
                    bgrd_block(index) = bgrd_block(index)+chip(m).map(k1,3)*chip(m).map(k1,4);
                end
            end
            for k1 = start:End
                if bgrd_block(k1) > 0
                    bgrd_block(k1) = chip(m).map(k1,5)/(chip(m).map(k1,3)*chip(m).map(k1,4) - bgrd_block(k1));
                else
                    bgrd_block(k1) = chip(m).map(k1,5)/(chip(m).map(k1,3)*chip(m).map(k1,4));
                end   
            end
            map = chip(m).map(start:End,:);
            blk_name = chip(m).blk_name(start:End);
            for k=start:1:End
                if chip(m).map(k, 6) == 0
                    density(k) = bgrd_block(k) - background;
                else
                    index = chip(m).map(k, 6);
                    density(k) = bgrd_block(k) - bgrd_block(index);
                end
                xl = sum(sum(chip_xmesh<chip(m).map(k,1)));
                if xl <= 0
                    xl = 1;
                end
                xr = sum(sum(chip_xmesh<chip(m).map(k,3)+chip(m).map(k,1)))+1;
                if xr >= gridNx_chip
                    xr = gridNx_chip;
                end 
                yb = sum(sum(chip_ymesh<chip(m).map(k,2)));
                 if yb <= 0
                    yb = 1;
                end           
                yt = sum(sum(chip_ymesh<chip(m).map(k,4)+chip(m).map(k,2)))+1;
                if yt >= gridNy_chip
                    yt = gridNy_chip;
                end            
                boundary_blk = [chip(m).map(k,1) chip(m).map(k,3)+chip(m).map(k,1) ...
                                chip(m).map(k,2) chip(m).map(k,4)+chip(m).map(k,2)];
                for i=xl:1:xr
                    for j=yb:1:yt                    
                        id = i+(j-1)*gridNx_chip+const*(1+(DIE-1)*chip(m).model)+index_start;
                        if i == 1
                            boundary(1) = chip_xmesh(i);
                        else
                            boundary(1) = (chip_xmesh(i-1)+chip_xmesh(i))/2;
                        end
                        if i == gridNx_chip
                            boundary(2) = chip_xmesh(i);
                        else                        
                            boundary(2) = (chip_xmesh(i)+chip_xmesh(i+1))/2;
                        end

                        if j == 1
                            boundary(3) = chip_ymesh(j);
                        else
                            boundary(3) = (chip_ymesh(j-1)+chip_ymesh(j))/2;
                        end
                        if j == gridNy_chip
                            boundary(4) = chip_ymesh(j);
                        else
                            boundary(4) = (chip_ymesh(j)+chip_ymesh(j+1))/2;
                        end                                        
                        grid_area = cal_overlap (boundary, boundary_blk);
                        heat(id) = heat(id) + density(k)*grid_area*(1-metal);    
                        heat(id+const) = heat(id+const) + density(k)*grid_area*metal;    
                        if draw_P == 1
                            gridx = boundary(2) - boundary(1);
                            gridy = boundary(4) - boundary(3);                   
                            drawP_die(j,i,DIE) = heat(id)/(gridx*gridy*10^4);
                            if DIE == die_num
                                drawP_pack(j-1+Chip_yb, i-1+Chip_xl) = heat(id)/(gridx*gridy*10^4);
                            end
                        end
                    end
                end
            end
            if draw_P == 1                        
                if m == 1
                    start = 0;
                else
                    start = 0;
                    for l= 1 : 1 : m-1
                        start = start + chip(l).N;
                    end
                end
                if chip(m).blk_num(DIE) > 0
                    figure(start+DIE+10);
                    h = pcolor(chip(m).Xmesh*100, chip(m).Ymesh*100,drawP_die(:,:,DIE));%, max(chip(m).blk_num)*2,'Linestyle','none');

                    grid off;
                    set(h, 'EdgeColor', 'none');
                    h=colorbar;
                    set(get(h,'Title'),'string','W/cm2','FontSize',16)
                    set(gca,'FontSize',16);
                    xlabel('x(cm)');
                    ylabel('y(cm)');set(gca,'FontSize',16);
                    axis equal;
                    axis off;
                    string = ['chip', num2str(m),': die',num2str(DIE),' power map'];
                    title(string);                    
                    hold on;
                    len = size(map, 1);
                    for i=1:len
                        xl = map(i,1);
                        width = map(i,3);
                        yb = map(i,2);
                        height = map(i,4);                    
                        rectangle('Position',[xl yb width height]*100);
                        if isempty(char(blk_name(i))) == 0                            
                            text((xl+width/2)*100, (yb+height/2)*100, char(blk_name(i)), 'HorizontalAlignment','center', 'FontSize', 12)
                        end
                        hold on;
                    end                    
                end
            end
        end
    end
    toc;

    if draw_P == 1 && system.chip.N > 1
        figure(system.Ndie+10+1);
        h = pcolor(connect_plane.Xmesh*100, connect_plane.Ymesh*100, drawP_pack);
        set(gca,'position',[0 0 1 1],'units','normalized')           
        grid off;
        set(h, 'EdgeColor', 'none'); 
        h=colorbar;        
        set(get(h,'Title'),'string','W/cm2','FontSize',16)
        set(gca,'FontSize',16);
        xlabel('x(cm)');
        ylabel('y(cm)');set(gca,'FontSize',16);
        caxis([0, 130]);
        axis equal;
        axis off;
        %title('power density of bottom dice in package view');
        
        hold on;
        [accu_map, blk_name] = chip_map_stack(chip, system.chip.N);
        len = size(accu_map, 1);
        for i=1:len
            xl = accu_map(i,1);
            width = accu_map(i,3);
            yb = accu_map(i,2);
            height = accu_map(i,4);                    
            rectangle('Position',[xl yb width height]*100);
            if isempty(char(blk_name(i))) == 0                            
                text((xl+width/2)*100, (yb+height/2)*100, char(blk_name(i)), 'HorizontalAlignment','center', 'FontSize', 12)
            end
            hold on;
        end   
    end

    for m = 1:system.bridge.N
        xl = find(abs(connect_plane.Xmesh-system.bridge.xl(m))<1e-12);
        xr = find(abs(connect_plane.Xmesh-system.bridge.xr(m))<1e-12);
        yb = find(abs(connect_plane.Ymesh-system.bridge.yb(m))<1e-12);   
        yt = find(abs(connect_plane.Ymesh-system.bridge.yt(m))<1e-12);

        density = system.bridge.power / ((system.bridge.xr(m)-system.bridge.xl(m))*(system.bridge.yt(m)-system.bridge.yb(m)));

        for i = xl:xr
            for j = yb:yt
                id = system.Nchip + i + (j-1)*connect_plane.Nx;
                dx = connect_plane.Xmesh(i+1) - connect_plane.Xmesh(i);
                dy = connect_plane.Ymesh(j+1) - connect_plane.Ymesh(j);
                heat(id) = density * dx * dy + heat(id);
            end
        end
    end
    
    fprintf('total power: %f\n', sum(heat));
    fprintf('finish assigning power map\n');
    fprintf('\n');
end

