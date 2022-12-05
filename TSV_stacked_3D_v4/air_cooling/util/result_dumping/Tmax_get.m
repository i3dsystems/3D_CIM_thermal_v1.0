function T = Tmax_get(x, system, chip, Layer, draw, t)
    Nt = 0;
    for l =1:1:system.chip.N
        Nt = Nt + chip(l).N;
    end    
    T = zeros(Nt*2,1);
    pointer = 1;
    for i = 1:1:system.chip.N
        const = system.Nhsp;
        for l =1:1:i-1
            const = const + chip(l).Nx*chip(l).Ny*(Layer.N(l)-1);
        end
        for k=1:1:chip(i).N
            left = chip(i).Nx*chip(i).Ny*(1+(k-1)*chip(i).model)+1+const;
            right = left + chip(i).Nx*chip(i).Ny - 1;
            T(pointer) = max(x(left:right))-273;
            T(pointer+Nt) = min(x(left:right))-273;
            if draw.displayT == 1
                fprintf('chip%i die%i Min ~ Max: %.2f ~ %.2f\n', ...
                              i,    k,  T(pointer+Nt), T(pointer));
            end
            %% always output min and max temperature
            file_name_max = ['./result/chip', num2str(i), '_', 'die', num2str(k), '_max.txt'];
            file_name_min = ['./result/chip', num2str(i), '_', 'die', num2str(k), '_min.txt'];
            if t == 0
                fid=fopen(file_name_max,'w+');
                fprintf(fid,'%g\r\n',T(pointer));
                fclose(fid);
                fid=fopen(file_name_min,'w+');
                fprintf(fid,'%g\r\n',T(pointer+Nt));
                fclose(fid);
            else
                fid=fopen(file_name_max,'a+');
                fprintf(fid,'%g\r\n',T(pointer));
                fclose(fid);
                fid=fopen(file_name_min,'a+');
                fprintf(fid,'%g\r\n',T(pointer+Nt));
                fclose(fid);
            end
            if draw.write == 1
                file_name = ['./result/chip', num2str(i), '_', 'die', num2str(k), '.txt'];
                if t == 0
                    fid=fopen(file_name,'w+');
                    fprintf(fid,'%g\r\n',x(left : right)-273);
                    fclose(fid);                 
                else
                    fid=fopen(file_name,'a+');
                    fprintf(fid,'%g\r\n',x(left : right)-273);
                    fclose(fid);
                end
            end
            pointer = pointer+1;
        end
    end
end

