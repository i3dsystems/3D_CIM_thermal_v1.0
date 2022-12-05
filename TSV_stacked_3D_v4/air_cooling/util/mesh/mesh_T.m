function [Len, mesh] = mesh_T(chip, system)
%this function is to generate time domain non-uniform meshes
    simtime = tic;
    fprintf('Time domain meshing: ');
    Len = system.T/system.dt+1;
    mesh = linspace(0, system.T, Len);
%     Nmax = floor(system.T/system.dt) + 1;
%     T = zeros(Nmax, 1);
%     l = 1;
%     for i = 1:system.chip.N
%         filename = ['chip', num2str(i), '.trace'];
%         if exist(filename, 'file') == 2
%             temp = load(filename);
%             unit_len = 3+sum(chip(i).blk_num);
%             len = length(temp)/unit_len;
%             T(l:l+len-1) = temp(1:unit_len:end);
%             l = l+len;
%         end
%     end
%     T = T(1:l-1);
%     T = my_unique(T, 1e-10);
%     
%     len = length(T);
%     
%     step = system.dt;
%     step_su = system.su;
%     threshold = system.dt * system.su_limit;
%     
%     Tfinal = zeros(Nmax, 1);
%     Tfinal(1) = T(1);
%     j = 2;
%     for i = 2:1:len
%         if T(i) >= system.T
%             Tfinal(j) = system.T;
%             j = j+1;
%             break;
%         end
%         if T(i) - T(i-1) > step
%             Tfinal(j) = T(i-1) + step;
%             Cur_step = min(step + step*step_su, threshold);
%             while(T(i) - Tfinal(j) > Cur_step)
%                 Tfinal(j+1) = Tfinal(j) + Cur_step;
%                 j = j+1;
%                 Cur_step = min(Cur_step + step_su*step, threshold);
%             end
%             if Tfinal(j) < T(i)
%                 Tfinal(j+1) = T(i);
%             else
%                 Tfinal(j) = T(i);
%             end
%             j = j+1; 
%         else
%             Tfinal(j) = T(i);
%             j = j+1;
%         end
%     end
%     mesh = Tfinal(1:j-1,1);
%     Len = length(mesh);
    fprintf('takes %.3f seconds\n', toc(simtime));
    fprintf('\n');
    fid=fopen('Tfinal.txt','w+');
    fprintf(fid,'%g\r\n',mesh);
    fclose(fid);    
end

