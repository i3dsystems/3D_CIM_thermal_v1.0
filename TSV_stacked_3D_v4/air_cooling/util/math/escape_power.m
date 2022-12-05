function escape_power(hsp, Tmap, Ta, h)
%given boundary condition, temperature distribution
%calculate how much power escaping from the boundary
    power = 0;
    for i = 1 : hsp.Nx
        for j=1 : hsp.Ny
            if(i>1) 
                x1 = hsp.Xmesh(i) - hsp.Xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<hsp.Nx)
                x2 = hsp.Xmesh(i+1) - hsp.Xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = hsp.Ymesh(j)- hsp.Ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<hsp.Ny)
                y2 = hsp.Ymesh(j+1) - hsp.Ymesh(j);
            else
                y2 = 0;
            end
            xgrid = (x1+x2)/2;
            ygrid = (y1+y2)/2;
            power = power + xgrid*ygrid*h*(Tmap(j,i)-Ta);
        end
    end
    fprintf('total power escaped from top surface: %f\n', power);
end

