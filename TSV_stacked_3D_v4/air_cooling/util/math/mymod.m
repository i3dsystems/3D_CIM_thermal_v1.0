function [f, r]= mymod(x, y)
    x_ab = abs(x);    
    f = 0;
    if (y==0)
        r = 0;
        return;
    end
    if x_ab < y - 1e-12
        r_ab = x_ab;
    else
        while(x_ab >= y - 1e-12)
            f = f+1;
            x_ab = x_ab - y;
        end
        r_ab = x_ab;
    end
    
    if x<0
        r = x;
        f = -f;
    else
        r = r_ab;
    end
end

