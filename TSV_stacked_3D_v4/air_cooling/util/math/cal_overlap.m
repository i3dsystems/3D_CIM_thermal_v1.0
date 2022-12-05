function [overlap, array]= cal_overlap( boundary1, boundary2 )
    xl1 = boundary1(1);
    xr1 = boundary1(2);
    yb1 = boundary1(3);
    yt1 = boundary1(4);
    
    xl2 = boundary2(1);
    xr2 = boundary2(2);
    yb2 = boundary2(3);
    yt2 = boundary2(4);    

    if xr1<=xl2 || xl1>=xr2 || yb1>=yt2 || yt1<=yb2
        overlap = 0;
        array = [];
    else
        x = sort([xl1 xr1 xl2 xr2]);
        y = sort([yb1 yt1 yb2 yt2]);
        array = [x(2), x(3), y(2), y(3)];
        overlap = (x(3)-x(2))*(y(3)-y(2));
    end
end

