function die = box2chip( idx, idy, box)
%four coordinates stored in box array
%determine whether (idx, idy) falls in any one from box
%box arrangement: #die id: xl xr yb yt
    N_die = size(box, 1);
    die = 0;
    for i= 1 : N_die
        a = box(i,:);
        xl = a(1); xr = a(2); yb = a(3); yt = a(4);
        if idx >= xl && idx <= xr && idy >= yb && idy <= yt
            die = i;
            break;
        end
    end
end

