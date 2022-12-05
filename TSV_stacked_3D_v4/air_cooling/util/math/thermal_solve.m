function x = thermal_solve(Y, C, H, heat, prev_x, prev_heat, Ta, var, dt)
%this is a thermal solver for transient state analysis
    disp('calculate temperature');
    tic;
    spparms('spumoni', 0);
    
    %this is a backward euler scheme
    A = Y + H + C/dt;
    Ha = H * ones(var, 1) * Ta;
    right = Ha + heat + (C/dt)*prev_x;
    
    %this is a trapezoid scheme
%     A = Y/2 + H/2 + C/dt;
%     Ha = H * ones(var, 1) * Ta;
%     right = Ha + (heat + prev_heat)/2 + (C/dt)*prev_x - ...
%             H/2*prev_x - Y/2*prev_x;
    
    A = (A+A')/2;
    x = A\right;
    toc;
    disp('calculatation done');
    disp(' ');

end

