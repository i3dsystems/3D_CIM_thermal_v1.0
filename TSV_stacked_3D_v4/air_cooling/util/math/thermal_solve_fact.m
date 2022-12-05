function x = thermal_solve_fact(R, R_tran, S, S_tran, Y, C, H, heat, prev_x, prev_heat, Ta, var, dt)
%this is a thermal solver for transient state analysis
    disp('calculate temperature');
    tic;
    spparms('spumoni', 0);
    
    %this is a backward euler scheme
    Ha = H * ones(var, 1) * Ta;
    b = Ha + heat + (C/dt)*prev_x;
    
    %this is a trapezoid scheme
%     A = Y/2 + H/2 + C/dt;
%     Ha = H * ones(var, 1) * Ta;
%     right = Ha + (heat + prev_heat)/2 + (C/dt)*prev_x - ...
%             H/2*prev_x - Y/2*prev_x;
    x = S*(R\(R_tran\(S_tran*b)));
    toc;
    fprintf('calculatation done, using %.2f \n\n', toc);
end

