function x = thermal_solve_ss(Y, H, heat, Ta, var)
%this is a thermal solver for steady state analysis
    disp('calculate temperature');
    tic;
    spparms('spumoni', 0);
    A = Y + H;
    A = 0.5* ( A + A' );
    Ha = H * ones(var, 1) * Ta;
    right = Ha + heat;
            
    x = A\right;
    toc;
    disp('calculatation done');
    disp(' ');
end

