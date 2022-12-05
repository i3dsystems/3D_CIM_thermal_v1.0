function [R, R_tran, S, S_tran] = thermal_fact(Y, C, H, dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    tic;
    fprintf('Choleskey factorization\n');
    A = Y + H + C/dt;
    A = (A+A')/2;
    [R, p, S] = chol(A);
    R_tran = R';
    S_tran = S';
    if (p~=0)
        error('Matrix is not SPD, please check');
    end
    fprintf('Choleskey factorization done, using %.2f seconds\n', toc);
end

