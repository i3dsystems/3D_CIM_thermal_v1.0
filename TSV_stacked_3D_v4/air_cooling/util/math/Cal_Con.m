function [K_V, K_L ] = Cal_Con(K_eff, t)
    t_total = sum(t);
    K_L = K_eff*t'/t_total;
    K_V = t_total/(t*(1./K_eff)');
end

