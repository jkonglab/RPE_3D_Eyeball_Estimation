function gap = func_k(k, indep_var)
    l = indep_var(:, 1);
    R = indep_var(:, 2);
    gap = 2*pi * l - 2*pi * (R-k) .* sin(pi - l./R);
end

