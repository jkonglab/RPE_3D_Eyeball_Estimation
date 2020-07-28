function gap = func_r(R, indep_var)
    l = indep_var(:, 1);
    k = indep_var(:, 2);
    gap = 2*pi * l - 2*pi * (R-k) .* sin(pi - l./R);
end

