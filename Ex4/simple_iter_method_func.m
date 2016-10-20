function [ t, t_max, k, u ] = simple_iter_method_func( A, f, epsilon, u, t )

    u_sol = A \ f;
    e = eig(A);
    last_elem = length(e);
    
    if nargin < 5 || isempty(t)
        t = 2 / (e(1) + e(last_elem));
    end;

    t_max = 2 / e(last_elem);
    
    E = eye(size(A));
    k = 0;
    while abs(norm(u - u_sol)) > epsilon
        u = (E - t * A) * u + t * f;
        k = k + 1;
    end;
end
