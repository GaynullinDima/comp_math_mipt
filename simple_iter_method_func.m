function [ t, t_max, iter, u ] = simple_iter_method_func( A, f, epsilon, u, t )

    u_sol = A \ f;
    A_new = A' * A;
    f_new = A' * f;
    e = eig(A_new);
    last_elem = length(e);
    
    if nargin < 5 || isempty(t)
       
        t = 2 / (e(1) + e(last_elem));
        
    end;

    t_max = 2 / e(last_elem);
    
    E = eye(size(A)(1));
    for iter = 1:100000
        if abs(norm(u - u_sol)) <= epsilon
            return;
        end;
        u = (E - t * A_new) * u + t * f_new;
    end;

    u = NaN;
    error('Too many iterations');
end
