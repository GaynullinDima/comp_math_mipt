function [ min, vect_min, iter ] = min_eig_value_func( A, v, epsilon )
%   Calculating precise eigenvalue 
    sol = eig(A);
    min_sol = sol(1);

%   Calculating approximate values
    for iter = 1:1000
        z = A \ v;
        z_norm = norm(z);
        v = z ./ z_norm;
        if abs( 1 / z_norm - min_sol ) <= epsilon
            vect_min = z ./ z_norm;
            min = 1 / z_norm;
            return;
        end;
    end;

%   Method doesn't converge 
    min = NaN;
    vect_min = NaN;
    error('Too many iterations');
end