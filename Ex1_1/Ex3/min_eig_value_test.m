A = [5, 1; 3, 2];
epsilon = 1e-2;
v = [1; 0];
[ min, vect_min, iter ] = min_eig_value_func( A, v, epsilon )