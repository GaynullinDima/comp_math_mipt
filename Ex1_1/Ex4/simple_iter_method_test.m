A = [3, 2; 4, 3];
f = [1; 1];
acc = 1e-2;
u0 = [0; 0];
[t_opt, t_max, k, u_iter] = simple_iter_method_func(A, f, acc, u0);
u_sol = A \ f;


t_vect = linspace(0.005, t_max-0.00001, 100);
t_len = size(t_vect);
k_vect = zeros(t_len);
for i = 1:t_len(2)
    [t_opt_var, t_max_var, k_var, u_iter_var] = simple_iter_method_func(A, f, acc, u0, t_vect(i));
    k_vect(i) = k_var;
end;

plot(t_vect,k_vect)
grid on