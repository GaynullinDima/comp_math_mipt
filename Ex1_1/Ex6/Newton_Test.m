%x = linspace ( -1, 1, 1000 );
func = @(x) -x.^2;
dfunc = @(x) -2*x;
x0 = -0.25;
acc = 1e-3;

[ sol, iter ] = NewtonWithoutDiff ( func, acc, x0, -0.5, 0.5, 300, 2 )
[sol, iter] = Newton ( func, dfunc, acc, x0, -0.5, 0.5, 300, 2 )