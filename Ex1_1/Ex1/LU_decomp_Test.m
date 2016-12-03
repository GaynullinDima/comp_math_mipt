% First test.

A1 = [ 0, 2, 0; 1, 0, 1; 0, 0, 1 ];
b1 = [ 1; 2; 3 ];

sol1 = Solver ( A1, b1 )

% Compare first results to true.

[ L1, U1 ] = lu ( A1 );
y = L1 \ b1;
x = U1 \ y

% Compare first results to classical.

Alu1 = ludec ( A1 );
sol_cl1 = lusol ( Alu1, b1 )


% =========================================================================


% Second test.

A2 = [ 0, 1, 0; 0, 0, 1; 1, 0, 0 ];
b2 = [ -1; 2; -5 ];

sol2 = Solver ( A2, b2 )

% Compare second results to true.

[ L2, U2 ] = lu ( A2 );
y = L2 \ b2;
x = U2 \ y

% Compare second results to classical.

Alu2 = ludec ( A2 );
sol_cl2 = lusol ( Alu2, b2 )


% =========================================================================


% Third test.

n = 10;

b3 = 1:n;

A3 = eye(n);
A3(1,:) = [];
A3(n,1) = 2;
A3(n,2:n) = b3(2:n);

b3 = b3';

sol3 = Solver ( A3, b3 )

% Compare third results to true.

[ L3, U3 ] = lu ( A3 );
y = L3 \ b3;
x = U3 \ y

% Compare second results to classical.

Alu3 = ludec ( A3 );
sol_cl3 = lusol ( Alu3, b3 )


% =========================================================================


% Check good matrix.

A4 = [ 1, 0, 0; 0, 1, 0; 0, 0, 1 ];
b4 = [ -1; 2; -5 ];

sol4 = Solver ( A4, b4 )

% Compare second results to true.

[ L4, U4 ] = lu ( A4 );
y = L4 \ b4;
x = U4 \ y

% Compare second results to classical.

Alu4 = ludec ( A4 );
sol_cl4 = lusol ( Alu4, b4 )