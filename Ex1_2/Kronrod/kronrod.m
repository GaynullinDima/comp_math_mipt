function [ x, w1, w2 ] = kronrod ( n, tol )


  m = floor ( ( n + 1 ) / 2 );
  even = ( 2 * m == n );

  d = 2.0;
  an = 0.0;
  for k = 1 : n
    an = an + 1.0;
    d = d * an / ( an + 0.5 );
  end
%
%  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
%
  tau(1) = ( an + 2.0 ) / ( an + an + 3.0 );
  b(m) = tau(1) - 1.0;
  ak = an;

  for l = 1 : m - 1

    ak = ak + 2.0;
    tau(l+1) = ( ( ak - 1.0 ) * ak ...
      - an * ( an + 1.0 ) ) * ( ak + 2.0 ) * tau(l) ...
      / ( ak * ( ( ak + 3.0 ) * ( ak + 2.0 ) ...
      - an * ( an + 1.0 ) ) );
    b(m-l) = tau(l+1);

    for ll = 1 : l
      b(m-l) = b(m-l) + tau(ll) * b(m-l+ll);
    end

  end

  b(m+1) = 1.0;
%
%  Calculation of approximate values for the abscissas.
%
  bb = sin ( 1.570796 / ( an + an + 1.0 ) );
  x1 = sqrt ( 1.0 - bb * bb );
  s = 2.0 * bb * x1;
  c = sqrt ( 1.0 - s * s );
  coef = 1.0 - ( 1.0 - 1.0 / an ) / ( 8.0 * an * an );
  xx = coef * x1;
%
%  Coefficient needed for weights.
%
%  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
%
  coef2 = 2.0 / ( 2 * n + 1 );
  for i = 1 : n
    coef2 = coef2 * 4.0 * i / ( n + i );
  end
%
%  Calculation of the K-th abscissa (a Kronrod abscissa) and the
%  corresponding weight.
%
  for k = 1 : 2 : n

    [ xx, w1(k) ] = kronrod_gauss ( n, m, tol, coef2, even, b, xx );
    w2(k) = 0.0;

    x(k) = xx;
    y = x1;
    x1 = y * c - bb * s;
    bb = y * s + bb * c;

    if ( k == n )
      xx = 0.0;
    else
      xx = coef * x1;
    end
%
%  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
%  corresponding weights.
%
    [ xx, w1(k+1), w2(k+1) ] = gauss ( n, m, tol, coef2, even, b, xx );

    x(k+1) = xx;
    y = x1;
    x1 = y * c - bb * s;
    bb = y * s + bb * c;
    xx = coef * x1;

  end
%
%  If N is even, we have one more Kronrod abscissa to compute,
%  namely the origin.
%
  if ( even )
    xx = 0.0;
    [ xx, w1(n+1) ] = kronrod_gauss ( n, m, tol, coef2, even, b, xx );
    w2(n+1) = 0.0;
    x(n+1) = xx;
  end

  return
end
