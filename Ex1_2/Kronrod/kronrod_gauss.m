function [ x, w ] = kronrod_gauss ( n, m, tol, coef2, even, b, x )

  if ( x == 0.0 )
    ka = 1;
  else
    ka = 0;
  end
%
%  Iterative process for the computation of a Kronrod abscissa.
%
  for iter = 1 : 50

    b1 = 0.0;
    b2 = b(m+1);
    yy = 4.0 * x * x - 2.0;
    d1 = 0.0;

    if ( even )
      ai = m + m + 1;
      d2 = ai * b(m+1);
      dif = 2.0;
    else
      ai = m + 1;
      d2 = 0.0;
      dif = 1.0;
    end

    for k = 1 : m
      ai = ai - dif;
      i = m - k + 1;
      b0 = b1;
      b1 = b2;
      d0 = d1;
      d1 = d2;
      b2 = yy * b1 - b0 + b(i);
      if ( ~even )
        i = i + 1;
      end
      d2 = yy * d1 - d0 + ai * b(i);
    end

    if ( even )
      f = x * ( b2 - b1 );
      fd = d2 + d1;
    else
      f = 0.5 * ( b2 - b0 );
      fd = 4.0 * x * d2;
    end
%
%  Newton correction.
%
    delta = f / fd;
    x = x - delta;

    if ( ka == 1 )
      break
    end

    if ( abs ( delta ) <= tol )
      ka = 1;
    end

  end
%
%  Catch non-convergence.
%
  if ( ka ~= 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'ABWE1 - Fatal error!\n' );
    fprintf ( 1, '  Iteration limit reached.\n' );
    fprintf ( 1, '  Last DELTA was %e\n', delta );
    error ( 'ABWE1 - Fatal error!' );
  end
%
%  Computation of the weight.
%
  d0 = 1.0;
  d1 = x;
  ai = 0.0;
  for k = 2 : n
    ai = ai + 1.0;
    d2 = ( ( ai + ai + 1.0 ) * x * d1 - ai * d0 ) / ( ai + 1.0 );
    d0 = d1;
    d1 = d2;
  end

  w = coef2 / ( fd * d2 );

  return
end
