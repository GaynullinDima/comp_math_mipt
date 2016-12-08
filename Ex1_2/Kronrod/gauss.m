function [ x, w1, w2 ] = gauss ( n, m, tol, coef2, even, b, x )

  if ( x == 0.0 )
    ka = 1;
  else
    ka = 0;
  end
%
%  Iterative process for the computation of a Gaussian abscissa.
%
  for iter = 1: 50

    p0 = 1.0;
    p1 = x;
    pd0 = 0.0;
    pd1 = 1.0;
%
%  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
%
    if ( n <= 1 )
      if ( eps < abs ( x ) )
        p2 = ( 3.0 * x * x - 1.0 ) / 2.0;
        pd2 = 3.0 * x;
      else
        p2 = 3.0 * x;
        pd2 = 3.0;
      end
    end

    ai = 0.0;
    for k = 2 : n
      ai = ai + 1.0;
      p2 = ( ( ai + ai + 1.0 ) * x * p1 - ai * p0 ) / ( ai + 1.0 );
      pd2 = ( ( ai + ai + 1.0 ) * ( p1 + x * pd1 ) - ai * pd0 ) / ( ai + 1.0 );
      p0 = p1;
      p1 = p2;
      pd0 = pd1;
      pd1 = pd2;
    end
%
%  Newton correction.
%
    delta = p2 / pd2;
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
    fprintf ( 1, 'ABWE2 - Fatal error!\n' );
    fprintf ( 1, '  Iteration limit reached.\n' );
    fprintf ( 1, '  Last DELTA was %e\n', delta );
    error ( 'ABWE2 - Fatal error!' );
  end
%
%  Computation of the weight.
%
  an = n;

  w2 = 2.0 / ( an * pd2 * p0 );

  p1 = 0.0;
  p2 = b(m+1);
  yy = 4.0 * x * x - 2.0;
  for k = 1 : m
    i = m - k + 1;
    p0 = p1;
    p1 = p2;
    p2 = yy * p1 - p0 + b(i);
  end

  if ( even )
    w1 = w2 + coef2 / ( pd2 * x * ( p2 - p1 ) );
  else
    w1 = w2 + 2.0 * coef2 / ( pd2 * ( p2 - p0 ) );
  end

  return
end
