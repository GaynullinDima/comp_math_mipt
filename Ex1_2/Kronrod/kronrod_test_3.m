function kronrod_test_3 ( )


  exact = 99999999999;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'KRONROD_TEST_3\n' );

%
%  TOL just tells KRONROD how carefully it must compute X, W1 and W2.
%  It is NOT a statement about the accuracy of your integral estimate!
%
  tol = 0.000001;
%
%  Start the process with a 1 point rule.
%
  n = 1;

  while ( 1 )

    [ x, w1, w2 ] = kronrod ( n, tol );

    i1 = w1(n+1) * f ( x(n+1) );
    i2 = w2(n+1) * f ( x(n+1) );

    for i = 1 : n
      i1 = i1 + w1(i) * ( f ( - x(i) ) + f ( x(i) ) );
      i2 = i2 + w2(i) * ( f ( - x(i) ) + f ( x(i) ) );
    end

    if ( abs ( i1 - i2 ) < 0.0001 )
      fprintf ( 1, '\n' );
      fprintf ( 1, '  Error tolerance satisfied with N = %d\n', n );
      fprintf ( 1, '  Coarse integral estimate = %14.6g\n', i1 );
      fprintf ( 1, '  Fine   integral estimate = %14.6g\n', i2 );
      fprintf ( 1, '  Error estimate = %g\n', abs ( i2 - i1 ) );
      fprintf ( 1, '  Actual error =   %g\n', abs ( exact - i2 ) );
      break;
    end

    if ( 25 < n )
      fprintf ( 1, '\n' );
      fprintf ( 1, '  Error tolerance failed even for n = %d\n', n );
      fprintf ( 1, '  Canceling iteration, and accepting bad estimates!\n' );
      fprintf ( 1, '  Coarse integral estimate = %14.6g\n', i1 );
      fprintf ( 1, '  Fine   integral estimate = %14.6g\n', i2 );
      fprintf ( 1, '  Error estimate = %g\n', abs ( i2 - i1 ) );
      fprintf ( 1, '  Actual error =   %g\n', abs ( exact - i2 ) );
      break;
    end

    n = 2 * n + 1;

  end

  return
end


function value = f ( x )
  value = 1.0 /  x ;

  return;
end
