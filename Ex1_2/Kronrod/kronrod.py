#! /usr/bin/env python
#
def kronrad_gauss ( n, m, tol, coef2, even, b, x ) :
  from sys import exit

  if ( x == 0.0 ):
    ka = 1
  else:
    ka = 0

  for iter in range ( 1, 51 ):

    b1 = 0.0
    b2 = b[m+1-1]
    yy = 4.0 * x * x - 2.0
    d1 = 0.0

    if ( even ):
      ai = m + m + 1
      d2 = ai * b[m+1-1]
      dif = 2.0
    else:
      ai = m + 1
      d2 = 0.0
      dif = 1.0

    for k in range ( 1, m + 1 ):
      ai = ai - dif
      i = m - k + 1
      b0 = b1
      b1 = b2
      d0 = d1
      d1 = d2
      b2 = yy * b1 - b0 + b[i-1]
      if ( not even ):
        i = i + 1
      d2 = yy * d1 - d0 + ai * b[i-1]

    if ( even ):
      f = x * ( b2 - b1 )
      fd = d2 + d1
    else:
      f = 0.5 * ( b2 - b0 )
      fd = 4.0 * x * d2
#
#  Newton correction.
#
    delta = f / fd
    x = x - delta

    if ( ka == 1 ):
      break

    if ( abs ( delta ) <= tol ):
      ka = 1
#
#  Catch non-convergence.
#
  if ( ka != 1 ):
    print ( '' )
    print ( 'ABWE1 - Fatal error!' )
    print ( '  Iteration limit reached.' )
    print ( '  Last DELTA was %e' % ( delta ) )
    exit ( 'ABWE1 - Fatal error!' )
#
#  Computation of the weight.
#
  d0 = 1.0
  d1 = x
  ai = 0.0
  for k in range ( 2, n + 1 ):
    ai = ai + 1.0
    d2 = ( ( ai + ai + 1.0 ) * x * d1 - ai * d0 ) / ( ai + 1.0 )
    d0 = d1
    d1 = d2

  w = coef2 / ( fd * d2 )

  return x, w

def gauss ( n, m, tol, coef2, even, b, x ):
  from sys import exit

  if ( x == 0.0 ):
    ka = 1
  else:
    ka = 0
#
#  Iterative process for the computation of a Gaussian abscissa.
#
  for iter in range ( 1, 51 ):

    p0 = 1.0
    p1 = x
    pd0 = 0.0
    pd1 = 1.0
#
#  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
#
    if ( n <= 1 ):
      if ( x != 0.0 ):
        p2 = ( 3.0 * x * x - 1.0 ) / 2.0
        pd2 = 3.0 * x
      else:
        p2 = 3.0 * x
        pd2 = 3.0

    ai = 0.0
    for k in range ( 2, n + 1 ):
      ai = ai + 1.0
      p2 = ( ( ai + ai + 1.0 ) * x * p1 - ai * p0 ) / ( ai + 1.0 )
      pd2 = ( ( ai + ai + 1.0 ) * ( p1 + x * pd1 ) - ai * pd0 ) / ( ai + 1.0 )
      p0 = p1
      p1 = p2
      pd0 = pd1
      pd1 = pd2
#
#  Newton correction.
#
    delta = p2 / pd2
    x = x - delta

    if ( ka == 1 ):
      break

    if ( abs ( delta ) <= tol ):
      ka = 1
#
#  Catch non-convergence.
#
  if ( ka != 1 ):
    print ( '' )
    print ( 'ABWE2 - Fatal error!' )
    print ( '  Iteration limit reached.' )
    print ( '  Last DELTA was %e' % ( delta ) )
    exit ( 'ABWE2 - Fatal error!' )
#
#  Computation of the weight.
#
  an = n

  w2 = 2.0 / ( an * pd2 * p0 )

  p1 = 0.0
  p2 = b[m+1-1]
  yy = 4.0 * x * x - 2.0
  for k in range ( 1, m + 1 ):
    i = m - k + 1
    p0 = p1
    p1 = p2
    p2 = yy * p1 - p0 + b[i-1]

  if ( even ):
    w1 = w2 + coef2 / ( pd2 * x * ( p2 - p1 ) )
  else:
    w1 = w2 + 2.0 * coef2 / ( pd2 * ( p2 - p0 ) )

  return x, w1, w2



def kronrod ( n, tol ):
  import numpy as np

  m = ( ( n + 1 ) // 2 )
  even = ( 2 * m == n )

  b = np.zeros ( m + 1 )
  tau = np.zeros ( m )
  w1 = np.zeros ( n + 1 )
  w2 = np.zeros ( n + 1 )
  x = np.zeros ( n + 1 )

  d = 2.0
  an = 0.0
  for k in range ( 1, n + 1 ):
    an = an + 1.0
    d = d * an / ( an + 0.5 )
#
#  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
#
  tau[1-1] = ( an + 2.0 ) / ( an + an + 3.0 )
  b[m-1] = tau[1-1] - 1.0
  ak = an

  for l in range ( 1, m ):
    ak = ak + 2.0
    tau[l+1-1] = ( ( ak - 1.0 ) * ak \
      - an * ( an + 1.0 ) ) * ( ak + 2.0 ) * tau[l-1] \
      / ( ak * ( ( ak + 3.0 ) * ( ak + 2.0 ) \
      - an * ( an + 1.0 ) ) )
    b[m-l-1] = tau[l+1-1]

    for ll in range ( 1, l + 1 ):
      b[m-l-1] = b[m-l-1] + tau[ll-1] * b[m-l+ll-1]

  b[m+1-1] = 1.0
#
#  Calculation of approximate values for the abscissas.
#
  bb = np.sin ( 0.5 * np.pi / ( an + an + 1.0 ) )
  x1 = np.sqrt ( 1.0 - bb * bb )
  s = 2.0 * bb * x1
  c = np.sqrt ( 1.0 - s * s )
  coef = 1.0 - ( 1.0 - 1.0 / an ) / ( 8.0 * an * an )
  xx = coef * x1
#
#  Coefficient needed for weights.
#
#  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
#
  coef2 = 2.0 / ( 2 * n + 1 )
  for i in range ( 1, n + 1 ):
    coef2 = coef2 * 4.0 * i / ( n + i )
#
#  Calculation of the K-th abscissa (a Kronrod abscissa) and the
#  corresponding weight.
#
  for k in range ( 1, n + 1, 2 ):

    [ xx, w1[k-1] ] = kronrad_gauss ( n, m, tol, coef2, even, b, xx )
    w2[k-1] = 0.0

    x[k-1] = xx
    y = x1
    x1 = y * c - bb * s
    bb = y * s + bb * c

    if ( k == n ):
      xx = 0.0
    else:
      xx = coef * x1
#
#  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
#  corresponding weights.
#
    [ xx, w1[k+1-1], w2[k+1-1] ] = gauss ( n, m, tol, coef2, even, b, xx )

    x[k+1-1] = xx
    y = x1
    x1 = y * c - bb * s
    bb = y * s + bb * c
    xx = coef * x1
#
#  If N is even, we have one more Kronrod abscissa to compute,
#  namely the origin.
#
  if ( even ):
    xx = 0.0
    [ xx, w1[n+1-1] ] = kronrad_gauss ( n, m, tol, coef2, even, b, xx )
    w2[n+1-1] = 0.0
    x[n+1-1] = xx

  return x, w1, w2

def kronrod_adjust ( a, b, n, x, w1, w2 ):
  for i in range ( 0, n + 1 ):

    x[i] = ( ( 1.0 - x[i] ) * a   \
           + ( 1.0 + x[i] ) * b ) \
             / 2.0

    w1[i] = ( ( b - a ) / 2.0 ) * w1[i]
    w2[i] = ( ( b - a ) / 2.0 ) * w2[i]

  return x, w1, w2




def kronrod_test_1 ( ):

  def f ( x ):
    value = 1.0 / ( x * x + 1.005 )

    return value


  import platform

  exact = 1.5643964440690497731

  print ( '' )
  print ( 'KRONROD_TEST_1' )
#
#  TOL just tells KRONROD how carefully it must compute X, W1 and W2.
#  It is NOT a statement about the accuracy of your integral estimate!
#
  tol = 0.000001
#
#  Start the process with a 1 point rule.
#
  n = 1

  while ( 1 ):

    x, w1, w2 = kronrod ( n, tol )

    i1 = w1[n+1-1] * f ( x[n+1-1] )
    i2 = w2[n+1-1] * f ( x[n+1-1] )

    for i in range ( 1, n + 1 ):
      i1 = i1 + w1[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )
      i2 = i2 + w2[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )

    if ( abs ( i1 - i2 ) < 0.0001 ):
      print ( '' )
      print ( '  Error tolerance satisfied with N = %d' % ( n ) )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    if ( 25 < n ):
      print ( '' )
      print ( '  Error tolerance failed even for n = %d' % ( n ) )
      print ( '  Canceling iteration, and accepting bad estimates!' )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    n = 2 * n + 1

  return


def kronrod_test_2 ( ):
  import math

  def f ( x ):
    value = math.sin(x)

    return value


  import platform

  exact = 0.0000000

  print ( '' )
  print ( 'KRONROD_TEST_2' )
#
#  TOL just tells KRONROD how carefully it must compute X, W1 and W2.
#  It is NOT a statement about the accuracy of your integral estimate!
#
  tol = 0.000001
#
#  Start the process with a 1 point rule.
#
  n = 1

  while ( 1 ):

    x, w1, w2 = kronrod ( n, tol )

    i1 = w1[n+1-1] * f ( x[n+1-1] )
    i2 = w2[n+1-1] * f ( x[n+1-1] )

    for i in range ( 1, n + 1 ):
      i1 = i1 + w1[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )
      i2 = i2 + w2[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )

    if ( abs ( i1 - i2 ) < 0.0001 ):
      print ( '' )
      print ( '  Error tolerance satisfied with N = %d' % ( n ) )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    if ( 25 < n ):
      print ( '' )
      print ( '  Error tolerance failed even for n = %d' % ( n ) )
      print ( '  Canceling iteration, and accepting bad estimates!' )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    n = 2 * n + 1

  return


def kronrod_test_3 ( ):
  import math

  def f ( x ):
    value = 1.0 / x

    return value


  import platform

  exact = math.inf

  print ( '' )
  print ( 'KRONROD_TEST_3' )
#
#  TOL just tells KRONROD how carefully it must compute X, W1 and W2.
#  It is NOT a statement about the accuracy of your integral estimate!
#
  tol = 0.000001
#
#  Start the process with a 1 point rule.
#
  n = 1

  while ( 1 ):

    x, w1, w2 = kronrod ( n, tol )

    i1 = w1[n+1-1] * f ( x[n+1-1] )
    i2 = w2[n+1-1] * f ( x[n+1-1] )

    for i in range ( 1, n + 1 ):
      i1 = i1 + w1[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )
      i2 = i2 + w2[i-1] * ( f ( - x[i-1] ) + f ( x[i-1] ) )

    if ( abs ( i1 - i2 ) < 0.0001 ):
      print ( '' )
      print ( '  Error tolerance satisfied with N = %d' % ( n ) )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    if ( 25 < n ):
      print ( '' )
      print ( '  Error tolerance failed even for n = %d' % ( n ) )
      print ( '  Canceling iteration, and accepting bad estimates!' )
      print ( '  Coarse integral estimate = %14.6g' % ( i1 ) )
      print ( '  Fine   integral estimate = %14.6g' % ( i2 ) )
      print ( '  Error estimate = %g' % ( abs ( i2 - i1 ) ) )
      print ( '  Actual error =   %g' % ( abs ( exact - i2 ) ) )
      break

    n = 2 * n + 1

  return


if ( __name__ == '__main__' ):
  kronrod_test_1 ( )
  kronrod_test_2 ( )
  kronrod_test_3 ( )
