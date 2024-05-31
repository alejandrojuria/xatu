# include "xatu/asa239.hpp"

namespace asa239 {

//****************************************************************************80

double alnorm ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    alnorm() computes the cumulative density of the standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the MIT license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
  constexpr double a1 = 5.75885480458;
  constexpr double a2 = 2.62433121679;
  constexpr double a3 = 5.92885724438;
  constexpr double b1 = -29.8213557807;
  constexpr double b2 = 48.6959930692;
  constexpr double c1 = -0.000000038052;
  constexpr double c2 = 0.000398064794;
  constexpr double c3 = -0.151679116635;
  constexpr double c4 = 4.8385912808;
  constexpr double c5 = 0.742380924027;
  constexpr double c6 = 3.99019417011;
  constexpr double con = 1.28;
  constexpr double d1 = 1.00000615302;
  constexpr double d2 = 1.98615381364;
  constexpr double d3 = 5.29330324926;
  constexpr double d4 = -15.1508972451;
  constexpr double d5 = 30.789933034;
  constexpr double ltone = 7.0;
  constexpr double p = 0.398942280444;
  constexpr double q = 0.39990348504;
  constexpr double r = 0.398942280385;
  bool up = false;
  constexpr double utzero = 18.66;
  double value;
  double y;
  double z;

  z = x;

  if ( z < 0.0 )
  {
    up = true;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
      / ( y + a1 + b1 
      / ( y + a2 + b2 
      / ( y + a3 ))));
  }
  else
  {
    value = r * std::exp ( - y ) 
      / ( z + c1 + d1 
      / ( z + c2 + d2 
      / ( z + c3 + d3 
      / ( z + c4 + d4 
      / ( z + c5 + d5 
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}

//****************************************************************************80

double gammad ( double x, double p )

//****************************************************************************80
//
//  Purpose:
//
//    gammad() computes the Incomplete Gamma Integral
//
//  Licensing:
//
//    This code is distributed under the MIT license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by B Shea.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    B Shea,
//    Algorithm AS 239:
//    Chi-squared and Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete 
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int IFAULT, error flag.
//    0, no error.
//    1, X < 0 or P <= 0.
//
//    Output, double GAMMAD, the value of the incomplete 
//    Gamma integral.
//
{
  double a;
  double an;
  double arg;
  double b;
  double c;
  constexpr double elimit = - 88.0;
  constexpr double oflo = 1.0E+37;
  constexpr double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  constexpr double tol = 1.0E-14;
  // bool upper;
  double value;
  constexpr double xbig = 1.0E+08;

  value = 0.0;
//
//  Check the input.
//
//   if ( x < 0.0 )
//   {
//     *ifault = 1;
//     return value;
//   }

//   if ( p <= 0.0 )
//   {
//     *ifault = 1;
//     return value;
//   }

//   *ifault = 0;

  if ( x == 0.0 )
  {
    value = 0.0;
    return value;
  }
//
//  If P is large, use a normal approximation.
//
  if ( plimit < p )
  {
    pn1 = 3.0 * std::sqrt ( p ) * ( std::pow ( x / p, 1.0 / 3.0 ) 
    + 1.0 / ( 9.0 * p ) - 1.0 );

    // upper = false;
    // value = alnorm ( pn1, upper );
    value = alnorm ( pn1 );
    return value;
  }
//
//  If X is large set value = 1.
//
  if ( xbig < x )
  {
    value = 1.0;
    return value;
  }
//
//  Use Pearson's series expansion.
//
  if ( x <= 1.0 || x < p )
  {
    arg = p * std::log ( x ) - x - std::lgamma ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + std::log ( value );

    if ( elimit <= arg )
    {
      value = std::exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
//
//  Use a continued fraction expansion.
//
  else 
  {
    arg = p * std::log ( x ) - x - std::lgamma ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if ( pn6 != 0.0 )
      {
        rn = pn5 / pn6;

        if ( std::fabs ( value - rn ) <= r8_min ( tol, tol * rn ) )
        {
          break;
        }
        value = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
//
//  Re-scale terms in continued fraction if terms are large.
//
      if ( oflo <= std::abs ( pn5 ) )
      {
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }

    arg = arg + std::log ( value );

    if ( elimit <= arg )
    {
      value = 1.0 - std::exp ( arg );
    }
    else
    {
      value = 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    r8_min() returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
    return (y < x)? y : x;
//   double value;

//   if ( y < x )
//   {
//     value = y;
//   } 
//   else
//   {
//     value = x;
//   }
//   return value;
}


}