// Copyright Eric Chauvin 2022



// This is licensed under the GNU General
// Public License (GPL).  It is the
// same license that Linux has.
// https://www.gnu.org/licenses/gpl-3.0.html



#include "MathC.h"
#include "../CppBase/Casting.h"


// Always include system files at the bottom
// under other include files.
#include <cmath>


// Linux has a long double but not Windows?


Float64 MathC::sqroot( Float64 x )
{
return std::sqrt( x );
}


Float64 MathC::roundF64( Float64 x )
{
return std::round( x );
}

Int64 MathC::roundI64( Float64 x )
{
return std::lround( x );
}


Int32 MathC::roundI32( Float64 x )
{
return Casting::i64ToI32( std::lround( x ));
}



/*
// These have long double versions too.
  acos(float __x)
  asin(float __x)
  atan(float __x)
  atan2(float __y, float __x)
  ceil(float __x)
  cos(float __x)
  cosh(float __x)
  exp(float __x)
  fabs(float __x)
  floor(float __x)
  fmod(float __x, float __y)
  frexp(float __x, int* __exp)
  ldexp(float __x, int __exp)
  log(float __x)
  log10(float __x)
  modf(float __x, float* __iptr)
  pow(float __x, float __y)
  pow(float __x, int __n)
  sin(float __x)
  sinh(float __x)
  tan(float __x)
  tanh(long double __x)

  // Is infinite.
  isinf(double __x)

  // Is not a number.
  isnan(double __x)

*/
