// Copyright Eric Chauvin 2022



// This is licensed under the GNU General
// Public License (GPL).  It is the
// same license that Linux has.
// https://www.gnu.org/licenses/gpl-3.0.html



#pragma once



#include "../CppBase/BasicTypes.h"


// This class is mainly for things from
// #include <cmath>.  Which is included in the
// .cpp file.


class MathC
  {
  private:

  public:
  static Float64 roundF64( Float64 x );
  static Int64 roundI64( Float64 x );
  static Int32 roundI32( Float64 x );
  static Float64 sqroot( Float64 x );


  };
