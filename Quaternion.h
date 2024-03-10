// Copyright Eric Chauvin 2022.


/*
This is only partially translated to C++ from
my C# code.

// This is licensed under the GNU General
// Public License (GPL).  It is the
// same license that Linux has.
// https://www.gnu.org/licenses/gpl-3.0.html


#pragma once



// There is a good book about the history of
// vectors and quaternions called A History of
// Vector Analysis, by Michael Crowe.  It
// describes how modern vector analysis evolved,
// and how Hamilton's ideas about quaternions
// evolved.


#include "BasicTypes.h"
#include "MathC.h"


class Quaternion
  {
  // a + bi + cj + dk
  // w + bx + cy + dz

  private:
  Float64 X;
  Float64 Y;
  Float64 Z;
  Float64 W;

  inline void setToOne( void )
    {
    X = 0;
    Y = 0;
    Z = 0;
    W = 1;
    }



  inline void setToZero( void )
    {
    X = 0;
    Y = 0;
    Z = 0;
    W = 0;
    }



  inline void negate( void  )
    {
    X = -X;
    Y = -Y;
    Z = -Z;
    W = -W;
    }




  inline void add( Quaternion& toAdd )
    {
    X += toAdd.X;
    Y += toAdd.Y;
    Z += toAdd.Z;
    W += toAdd.W;
    }



  inline Float64 normSquared( Quaternion& in )
    {
    Float64 result = (in.X * in.X) +
                     (in.Y * in.Y) +
                     (in.Z * in.Z) +
                     (in.W * in.W);

    return result;
    }



  inline Float64 norm( Quaternion& in )
    {
    Float64 nSquared = normSquared( in );
    return MathC::sqroot( nSquared );
    }



  static inline void normalize( Quaternion& in )
    {
    Float64 length = norm( In );
    if( length < 0.000000000000000000000000000001d )
      {
      in.setToZero();
      return;
      }

    Float64 inverse = 1.0d / length;

    in.X = in.X * inverse;
    in.Y = in.Y * inverse;
    in.Z = in.Z * inverse;
    in.W = in.W * inverse;
    }


=========

  internal static bool DoubleIsAlmostEqual( double A, double B, double SmallNumber )
    {
    // How small can this be?

    if( A + SmallNumber < B )
      return false;

    if( A - SmallNumber > B )
      return false;

    return true;
    }



  internal static bool IsAlmostEqual( QuaternionRec Left, QuaternionRec Right, double SmallNumber )
    {
    if( !DoubleIsAlmostEqual( Left.X, Right.X, SmallNumber ))
      return false;

    if( !DoubleIsAlmostEqual( Left.Y, Right.Y, SmallNumber ))
      return false;

    if( !DoubleIsAlmostEqual( Left.Z, Right.Z, SmallNumber ))
      return false;

    if( !DoubleIsAlmostEqual( Left.W, Right.W, SmallNumber ))
      return false;

    return true;
    }



  internal static QuaternionRec Conjugate( QuaternionRec In )
    {
    QuaternionRec Result;
    Result.X = -In.X;
    Result.Y = -In.Y;
    Result.Z = -In.Z;
    Result.W = In.W;
    return Result;
    }



  internal static QuaternionRec Inverse( QuaternionRec In )
    {
    // QX = 1, so X is the multiplicative inverse
    // of Q.  So X = 1 / Q, or Q^(-1).

    double NSquared = NormSquared( In );
    if( NSquared < 0.0000000000001 )
      return SetZero();
      // throw( new Exception( "NSquared is too small in QuaternionEC.Inverse()." ));

    double InverseNS = 1.0d / NSquared;

    // The negative parts are to make it the
    // conjugate.  So the Result is the conjugate
    // divided by the norm squared.
    QuaternionRec Result;
    Result.X = -In.X * InverseNS;
    Result.Y = -In.Y * InverseNS;
    Result.Z = -In.Z * InverseNS;
    Result.W = In.W * InverseNS;
    return Result;
    }


// Partial Differential equations were around for
// about a couple of centuries before quaternions.
// The original meaining of the phrase Cross
// Product is with differentials like dxdy as
// opposed to dx squared.  And it has to do with
// the Metric ds squared.
// (a + b + c)(a + b + c)
// aa + ab + ac + and so on...
// aa is not a cross product.
// ab  and ac are cross products.


  internal static QuaternionRec CrossProduct(
                                 QuaternionRec Left,
                                 QuaternionRec Right )
    {
    // These are cross products:
    // i x j = k
    // j x k = i
    // k x i = j

    QuaternionRec Result;
    // W is not used.
    Result.W = 0;

    Result.X = (Left.Y * Right.Z) -
               (Left.Z * Right.Y);

    Result.Y = (Left.Z * Right.X) -
               (Left.X * Right.Z);

    Result.Z = (Left.X * Right.Y) -
               (Left.Y * Right.X);

    return Result;
    }




  ///////////////////////////
  // Notes on Multiplication:

  // It is a right-handed coordinate system.  Positive
  // Z values go toward the viewer.  X goes to the
  // right, Y goes up.

  // ij = k    ji = -k
  // jk = i    kj = -i
  // ki = j    ik = -j

  // xy = z    yx = -z
  // yz = x    zy = -x
  // zx = y    xz = -y

  //   i^2 = j^2 = k^2 = ijk = -1

  // With two regular complex numbers you do:
  //      (a + bi)(c + di) =
  //      ac + adi + bic + bidi

  // a1 is like a with subscript 1.

  // Multiply two quaternions:
  // Left times Right.
  // (a1x + b1y + c1z + w1)(a2x + b2y + c2z + w2)

  // Distributive Property:
  // a1x(a2x + b2y + c2z + w2) +
  // b1y(a2x + b2y + c2z + w2) +
  // c1z(a2x + b2y + c2z + w2) +
  // w1(a2x + b2y + c2z + w2)

  // Distributive Property again:
  // a1xa2x + a1xb2y + a1xc2z + a1xw2 +
  // b1ya2x + b1yb2y + b1yc2z + b1yw2 +
  // c1za2x + c1zb2y + c1zc2z + c1zw2 +
  // w1a2x + w1b2y + w1c2z + w1w2

  // a1a2xx + a1b2xy + a1c2xz + a1w2x +
  // b1a2yx + b1b2yy + b1c2yz + b1w2y +
  // c1a2zx + c1b2zy + c1c2zz + c1w2z +
  // w1a2x + w1b2y + w1c2z + w1w2

  // xy = z    yx = -z
  // yz = x    zy = -x
  // zx = y    xz = -y
  // a1a2-1 + a1b2z + a1c2-y + a1w2x +
  // b1a2-z + b1b2-1 + b1c2x + b1w2y +
  // c1a2y + c1b2-x + c1c2-1 + c1w2z +
  // w1a2x + w1b2y + w1c2z + w1w2

  // Rearrange it so the components are together:
  // a1w2x + w1a2x + b1c2x + c1b2-x +
  // a1c2-y + b1w2y + c1a2y + w1b2y +
  // a1b2z + b1a2-z + c1w2z + w1c2z +
  // -a1a2 + -b1b2 + -c1c2 + w1w2

  // a1w2x + w1a2x + b1c2x + -c1b2x +
  // -a1c2y + b1w2y + c1a2y + w1b2y +
  // a1b2z + -b1a2z + c1w2z + w1c2z +
  // -a1a2 + -b1b2 + -c1c2 + w1w2

  // x(a1w2 + w1a2 + b1c2 + -c1b2) +
  // y(-a1c2 + b1w2 + c1a2 + w1b2) +
  // z(a1b2 + -b1a2 + c1w2 + w1c2) +
  // The W parts:
  // -a1a2 + -b1b2 + -c1c2 + w1w2
  ///////////////////////////




  internal static QuaternionRec Multiply(
                                     QuaternionRec L,
                                     QuaternionRec R )
    {
    // Result.X = a1w2 + w1a2 + b1c2 + -c1b2;
    // Result.Y = -a1c2 + b1w2 + c1a2 + w1b2;
    // Result.Z = a1b2 + -b1a2 + c1w2 + w1c2;
    // Result.W = -a1a2 + -b1b2 + -c1c2 + w1w2;

    QuaternionRec Result;

    // The vector Cross Product part:
    Result.X =  (L.X * R.W) +  (L.W * R.X) +  (L.Y * R.Z) + (-L.Z * R.Y);
    Result.Y = (-L.X * R.Z) +  (L.Y * R.W) +  (L.Z * R.X) +  (L.W * R.Y);
    Result.Z =  (L.X * R.Y) + (-L.Y * R.X) +  (L.Z * R.W) +  (L.W * R.Z);

    // Almost the same as the vector Dot Product.
    Result.W = (-L.X * R.X) + (-L.Y * R.Y) + (-L.Z * R.Z) +  (L.W * R.W);
    return Result;
    }



  internal static QuaternionRec MultiplyWithLeftVector3(
                                 Vector3.Vector L,
                                 QuaternionRec R )
    {
    QuaternionRec Result;

    // Result.X =  (L.X * R.W) +  (0 * R.X) +  (L.Y * R.Z) + (-L.Z * R.Y);
    // Result.Y = (-L.X * R.Z) +  (L.Y * R.W) +  (L.Z * R.X) +  (0 * R.Y);
    // Result.Z =  (L.X * R.Y) + (-L.Y * R.X) +  (L.Z * R.W) +  (0 * R.Z);
    // Result.W = (-L.X * R.X) + (-L.Y * R.Y) + (-L.Z * R.Z) +  (0 * R.W);

    Result.X =  (L.X * R.W) +  (L.Y * R.Z) + (-L.Z * R.Y);
    Result.Y = (-L.X * R.Z) +  (L.Y * R.W) +  (L.Z * R.X);
    Result.Z =  (L.X * R.Y) + (-L.Y * R.X) +  (L.Z * R.W);
    Result.W = (-L.X * R.X) + (-L.Y * R.Y) + (-L.Z * R.Z);
    return Result;
    }



  internal static Vector3.Vector MultiplyWithResultVector3(
                                  QuaternionRec L,
                                  QuaternionRec R )
    {
    Vector3.Vector Result;
    Result.X =  (L.X * R.W) +  (L.W * R.X) +  (L.Y * R.Z) + (-L.Z * R.Y);
    Result.Y = (-L.X * R.Z) +  (L.Y * R.W) +  (L.Z * R.X) +  (L.W * R.Y);
    Result.Z =  (L.X * R.Y) + (-L.Y * R.X) +  (L.Z * R.W) +  (L.W * R.Z);

    // It doesn't need this calculation:
    // Result.W = (-L.X * R.X) + (-L.Y * R.Y) + (-L.Z * R.Z) +  (L.W * R.W);

    return Result;
    }

Since it's taking the sin and cos of the angle
... what about Earth and 24 hour time?

  internal static QuaternionRec SetAsRotation(
                                 QuaternionRec Axis,
                                 double Angle )
    {
    // Make sure it's a unit quaternion.
    Axis.W = 0;
    Axis = Normalize( Axis );

    // If Angle was Pi / 2 then this would be
    // Pi / 4.
    double HalfAngle = Angle * 0.5d;
    double SineHalfAngle = Math.Sin( HalfAngle );
    double CosineHalfAngle = Math.Cos( HalfAngle );

    QuaternionRec Result;
    Result.X = Axis.X * SineHalfAngle;
    Result.Y = Axis.Y * SineHalfAngle;
    Result.Z = Axis.Z * SineHalfAngle;
    Result.W = CosineHalfAngle;
    return Result;
    }



  internal static QuaternionRec Rotate(
                       QuaternionRec RotationQ,
                       QuaternionRec InverseRotationQ,
                       QuaternionRec StartPoint )
    {
    QuaternionRec MiddlePoint = Multiply( StartPoint, InverseRotationQ );
    QuaternionRec Result = Multiply( RotationQ, MiddlePoint );
    return Result;
    }



  internal static Vector3.Vector RotateVector3(
                     QuaternionRec RotationQ,
                     QuaternionRec InverseRotationQ,
                     Vector3.Vector StartPoint )
    {
    // A quaternion rotation is clockwise around a vector
    // if you are looking down the vector from the origin point.
    // Like an archer with an arrow in the bow, you are sighting
    // down the arrow.
    // Compare this with representing a rotation or a moment
    // of inertia, or a torque, by using a vector cross product
    // in a right-handed coordinate system.  It goes in the
    // right direction like it should.  If the Z axis is pointing
    // straight toward you then it is the opposite
    // point of view from what an archer would see when he
    // is sighting down an arrow in his bow.  That opposite
    // point of view is a counter-clockwise rotation.


    QuaternionRec MiddlePoint = MultiplyWithLeftVector3(
                       StartPoint, InverseRotationQ );
    Vector3.Vector Result = MultiplyWithResultVector3(
                              RotationQ, MiddlePoint );
    return Result;
    }




  internal static Vector3.Vector RotationWithSetupDegrees(
                                  double AngleDegrees,
                                  QuaternionRec Axis,
                                  Vector3.Vector InVector )
    {
    double Angle = NumbersEC.DegreesToRadians( AngleDegrees );

    QuaternionRec RotationQ = SetAsRotation( Axis, Angle );
    QuaternionRec InverseRotationQ = Inverse( RotationQ );

    Vector3.Vector ResultPoint = RotateVector3(
                     RotationQ,
                     InverseRotationQ,
                     InVector );

    return ResultPoint;
    }



  }
}

*/
