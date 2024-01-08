#ifndef _FORMULA_TRIANGLES_DEFORM_H_
#define _FORMULA_TRIANGLES_DEFORM_H_

#include "MMVII_TplSymbTriangle.h"
// #include "MMVII_TplSymbImage.h"
#include "MMVII_util_tpl.h"

#include "SymbDer/SymbolicDerivatives.h"
#include "SymbDer/SymbDer_MACRO.h"

/**
    \file   Formulas_TrianglesDeform.h
    \brief  class to generate code for triangle transformation by minimization
**/

using namespace NS_SymbolicDerivative;

namespace MMVII
{

  class cTriangleDeformation
  {
  public:
    cTriangleDeformation()
    {
    }

    static const std::vector<std::string> VNamesUnknowns() { return {"GeomTrXPointA", "GeomTrYPointA", "GeomTrXPointB", "GeomTrYPointB", "GeomTrXPointC", "GeomTrYPointC"}; }
    static const std::vector<std::string> VNamesObs() 
    {
      return Append
                  (
                    std::vector<std::string>{"PixelCoordinatesX", "PixelCoordinatesY", "AlphaCoordPixel", "BetaCoordPixel", "GammaCoordPixel", "IntensityImPre"},
                    FormalBilinIm2D_NameObs("T")  // 6 obs for bilinear interpol of Im
                     // std::vector<std::string>{"xMod","yMod","ValueMod"} // x,y of point, value of modele
                  );
    }

    std::string FormulaName() const { return "TriangleDeformation"; }

    template <typename tUk, typename tObs>
    static std::vector<tUk> formula(
        const std::vector<tUk> &aVUk,
        const std::vector<tObs> &aVObs) // const
    {
      size_t IndTri = 0;
      // size_t IndX = IndTri + TriangleDisplacement_NbObs;

      // extract observation on model

      const auto &aXCoordinates = aVObs[IndTri];
      const auto &aYCoordinates = aVObs[IndTri + 1];
      const auto &aAlphaCoordinate = aVObs[IndTri + 2];
      const auto &aBetaCoordinate = aVObs[IndTri + 3];
      const auto &aGammaCoordinate = aVObs[IndTri + 4]; // change here
      const auto &aIntensityImPre = aVObs[IndTri + 5];

      // const auto &aComputedXCoordinate = aVObs[IndTri + 6];
      // const auto &aComputedYCoordinate = aVObs[IndTri + 7];  // Only these variables are needed for application of bilinear formula.

      // extract unknowns
      // const auto &aRadSc = aVUk[0];
      // const auto &aRadTr = aVUk[0];
      // const auto &aGeomScale = aVUk[1];

      const auto &aGeomTrXPointA = aVUk[0]; // change here
      const auto &aGeomTrYPointA = aVUk[1];
      const auto &aGeomTrXPointB = aVUk[2];
      const auto &aGeomTrYPointB = aVUk[3];
      const auto &aGeomTrXPointC = aVUk[4];
      const auto &aGeomTrYPointC = aVUk[5];

      auto aXTri = aXCoordinates + aAlphaCoordinate * aGeomTrXPointA + aBetaCoordinate * aGeomTrXPointB + aGammaCoordinate * aGeomTrXPointC; // change
      auto aYTri = aYCoordinates + aAlphaCoordinate * aGeomTrYPointA + aBetaCoordinate * aGeomTrYPointB + aGammaCoordinate * aGeomTrYPointC; // change

      // compute formula of bilinear interpolation
      auto aEstimatedValueTri = FormalBilinTri_Formula(aVObs, TriangleDisplacement_NbObs, aXTri, aYTri);

      // residual is simply the difference between both values
      return {aIntensityImPre - aEstimatedValueTri};
    }
  };

}; // MMVII

#endif // _FORMULA_TRIANGLES_DEFORM_H_