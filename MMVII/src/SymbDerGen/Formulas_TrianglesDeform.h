#ifndef _FORMULA_TRIANGLES_DEFORM_H_
#define _FORMULA_TRIANGLES_DEFORM_H_

#include "MMVII_TplSymbTriangle.h"
#include "MMVII_util_tpl.h"

#include "SymbDer/SymbolicDerivatives.h"
#include "SymbDer/SymbDer_MACRO.h"

/**
    \file   Formulas_TrianglesDeform.h
    \brief  class to generate code for triangle deformation by minimization
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

    static const std::vector<std::string> VNamesUnknowns() {return Append(std::vector<std::string>{"GeomTrXPointA", "GeomTrYPointA"}, 
                                                                          std::vector<std::string>{"GeomTrXPointB", "GeomTrYPointB"},
                                                                          std::vector<std::string>{"GeomTrXPointC", "GeomTrYPointC"});}
    static const std::vector<std::string> VNamesObs() 
    {
      return Append
                  (
                    std::vector<std::string>{"PixelCoordinatesX", "PixelCoordinatesY", "AlphaCoordPixel", "BetaCoordPixel", "GammaCoordPixel", "IntensityImPre"},
                    FormalBilinIm2D_NameObs("T") // 6 obs for bilinear interpol of Im
                  );
    }

    std::string FormulaName() const {return "TriangleDeformation";}

    template <typename tUk, typename tObs>
    static std::vector<tUk> formula(
        const std::vector<tUk> &aVUnk,
        const std::vector<tObs> &aVObs)
    {
      // size_t IndTri = 0;
      // size_t IndX = IndTri + TriangleDisplacement_NbObs;

      // extract observation on model
      const auto &aXCoordinate = aVObs[0];
      const auto &aYCoordinate = aVObs[1];
      const auto &aAlphaCoordinate = aVObs[2];
      const auto &aBetaCoordinate = aVObs[3];
      const auto &aGammaCoordinate = aVObs[4];
      const auto &aIntensityImPre = aVObs[5];

      // extract unknowns
      // const auto &aRadSc = aVUk[0];
      // const auto &aRadTr = aVUnk[0];
      // const auto &aGeomScale = aVUnk[1];

      const auto &aGeomTrXPointA = aVUnk[0];
      const auto &aGeomTrYPointA = aVUnk[1];
      const auto &aGeomTrXPointB = aVUnk[2];
      const auto &aGeomTrYPointB = aVUnk[3];
      const auto &aGeomTrXPointC = aVUnk[4];
      const auto &aGeomTrYPointC = aVUnk[5];

      auto aXTri = aXCoordinate + aAlphaCoordinate * aGeomTrXPointA + aBetaCoordinate * aGeomTrXPointB + aGammaCoordinate * aGeomTrXPointC;
      auto aYTri = aYCoordinate + aAlphaCoordinate * aGeomTrYPointA + aBetaCoordinate * aGeomTrYPointB + aGammaCoordinate * aGeomTrYPointC;

      // compute formula of bilinear interpolation
      auto aEstimatedValueTri = FormalBilinTri_Formula(aVObs, TriangleDisplacement_NbObs, aXTri, aYTri);

      // auto newVal = mul * aEstimatedValueTri + add

      // residual is simply the difference between values in before image and estimated value in new image.
      return {aIntensityImPre - aEstimatedValueTri};
    }
  };

}; // MMVII

#endif // _FORMULA_TRIANGLES_DEFORM_H_