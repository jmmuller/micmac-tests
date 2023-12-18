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
      size_t IndX = TriangleDisplacement_NbObs + IndTri;

      // extract observation on model
      const auto &aXCoordinates = aVObs[IndX];
      const auto &aYCoordinates = aVObs[IndX + 1];
      const auto &aAlphaCoordinates = aVObs[IndX + 2];
      const auto &aBetaCoordinates = aVObs[IndX + 3];
      const auto &aGammaCoordinates = aVObs[IndX + 4];
      const auto &aIntensityImPre = aVObs[IndX + 5];

      // extract unknowns
      // const auto &aRadSc = aVUk[0];
      // const auto &aRadTr = aVUk[0];
      // const auto &aGeomScale = aVUk[1];
      const auto &aGeomTrXPointA = aVUk[0];
      const auto &aGeomTrYPointA = aVUk[1];
      const auto &aGeomTrXPointB = aVUk[2];
      const auto &aGeomTrYPointB = aVUk[3];
      const auto &aGeomTrXPointC = aVUk[4];
      const auto &aGeomTrYPointC = aVUk[5];

      // auto xIm = aGeomTrx + aGeomScale * xModele;
      // auto yIm = aGeomTry + aGeomScale * yModele;
      auto aXTri = aXCoordinates + aAlphaCoordinates * aGeomTrXPointA + aBetaCoordinates * aGeomTrXPointB + aGammaCoordinates * aGeomTrXPointC;
      auto aYTri = aYCoordinates + aAlphaCoordinates * aGeomTrYPointA + aBetaCoordinates * aGeomTrYPointB + aGammaCoordinates * aGeomTrYPointC;

      // compute formula of bilinear interpolation
      // auto aValueTri = Apply_TriangleMeshDisplacement_Formula(aVObs, IndTri, xTri, yTri);
      auto aEstimatedValueTri = FormalBilinTri_Formula(aVObs, IndTri, aXTri, aYTri);
      // take into account radiometric transform
      // auto aValueModele = aRadTr + aRadSc * vModelInit;
      // auto aValueModele = vModelInit;

      // residual is simply the difference between both values
      return {aEstimatedValueTri - aIntensityImPre};
    }
  };

}; // MMVII

#endif // _FORMULA_TRIANGLES_DEFORM_H_