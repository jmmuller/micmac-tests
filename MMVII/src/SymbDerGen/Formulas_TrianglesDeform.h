#ifndef _FORMULA_TRIANGLES_DEFORM_H_
#define _FORMULA_TRIANGLES_DEFORM_H_

#include "MMVII_TplSymbTriangle.h"
// #include "MMVII_util_tpl.h"

#include "SymbDer/SymbolicDerivatives.h"
#include "SymbDer/SymbDer_MACRO.h"

/**
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

    static const std::vector<std::string> VNamesUnknowns() { return {"GeomTrX", "GeomTrY"}; }
    static const std::vector<std::string> VNamesObs() { return {"PixelCoordinatesX", "PixelCoordinatesY", "IntensitePreIm", "IntensitePostIm"}; }

    std::string FormulaName() const { return "TriangleDeformation"; }

    template <typename tUk, typename tObs>
    static std::vector<tUk> formula(
        const std::vector<tUk> &aVUk,
        const std::vector<tObs> &aVObs) // const
    {
      size_t IndBilin = 0;
      size_t IndX = FormalBilinTriangle_NbObs + IndBilin;

      // extract observation on model
      const auto &xModele = aVObs[IndX];
      const auto &yModele = aVObs[IndX + 1];
      const auto &vModelInit = aVObs[IndX + 2];

      // extract unknowns
      // const auto &aRadSc = aVUk[0];
      // const auto &aRadTr = aVUk[0];
      // const auto &aGeomScale = aVUk[1];
      const auto &aGeomTrx = aVUk[0];
      const auto &aGeomTry = aVUk[1];

      // compute pixel homologous to model in image
      // auto xIm = aGeomTrx + aGeomScale * xModele;
      // auto yIm = aGeomTry + aGeomScale * yModele;
      auto xIm = aGeomTrx + xModele;
      auto yIm = aGeomTry + yModele;

      // compute formula of bilinear interpolation
      auto aValueTri = Apply_TriangleMeshDisplacement_Formula(aVObs, IndBilin, xIm, yIm);
      // take into account radiometric transform
      // auto aValueModele = aRadTr + aRadSc * vModelInit;
      auto aValueModele = vModelInit;

      // residual is simply the difference between both values
      return {aValueModele - aValueTri};
    }
  };

}; // MMVII

#endif // _FORMULA_TRIANGLES_DEFORM_H_