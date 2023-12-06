#ifndef _FORMULA_TRIANGLES_DEFORM_H_
#define _FORMULA_TRIANGLES_DEFORM_H_

#include "MMVII_TplSymbImage.h"
#include "MMVII_util_tpl.h"

/**
    \brief  class to generate code for triangle transformation by minimization
**/

using namespace NS_SymbolicDerivative;


namespace MMVII
{

class cTriangleDeformation
{
  public :
    cTriangleDeformation() 
    {
    }

    static const std::vector<std::string> VNamesUnknowns() {return {"GeomTrX","GeomTrY"};}
    static const std::vector<std::string> VNamesObs()      { return {"D"}; }

    std::string FormulaName() const { return "TriangleDeformation";}

    template <typename tUk,typename tObs> 
             static std::vector<tUk> formula
                  (
                      const std::vector<tUk> & aVUk,
                      const std::vector<tObs> & aVObs
                  ) // const
    {
          cPtxd<tUk,2>  p1 = VtoP2(aVUk,0);
          cPtxd<tUk,2>  p2 = VtoP2(aVUk,2);
          cPtxd<tUk,2>  v  = p2-p1;
          const auto & ObsDist  = aVObs[0];  
	  const auto aCst1 = CreateCste(1.0, p1.x());  // create a symbolic formula for constant 1


          return { Norm2(v)/ObsDist - aCst1 } ;
     }
};


}; // MMVII

#endif // _FORMULA_TRIANGLES_DEFORM_H_