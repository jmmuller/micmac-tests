// File Automatically generated by eLiSe
#include "StdAfx.h"
#include "cREgDistDx_Polyn0.h"


cREgDistDx_Polyn0::cREgDistDx_Polyn0():
    cElCompiledFonc(2)
{
   AddIntRef (cIncIntervale("Intr",0,3));
   Close(false);
}



void cREgDistDx_Polyn0::ComputeVal()
{

  mVal[0] = mLocRegDistu1_x - mLocRegDistu2_x - mLocRegDistu3_x;

  mVal[1] = mLocRegDistu1_y - mLocRegDistu2_y - mLocRegDistu3_y;

}


void cREgDistDx_Polyn0::ComputeValDeriv()
{

  mVal[0] = mLocRegDistu1_x - mLocRegDistu2_x - mLocRegDistu3_x;

  mCompDer[0][0] = 0;
  mCompDer[0][1] = 0;
  mCompDer[0][2] = 0;
  mVal[1] = mLocRegDistu1_y - mLocRegDistu2_y - mLocRegDistu3_y;

  mCompDer[1][0] = 0;
  mCompDer[1][1] = 0;
  mCompDer[1][2] = 0;
}


void cREgDistDx_Polyn0::ComputeValDerivHessian()
{
  ELISE_ASSERT(false,"Foncteur cREgDistDx_Polyn0 Has no Der Sec");
}

void cREgDistDx_Polyn0::SetRegDistu1_x(double aVal){ mLocRegDistu1_x = aVal;}
void cREgDistDx_Polyn0::SetRegDistu1_y(double aVal){ mLocRegDistu1_y = aVal;}
void cREgDistDx_Polyn0::SetRegDistu2_x(double aVal){ mLocRegDistu2_x = aVal;}
void cREgDistDx_Polyn0::SetRegDistu2_y(double aVal){ mLocRegDistu2_y = aVal;}
void cREgDistDx_Polyn0::SetRegDistu3_x(double aVal){ mLocRegDistu3_x = aVal;}
void cREgDistDx_Polyn0::SetRegDistu3_y(double aVal){ mLocRegDistu3_y = aVal;}



double * cREgDistDx_Polyn0::AdrVarLocFromString(const std::string & aName)
{
   if (aName == "RegDistu1_x") return & mLocRegDistu1_x;
   if (aName == "RegDistu1_y") return & mLocRegDistu1_y;
   if (aName == "RegDistu2_x") return & mLocRegDistu2_x;
   if (aName == "RegDistu2_y") return & mLocRegDistu2_y;
   if (aName == "RegDistu3_x") return & mLocRegDistu3_x;
   if (aName == "RegDistu3_y") return & mLocRegDistu3_y;
   return 0;
}


cElCompiledFonc::cAutoAddEntry cREgDistDx_Polyn0::mTheAuto("cREgDistDx_Polyn0",cREgDistDx_Polyn0::Alloc);


cElCompiledFonc *  cREgDistDx_Polyn0::Alloc()
{  return new cREgDistDx_Polyn0();
}


