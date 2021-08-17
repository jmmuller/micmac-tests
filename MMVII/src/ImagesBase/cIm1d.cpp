#include "include/MMVII_all.h"

namespace MMVII
{

/** Class for "tabulating" a 1D fonction in a given interval.
    The value are memorized in a tab, and returned when needed
*/

/* ========================== */
/*       cTabulFonc1D         */
/* ========================== */

double cTabulFonc1D::F(double aX) const 
{
   int aK = ToIntCoord(aX);
   if (aK<0) return mValXMin;
   if (aK>mNbStep) return mValXMax;
   return mDIm->GetV(aK);
}

double cTabulFonc1D::ToRealCoord(int   aI) const
{
   return mXMin + aI * mStep;
}

int    cTabulFonc1D::ToIntCoord(double aX) const
{
   return round_ni((aX-mXMin)/mStep);
}


cTabulFonc1D::cTabulFonc1D
(
    const cFctrRR & aFctr,
    double aXMin,
    double aXMax,
    int aNbStep
) :
    mXMin    (aXMin),
    mXMax    (aXMax),
    mNbStep  (aNbStep),
    mStep    ((mXMax-mXMin)/mNbStep),
    mValXMin (aFctr.F(mXMin)),
    mValXMax (aFctr.F(mXMax)),
    mIm      (mNbStep+1),
    mDIm     (&mIm.DIm())
{
    for (int aK=0 ; aK<=mNbStep ; aK++)
        mDIm->SetV(aK,aFctr.F(ToRealCoord(aK)));
}

/* ========================== */
/*          cDataIm1D         */
/* ========================== */

template <class Type>  void cDataIm1D<Type>::PostInit()
{
    mRawData1D = tBI::mRawDataLin - X0();
}

template <class Type>  cDataIm1D<Type>::cDataIm1D
                       (
                            const cPt1di & aP0,
                            const cPt1di & aP1,
                            Type * aRawDataLin,
                            eModeInitImage aModeInit
                       ) : 
    cDataTypedIm<Type,1> (aP0,aP1,aRawDataLin,aModeInit)
{
    PostInit();
}

template <class Type>  void cDataIm1D<Type>::Resize ( const cPt1di & aP0, const cPt1di & aP1, eModeInitImage aModeInit)  
{
   cDataTypedIm<Type,1>::Resize(aP0,aP1,aModeInit);
   PostInit();
}

template <class Type>  void cDataIm1D<Type>::Resize ( int aSz, eModeInitImage aModeInit)  
{
    Resize(cPt1di(0),cPt1di(aSz),aModeInit);
}

template <class Type>  cDataIm1D<Type>::~cDataIm1D()
{
}

template <class Type> int     cDataIm1D<Type>::VI_GetV(const cPt1di& aP)  const
{
   return GetV(aP);
}
template <class Type> double  cDataIm1D<Type>::VD_GetV(const cPt1di& aP)  const 
{
   return GetV(aP);
}

template <class Type> void  cDataIm1D<Type>::VI_SetV(const cPt1di& aP,const int & aV)
{
   SetVTrunc(aP,aV);
}
template <class Type> void  cDataIm1D<Type>::VD_SetV(const cPt1di& aP,const double & aV)
{
   SetVTrunc(aP,aV);
}


/* ========================== */
/*          cIm1D             */
/* ========================== */

template <class Type>  cIm1D<Type>::cIm1D(const int & aP0,const int & aP1,Type * aRawDataLin,eModeInitImage aModeInit) :
   mSPtr(new cDataIm1D<Type>(cPt1di(aP0),cPt1di(aP1),aRawDataLin,aModeInit)),
   mPIm (mSPtr.get())
{
}

template <class Type>  cIm1D<Type>::cIm1D(const int & aSz,Type * aRawDataLin,eModeInitImage aModeInit) :
   cIm1D<Type> (0,aSz,aRawDataLin,aModeInit)
{
}

template <class Type>  cIm1D<Type>  cIm1D<Type>::Dup() const
{
   cIm1D<Type> aRes(DIm().P0().x(),DIm().P1().x());
   DIm().DupIn(aRes.DIm());
   return aRes;
}


/* ========================== */
/*        cHistoCumul         */
/* ========================== */

template <class TypeH,class TypeCumul>  cHistoCumul<TypeH,TypeCumul>::cHistoCumul(int aNbVal) :
    mNbVal  (aNbVal),
    mH      (aNbVal,nullptr,eModeInitImage::eMIA_Null),
    mDH     (&mH.DIm()),
    mHC     (aNbVal,nullptr,eModeInitImage::eMIA_Null),
    mDHC    (&mHC.DIm()),
    mHCOk   (false),
    mPopTot (0.0)
{
}

template <class TypeH,class TypeCumul>  cHistoCumul<TypeH,TypeCumul>::cHistoCumul() :
    cHistoCumul<TypeH,TypeCumul>(1)
{
}

template <class TypeH,class TypeCumul> void cHistoCumul<TypeH,TypeCumul>::AddV(const int & aP,const TypeH & aV2Add)
{
    mHCOk = false;
    mDH->AddV(aP,aV2Add);
}

template <class TypeH,class TypeCumul> void cHistoCumul<TypeH,TypeCumul>::MakeCumul()
{
    mHCOk = true;
    mPopTot = 0.0;
    const TypeH * aDataH = mDH->RawDataLin();
    TypeCumul * aDataHC = mDHC->RawDataLin();
    for (int aKV=0 ; aKV<mNbVal; aKV++)
    {
        mPopTot += aDataH[aKV];
        aDataHC[aKV] = mPopTot;
    }
}

template <class TypeH,class TypeCumul> tREAL8 cHistoCumul<TypeH,TypeCumul>::PropCumul(const int & aP) const
{
   MMVII_INTERNAL_ASSERT_tiny(mHCOk,"PropCumul  , Cumul Not Ok");
   return mDHC->GetV(aP) / mPopTot;
}

template <class TypeH,class TypeCumul>  const cDataIm1D<TypeH>& cHistoCumul<TypeH,TypeCumul>::H() const
{
    return *mDH;
}

    
#define INSTANTIATE_HCUMUL(TYPE_H,TYPE_C)\
template  class cHistoCumul<TYPE_H,TYPE_C>;


INSTANTIATE_HCUMUL(tREAL4,tREAL8);
INSTANTIATE_HCUMUL(tINT4,tREAL8);
INSTANTIATE_HCUMUL(tINT8,tINT8);




#define INSTANTIATE_IM1D(Type)\
template  class cIm1D<Type>;\
template  class cDataIm1D<Type>;

INSTANTIATE_IM1D(tINT1)
INSTANTIATE_IM1D(tINT2)
INSTANTIATE_IM1D(tINT4)
INSTANTIATE_IM1D(tINT8)

INSTANTIATE_IM1D(tU_INT1)
INSTANTIATE_IM1D(tU_INT2)
INSTANTIATE_IM1D(tU_INT4)


INSTANTIATE_IM1D(tREAL4)
INSTANTIATE_IM1D(tREAL8)
INSTANTIATE_IM1D(tREAL16)



};
