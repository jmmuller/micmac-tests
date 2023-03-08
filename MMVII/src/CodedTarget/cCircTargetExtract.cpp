#include "MMVII_Tpl_Images.h"
#include "MMVII_Linear2DFiltering.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Sensor.h"
#include "MMVII_TplImage_PtsFromValue.h"
#include "MMVII_ImageInfoExtract.h"
#include "CodedTarget.h"

/*   Modularistion
 *   Code extern tel que ellipse
 *   Ellipse => avec centre
 *   Pas de continue
 */

namespace MMVII
{

using namespace cNS_CodedTarget;

class cCCDecode
{
    public :
         cCCDecode(const cExtractedEllipse & anEE,const cDataIm2D<tREAL4> & aDIm,const cFullSpecifTarget &);

	 void Show(const std::string & aPrefix);

	 void ComputePhaseTeta() ;
	 void ComputeCode(bool Show) ;
    private :

	 tREAL8 Dev(int aK1,int aK2) const;
	 tREAL8 Avg(int aK1,int aK2) const;
	 tREAL8 DevOfPhase(int aK0) const;

         tREAL8 K2Rho (int aK) const;
         tREAL8 K2Teta(int aK) const;
         int Rho2K (tREAL8 aR) const;
	 cPt2dr  KTetaRho2Im(const cPt2di & aKTetaRho) const;

	 tREAL8 RhoOfWeight(const tREAL8 &) const;

	 const cExtractedEllipse & mEE;
         const cDataIm2D<tREAL4> & mDIm;
	 const cFullSpecifTarget & mSpec;


	 bool                      mOK;
	 const int                 mPixPerB; ///< number of pixel for each bit to decode
	 const int                 mNbRho;   ///< number of pixel for rho
	 int                       mNbB;     ///< number of bits in code
	 int                       mNbTeta;  ///< number of pixel for each bit to decode
	 tREAL8                    mRho0;
	 tREAL8                    mRho1;
         cIm2D<tREAL4>             mImPolar;
         cDataIm2D<tREAL4> &       mDIP;
         cIm1D<tREAL4>             mAvg;
         cDataIm1D<tREAL4> &       mDAvg;
	 int                       mKR0;
	 int                       mKR1;
	 int                       mPhase0;
	 tREAL8                    mBlack;
	 tREAL8                    mWhite;
	 tREAL8                    mBWAmpl;
	 tREAL8                    mBWAvg;
	 const cOneEncoding *      mEnCode;
};

tREAL8 cCCDecode::Dev(int aK1,int aK2) const
{
    cComputeStdDev<tREAL8> aCS;
    for (int aK=aK1 ; aK<aK2 ; aK++)
    {
         aCS.Add(mDAvg.GetV(aK%mNbTeta));
    }
    return aCS.StdDev(0);
}

tREAL8 cCCDecode::Avg(int aK1,int aK2) const
{
    tREAL8 aSom =0 ;
    for (int aK=aK1 ; aK<aK2 ; aK++)
    {
          aSom += mDAvg.GetV(aK%mNbTeta);
    }
    return aSom / (aK2-aK1);
}

tREAL8 cCCDecode::DevOfPhase(int aK0) const
{
    tREAL8 aSum=0;
    for (int aKBit=0 ; aKBit<mNbB ; aKBit++)
    {
        int aK1 = aK0+aKBit*mPixPerB;
        aSum +=  Dev(aK1+1,aK1+mPixPerB-1);
    }

    return aSum / mNbB;
}

void cCCDecode::ComputePhaseTeta() 
{
    cWhichMin<int,tREAL8> aMinDev;

    for (int aK0=0 ;aK0< mPixPerB ; aK0++)
	    aMinDev.Add(aK0,DevOfPhase(aK0));

    mPhase0 = aMinDev.IndexExtre();

    if (     (aMinDev.ValExtre() > 0.1 * Dev(0,mNbTeta))
          || (aMinDev.ValExtre() > 0.05 *  mBWAmpl)
       )
    {
        mOK = false;
	return;
    }
}

void cCCDecode::ComputeCode(bool Show)
{
    size_t aFlag=0;
    for (int aKBit=0 ; aKBit<mNbB ; aKBit++)
    {
        int aK1 = mPhase0+aKBit*mPixPerB;
        tREAL8 aMoy =  Avg(aK1+1,aK1+mPixPerB-1);

	if (mSpec.BitIs1(aMoy>mBWAvg))
           aFlag |= (1<<aKBit);
    }


    if (! mSpec.AntiClockWiseBit())
       aFlag = BitMirror(aFlag,1<<mSpec.NbBits());

    mEnCode = mSpec.EncodingFromCode(aFlag);


    if (Show)
    {
       StdOut() << "Adr=" << mEnCode << " ";
       if (mEnCode) 
            StdOut() << " Name=" << mEnCode->Name()  
		     << " Code=" <<  mEnCode->Code() 
		     << " BF=" << StrOfBitFlag(mEnCode->Code(), 1<<mNbB);
       StdOut() << "\n";
    }
}



cPt2dr cCCDecode::KTetaRho2Im(const cPt2di & aKTR) const
{
     return mEE.mEllipse.PtOfTeta(K2Teta(aKTR.x()),K2Rho(aKTR.y()));
}

tREAL8 cCCDecode::K2Rho(const int aK)  const {return mRho0+ ((mRho1-mRho0)*aK) / mNbRho;}
tREAL8 cCCDecode::K2Teta(const int aK) const {return  (2*M_PI*aK)/mNbTeta;}

int  cCCDecode::Rho2K(const tREAL8 aR)  const 
{
     return round_ni( ((aR-mRho0)/(mRho1-mRho0)) * mNbRho );
}

tREAL8 cCCDecode::RhoOfWeight(const tREAL8 & aW) const
{
	return (1-aW) * mSpec.Rho_1_BeginCode() + aW * mSpec.Rho_2_EndCode();
}

cCCDecode::cCCDecode(const cExtractedEllipse & anEE,const cDataIm2D<tREAL4> & aDIm,const cFullSpecifTarget & aSpec) :
	mEE        (anEE),
	mDIm       (aDIm),
	mSpec      (aSpec),
	mOK        (true),
	mPixPerB   (10),
	mNbRho     (20),
	mNbB       (mSpec.NbBits()),
	mNbTeta    (mPixPerB * mNbB),
	mRho0      ((mSpec.Rho_0_EndCCB()+mSpec.Rho_1_BeginCode()) /2.0),
	mRho1      (mSpec.Rho_2_EndCode() +0.2),
	mImPolar   (cPt2di(mNbTeta,mNbRho)),
        mDIP       (mImPolar.DIm()),
	mAvg       ( mNbTeta,nullptr,eModeInitImage::eMIA_Null ),
	mDAvg      ( mAvg.DIm()),
	mKR0       ( Rho2K(RhoOfWeight(0.25)) ) ,
	mKR1       ( Rho2K(RhoOfWeight(0.75)) ) ,
	mPhase0    (-1),
	mBlack     (mEE.mSeed.mBlack),
	mWhite     (mEE.mSeed.mWhite),
	mBWAmpl    (mWhite-mBlack),
	mBWAvg     ((mBlack+mWhite)/2.0),
	mEnCode    (nullptr)
{

    //  compute a polar image
    for (int aKTeta=0 ; aKTeta < mNbTeta; aKTeta++)
    {
        for (int aKRho=0 ; aKRho < mNbRho; aKRho++)
        {
		cPt2dr aPt = KTetaRho2Im(cPt2di(aKTeta,aKRho));
		tREAL8 aVal = mDIm.DefGetVBL(aPt,-1);
		if (aVal<0)
		{
                   mOK=false;
                   return;
		}

		mDIP.SetV(cPt2di(aKTeta,aKRho),aVal);
        }
    }

    if (!mOK)
       return;

    // compute an image
    for (int aKTeta=0 ; aKTeta < mNbTeta; aKTeta++)
    {
        std::vector<tREAL8> aVGray;
        for (int aKRho=mKR0 ; aKRho <= mKR1; aKRho++)
	{
            aVGray.push_back(mDIP.GetV(cPt2di(aKTeta,aKRho)));
	}
        mDAvg.SetV(aKTeta,NonConstMediane(aVGray));
    }

    ComputePhaseTeta() ;
    if (!mOK) return;

    ComputeCode(false);
    if (!mOK) return;
}



void  cCCDecode::Show(const std::string & aPrefix)
{
    static int aCpt=0; aCpt++;

    cRGBImage  aIm = RGBImFromGray(mImPolar.DIm(),1.0,9);

    if (mPhase0>=0)
    {
       for (int aKBit=0 ; aKBit<mNbB ; aKBit++)
       {
           tREAL8 aK1 = mPhase0+aKBit*mPixPerB -0.5;

	   aIm.DrawLine(cPt2dr(aK1,0),cPt2dr(aK1,mNbTeta),cRGBImage::Red);

       }
    }
    ComputeCode(true);


    aIm.ToFile(aPrefix + "_ImPolar_"+ToStr(aCpt)+".tif");
}

/*  *********************************************************** */
/*                                                              */
/*             cAppliExtractCodeTarget                          */
/*                                                              */
/*  *********************************************************** */

class cAppliExtractCircTarget : public cMMVII_Appli,
	                        public cAppliParseBoxIm<tREAL4>
{
     public :
        typedef tREAL4              tElemIm;
        typedef cDataIm2D<tElemIm>  tDataIm;
        typedef cImGrad<tElemIm>    tImGrad;


        cAppliExtractCircTarget(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

     private :
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

        int ExeOnParsedBox() override;

	void MakeImageLabel();

	std::string         mNameSpec;
	cFullSpecifTarget * mSpec;
        bool            mVisu;
        cExtract_BW_Ellipse * mExtrEll;
        cParamBWTarget  mPBWT;

        tImGrad         mImGrad;
        cRGBImage       mImVisu;
        cIm2D<tU_INT1>  mImMarq;

        std::vector<cSeedBWTarget> mVSeeds;
	cPhotogrammetricProject     mPhProj;

	std::string                 mPrefixOut;
	bool                        mHasMask;
	std::string                 mNameMask;
};



cAppliExtractCircTarget::cAppliExtractCircTarget
(
    const std::vector<std::string> & aVArgs,
    const cSpecMMVII_Appli & aSpec
) :
   cMMVII_Appli  (aVArgs,aSpec),
   cAppliParseBoxIm<tREAL4>(*this,true,cPt2di(10000,10000),cPt2di(300,300),false) ,
   mSpec    (nullptr),
   mVisu    (true),
   mExtrEll (nullptr),
   mImGrad  (cPt2di(1,1)),
   mImVisu  (cPt2di(1,1)),
   mImMarq  (cPt2di(1,1)),
   mPhProj  (*this)

{
}

        // cExtract_BW_Target * 
cCollecSpecArg2007 & cAppliExtractCircTarget::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   // Standard use, we put args of  cAppliParseBoxIm first
   return
             APBI_ArgObl(anArgObl)
        <<   Arg2007(mNameSpec,"XML name for bit encoding struct")

                   //  << AOpt2007(mDiamMinD, "DMD","Diam min for detect",{eTA2007::HDV})
   ;
}

cCollecSpecArg2007 & cAppliExtractCircTarget::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return APBI_ArgOpt
          (
                anArgOpt
             << mPhProj.DPMask().ArgDirInOpt("TestMask","Mask for selecting point used in detailed mesg/output")
             << AOpt2007(mPBWT.mMinDiam,"DiamMin","Minimum diameters for ellipse",{eTA2007::HDV})
             << AOpt2007(mPBWT.mMaxDiam,"DiamMax","Maximum diameters for ellipse",{eTA2007::HDV})
          );
}


void cAppliExtractCircTarget::MakeImageLabel()
{
    const cExtract_BW_Target::tDImMarq&     aDMarq =  mExtrEll->DImMarq();
    for (const auto & aPix : aDMarq)
    {
         if (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eTmp))
            mImVisu.SetRGBPix(aPix,cRGBImage::Green);

         if (     (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eBadZ))
               || (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eElNotOk))
            )
            mImVisu.SetRGBPix(aPix,cRGBImage::Blue);
         if (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eBadFr))
            mImVisu.SetRGBPix(aPix,cRGBImage::Cyan);
         if (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eBadEl))
            mImVisu.SetRGBPix(aPix,cRGBImage::Red);
         if (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eAverEl))
            mImVisu.SetRGBPix(aPix,cRGBImage::Orange);
         if (aDMarq.GetV(aPix)==tU_INT1(eEEBW_Lab::eBadTeta))
            mImVisu.SetRGBPix(aPix,cRGBImage::Yellow);
    }

    for (const auto & aSeed : mExtrEll->VSeeds())
    {
        if (aSeed.mOk)
        {
           mImVisu.SetRGBPix(aSeed.mPixW,cRGBImage::Red);
           mImVisu.SetRGBPix(aSeed.mPixTop,cRGBImage::Yellow);
        }
        else
        {
           mImVisu.SetRGBPix(aSeed.mPixW,cRGBImage::Yellow);
        }
    }
    mImVisu.ToFile(mPrefixOut + "_Label.tif");
}




int cAppliExtractCircTarget::ExeOnParsedBox()
{
   double aT0 = SecFromT0();

   mExtrEll = new cExtract_BW_Ellipse(APBI_Im(),mPBWT,mPhProj.MaskWithDef(mNameIm,CurBoxIn(),false));
   if (mVisu)
      mImVisu =  cRGBImage::FromFile(mNameIm,CurBoxIn());

   double aT1 = SecFromT0();
   mExtrEll->ExtractAllSeed();
   double aT2 = SecFromT0();
   mExtrEll->AnalyseAllConnectedComponents(mNameIm);
   double aT3 = SecFromT0();

   StdOut() << "TIME-INIT " << aT1-aT0 << "\n";
   StdOut() << "TIME-SEED " << aT2-aT1 << "\n";
   StdOut() << "TIME-CC   " << aT3-aT2 << "\n";

   for (const auto & anEE : mExtrEll->ListExtEl() )
   {
       if (anEE.mSeed.mMarked4Test)
       {
          // anEE.ShowOnFile(mNameIm,35,mPrefixOut);

          cCCDecode aCCD(anEE,APBI_DIm(),*mSpec);
	  aCCD.Show(mPrefixOut);
          // ShowCode(anEE,APBI_DIm(),14);
       }
   }

   if (mVisu)
      MakeImageLabel();

   delete mExtrEll;

   return EXIT_SUCCESS;
}



int  cAppliExtractCircTarget::Exe()
{
   mPrefixOut = "CircTarget_" +  Prefix(APBI_NameIm());

   mSpec = cFullSpecifTarget::CreateFromFile(mNameSpec);

   mPhProj.FinishInit();

   mHasMask =  mPhProj.ImageHasMask(APBI_NameIm()) ;
   if (mHasMask)
      mNameMask =  mPhProj.NameMaskOfImage(APBI_NameIm());


   APBI_ExecAll();  // run the parse file  SIMPL

   StdOut() << "MAK=== " <<   mHasMask << " " << mNameMask  << "\n";

   delete mSpec;
   return EXIT_SUCCESS;
}

/* =============================================== */
/*                                                 */
/*                       ::                        */
/*                                                 */
/* =============================================== */

tMMVII_UnikPApli Alloc_ExtractCircTarget(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliExtractCircTarget(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecExtractCircTarget
(
     "CodedTargetCircExtract",
      Alloc_ExtractCircTarget,
      "Extract coded target from images",
      {eApF::ImProc,eApF::CodedTarget},
      {eApDT::Image,eApDT::Xml},
      {eApDT::Xml},
      __FILE__
);


};

