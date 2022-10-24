#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"
#include "MMVII_Geom2D.h"

namespace MMVII
{

/**  This class make a conversion between pixel space and real space using a R^2->R^2 map,
 * frequently used wih homothety when we do sampling*/

	
template <class TypeMap>  class  cMapPixelization
{
        public :
            typedef typename TypeMap::tTypeElem   tTypeElem;
            typedef cPtxd<tTypeElem,2>            tPtR;
            typedef cPtxd<int,2>                  tPixel;  // use typedef, maybe later in 3D 

	    cMapPixelization(const TypeMap & aMap) : mMap (aMap) {}

            inline tPixel  ToPix(const tPtR & aPtR) const   {return  ToI(mMap.Value(aPtR));  }

        private :
	    TypeMap  mMap;
};

class cSetVisibility : public cDataBoundedSet<tREAL8,3> 
{
    public :
        cSetVisibility(cSensorImage * aSens) : 
            cDataBoundedSet<tREAL8,3>(cBox3dr::BigBox()),
            mSens (aSens) 
	{}
        tREAL8 Insideness(const tPt & aP) const {return mSens->IsVisible(aP) ? 1 : -1;}
    private :
	  cSensorImage * mSens;
};

/**   This abstract class is used to decribe an object containing many triangles.
 *
 * In Z-Buffer, we can use an explicit mesh, but also an implicit one if we parse an image
 * where each pixel is made from two triangle. This implicit class allow to maipulate the two
 * object in the same interface (an avoid converting the pixels in hundred million of triangles ...)
 */

class  cTri3DIterator
{
     public :
        virtual bool GetNextTri(cTri3dR &) = 0;
        virtual bool GetNextPoint(cPt3dr &) = 0;
        virtual void ResetTri()  = 0;
        virtual void ResetPts()  = 0;

        void ResetAll() ;
     private :
};

/** in many case, the implementation can be done by counters */

class cCountTri3DIterator : public cTri3DIterator
{
     public :
        cCountTri3DIterator(size_t aNbP,size_t aNbF);

	virtual cPt3dr  KthP(int aKP) const = 0;
	virtual cTri3dR KthF(int aKF) const = 0;

        bool GetNextTri(cTri3dR &) override;
        bool GetNextPoint(cPt3dr &) override;
        void ResetTri()  override;
        void ResetPts()  override;
     private :
	size_t  mNbP;
	size_t  mNbF;
	size_t  mIndexF;
	size_t  mIndexP;
};

class cMeshTri3DIterator : public cCountTri3DIterator
{
     public :
        cMeshTri3DIterator(cTriangulation3D<tREAL8> *);

	cPt3dr  KthP(int aKP) const override;
	cTri3dR KthF(int aKF) const override;
     private :
	cTriangulation3D<tREAL8> *  mTri;
};

/* =============================================== */
/*                                                 */
/*                 cTri3DIterator                  */
/*                                                 */
/* =============================================== */

void cTri3DIterator::ResetAll()
{
    ResetTri();
    ResetPts();
}

/* =============================================== */
/*                                                 */
/*            cCountTri3DIterator                  */
/*                                                 */
/* =============================================== */

cCountTri3DIterator::cCountTri3DIterator(size_t aNbP,size_t aNbF) :
    mNbP  (aNbP),
    mNbF  (aNbF)
{
   ResetPts();
   ResetTri();
}

void cCountTri3DIterator::ResetTri() { mIndexF=0;}
void cCountTri3DIterator::ResetPts() { mIndexP=0;}

bool cCountTri3DIterator::GetNextPoint(cPt3dr & aP )
{
    if (mIndexP>=mNbP) return false;
    aP = KthP(mIndexP);
    mIndexP++;
    return true;
}

bool cCountTri3DIterator::GetNextTri(cTri3dR & aTri)
{
    if (mIndexF>=mNbF) return false;
    aTri = KthF(mIndexF);
    mIndexF++;
    return true;
}

/* =============================================== */
/*                                                 */
/*              cMeshTri3DIterator                 */
/*                                                 */
/* =============================================== */

cMeshTri3DIterator::cMeshTri3DIterator(cTriangulation3D<tREAL8> * aTri) :
    cCountTri3DIterator(aTri->NbPts(),aTri->NbFace()),
    mTri (aTri)
{
}

cPt3dr  cMeshTri3DIterator::KthP(int aKP) const {return mTri->KthPts(aKP);}
cTri3dR cMeshTri3DIterator::KthF(int aKF) const {return mTri->KthTri(aKF);}


/* =============================================== */
/* =============================================== */
/* =============================================== */


class  cZBuffer
{
      public :

          typedef tREAL4                            tElem;
          typedef cDataIm2D<tElem>                  tDIm;
          typedef cIm2D<tElem>                      tIm;
	  typedef cIm2D<tINT1>                      tImSign;
	  typedef cDataIm2D<tINT1>                  tDImSign;

          typedef cDataInvertibleMapping<tREAL8,3>  tMap;
          typedef cDataBoundedSet<tREAL8,3>         tSet;

	  static constexpr tElem mInfty =  -1e20;

          cZBuffer(cTri3DIterator & aMesh,const tSet & aSetIn,const tMap & aMap,const tSet & aSetOut,double aResolOut);

          const cPt2di  SzPix() ;
          void MakeOneTri(const cTri3dR &);

	  void MakeZBuf
               (
	       );
      private :
          cZBuffer(const cZBuffer & ) = delete;

	  cPt2dr  ToPix(const cPt3dr&) const;
	  cTri3DIterator & mMesh;
	  const tMap &     mMapI2O;
          const tSet &     mSetIn;
          const tSet &     mSetOut;
	  double           mResolOut;

          cBox3dr          mBoxIn;     ///< Box in input space, not sure usefull, but ....
          cBox3dr          mBoxOut;    ///< Box in output space, usefull for xy, not sure for z , but ...
	  cHomot2D<tREAL8> mROut2Pix;  ///<  Mapping Out Coord -> Pix Coord
	  tIm              mZBuf;
	  tImSign          mImSign;   ///< sign of normal  1 or -1 , 0 if uninit
	  tDImSign  *      mDImSign;
          cPt2di           mSzPix;
};

cZBuffer::cZBuffer(cTri3DIterator & aMesh,const tSet &  aSetIn,const tMap & aMapI2O,const tSet &  aSetOut,double aResolOut) :
    mMesh     (aMesh),
    mMapI2O   (aMapI2O),
    mSetIn    (aSetIn),
    mSetOut   (aSetOut),
    mResolOut (aResolOut),

    mBoxIn    (cBox3dr::Empty()),
    mBoxOut   (cBox3dr::Empty()),
    mROut2Pix (),
    mZBuf     (cPt2di(1,1)),
    mImSign   (cPt2di(1,1)),
    mDImSign  (nullptr)
{
    cTplBoxOfPts<tREAL8,3> aBoxOfPtsIn;
    cTplBoxOfPts<tREAL8,3> aBoxOfPtsOut;

    //  compute the box in put and output space
    cPt3dr aPIn;

    mMesh.ResetAll();
    int aCptTot=0;
    int aCptIn=0;
    while (mMesh.GetNextPoint(aPIn))
    {
        aCptTot++;
	if (mSetIn.InsideWithBox(aPIn))
	{
            cPt3dr aPOut = mMapI2O.Value(aPIn);

	    if (mSetOut.InsideWithBox(aPOut))
	    {
               aCptIn++;
               aBoxOfPtsIn.Add(aPIn);
               aBoxOfPtsOut.Add(aPOut);
	    }
	}
    }
    // StdOut() << " cCCCCCCCCCCC " << aCptIn  << " " << aCptTot << "\n";
    mMesh.ResetPts();

    mBoxIn = aBoxOfPtsIn.CurBox();
    mBoxOut = aBoxOfPtsOut.CurBox();

    //   aP0/aResout + aTr -> 1,1
    cPt2dr aTr = cPt2dr(1,1) - Proj(mBoxOut.P0()) * (1.0/mResolOut);
    mROut2Pix = cHomot2D<tREAL8>(aTr,1.0/mResolOut);

    mSzPix =  Pt_round_up(ToPix(mBoxOut.P1()));

    StdOut() << "SZPIX " << mSzPix << " BOX=" << mBoxOut.P0() << " " << mBoxOut.P1() << "\n";

    mZBuf = tIm(mSzPix);
    mZBuf.DIm().InitCste(mInfty);
    mImSign = tImSign(mSzPix,nullptr,eModeInitImage::eMIA_Null);
    mDImSign =  &(mImSign.DIm());
}

cPt2dr  cZBuffer::ToPix(const cPt3dr & aPt) const {return mROut2Pix.Value(Proj(aPt));}

void cZBuffer::MakeZBuf
     (
     )
{
    cTri3dR  aTriIn = cTri3dR::Tri000();
    while (mMesh.GetNextTri(aTriIn))
    {
         //  not sure this us to test that, or the user to assure it give clean data ...
         if ((aTriIn.Regularity() >0)  && mSetIn.InsideWithBox(aTriIn))
	 {
             cTri3dR aTriOut = mMapI2O.TriValue(aTriIn);
	     
             if ((aTriOut.Regularity() >0)  && mSetOut.InsideWithBox(aTriOut))
	     {
                  MakeOneTri(aTriOut);
	     }
	 }
    }
}

void cZBuffer::MakeOneTri(const cTri3dR &aTri3)
{

    cTriangle2DCompiled<tREAL8>  aTri2(ToPix(aTri3.Pt(0)) , ToPix(aTri3.Pt(1)) ,ToPix(aTri3.Pt(2)));

    StdOut() << " tttt " << ToPix(aTri3.Pt(0)) <<  ToPix(aTri3.Pt(1)) << ToPix(aTri3.Pt(2)) << "\n";

    std::vector<cPt2di> aVPix;
    std::vector<cPt3dr> aVW;

     aTri2.PixelsInside(aVPix,1e-8,&aVW);

    StdOut() << " qqqqtttt " << aVPix.size() << "\n";

     cPt3dr aNorm = Normal(aTri3);
     int aSign = (aNorm.z() > 0) ? 1 : - 1;

     static int aCpt     =0;
     static int aCptPlus =0;
     aCpt++;
     if (aSign==1) aCptPlus++;
     StdOut() << " CPT Sign Normal " << aCpt << " +=" << aCptPlus << "\n";

     for (size_t aK=0 ; aK<aVPix.size() ; aK++)
     {
          
     }

}

/* =============================================== */
/*                                                 */
/*                 cAppliCloudClip                 */
/*                                                 */
/* =============================================== */

/** Application for projecting a mesh on image, for each triangle,
 * indicate if it is visible and what is its quality=resolution in lowest dim
*/

class cAppliProMeshImage : public cMMVII_Appli
{
     public :

        cAppliProMeshImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

     private :
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

     // --- Mandatory ----
	std::string mNameCloud3DIn;
	std::string mNameIm;
	std::string mNameOri;


     // --- Optionnal ----
	std::string mNameCloud2DIn;

     // --- constructed ---
        cPhotogrammetricProject   mPhProj;
        cTriangulation3D<tREAL8>* mTri3D;
        cSensorCamPC *            mCamPC;
};

cCollecSpecArg2007 & cAppliProMeshImage::ArgObl(cCollecSpecArg2007 & anArgObl) 
{
   return anArgObl
	  <<   Arg2007(mNameCloud3DIn,"Name of input cloud/mesh", {eTA2007::FileDirProj,eTA2007::FileCloud,eTA2007::Input})
	  <<   Arg2007(mNameIm,"Name of image", {eTA2007::FileImage,eTA2007::OptionalExist})
	  <<   mPhProj.OriInMand()

   ;
}

cAppliProMeshImage::cAppliProMeshImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
   cMMVII_Appli     (aVArgs,aSpec),
   mPhProj          (*this),
   mTri3D           (nullptr),
   mCamPC           (nullptr)
{
}


cCollecSpecArg2007 & cAppliProMeshImage::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return anArgOpt
           << AOpt2007(mNameCloud2DIn,"M2","Mesh 2D, dev of cloud 3D, to generate a visu of hiden part ", {eTA2007::FileCloud,eTA2007::Input})
   ;

}


int cAppliProMeshImage::Exe() 
{
   mPhProj.FinishInit();

   mTri3D = new cTriangulation3D<tREAL8>(DirProject()+mNameCloud3DIn);
   cMeshTri3DIterator  aTriIt(mTri3D);

   mCamPC = mPhProj.AllocCamPC(mNameIm,true);
   cSIMap_Ground2ImageAndProf aMapCamDepth(mCamPC);

   cSetVisibility aSetVis(mCamPC);

   double Infty =1e20;
   cPt2di aSzPix = mCamPC->SzPix();
   // StdOut() << "SZPIX " << aSzPix << "\n";
   cBox3dr  aBox(cPt3dr(0,0,-Infty),cPt3dr(aSzPix.x(),aSzPix.y(),Infty));
   cDataBoundedSet<tREAL8,3>  aSetCam(aBox);

   

   StdOut() << "FOCALE "  << mCamPC->InternalCalib()->F() << " " << &aMapCamDepth << "\n";

   cZBuffer aZBuf(aTriIt,aSetVis,aMapCamDepth,aSetCam,3.0);
   aZBuf.MakeZBuf();

   delete mTri3D;
   return EXIT_SUCCESS;
}

     /* =============================================== */
     /*                       ::                        */
     /* =============================================== */

tMMVII_UnikPApli Alloc_ProMeshImage(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliProMeshImage(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecProMeshImage
(
     "0_MeshProjImage",
      Alloc_ProMeshImage,
      "(internal) Project a mes on an image",
      {eApF::Cloud},
      {eApDT::Ply,eApDT::Orient},
      {eApDT::FileSys},
      __FILE__
);

#if (0)
#endif

}
