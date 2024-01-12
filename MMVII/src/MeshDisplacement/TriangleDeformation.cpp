#include "cMMVII_Appli.h"

#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"

#include "MMVII_PhgrDist.h"

#include "MMVII_TplSymbTriangle.h"

using namespace NS_SymbolicDerivative;

/**
   \file TriangleDeformation.cpp

   \brief file for generating random points distributed uniformely
   and applying 2D Delaunay triangulation.
**/

namespace MMVII
{
    /********************************************/
    /*                                          */
    /*            cTriangleDeformation          */
    /*                                          */
    /********************************************/

    class cAppli_cTriangleDeformation : public cMMVII_Appli
    {
    public:
        typedef cIm2D<tREAL8> tIm;
        typedef cDataIm2D<tREAL8> tDIm;

        cAppli_cTriangleDeformation(const std::vector<std::string> &aVArgs,
                                    const cSpecMMVII_Appli &aSpec);
        ~cAppli_cTriangleDeformation();

        int Exe() override;
        cCollecSpecArg2007 &ArgObl(cCollecSpecArg2007 &anArgObl) override;
        cCollecSpecArg2007 &ArgOpt(cCollecSpecArg2007 &anArgOpt) override;

        void ConstructUniformRandomVectorAndApplyDelaunay();
        void GeneratePointsForDelaunay();
        void OneIterationFitModele(bool aIsLastIter);
        cPt2dr ApplyBarycenterTranslationFormulaToFilledPixel(cHomot2D<tREAL8> & aCurrentTranslationPointA, cHomot2D<tREAL8> & aCurrentTranslationPointB, 
                                                              cHomot2D<tREAL8> & aCurrentTranslationPointC, std::vector<double> & aVObs); // new addition

    private:
        // ==  Mandatory args ====
        std::string mNamePreImage; //< Name of given pre-image.
        std::string mNamePostImage; //< Name of given post-image.
        int mNumberPointsToGenerate;	 // number of generated points
		int mRandomUniformLawUpperBound; // Uniform law generate numbers from [0, mRandomUniformLawUpperBound [
        int mNumberOfOptimisationIterations; // number of iterations in optimisation process

        // ==  Optionnal args ====
        bool mShow; // print result, export image ...
        bool mGenerateDisplacementImage;

        // ==  Internal variables ====

        cPt2di mSzImPre; ///<  size of image.
        tIm mImPre;      ///<  memory representation of the image.
        tDIm *mDImPre;   ///<  memory representation of the image.

        cPt2di mSzImPost; ///<  size of image.
        tIm mImPost;      ///<  memory representation of the image.
        tDIm *mDImPost;   ///<  memory representation of the image.

        cPt2di mSzImOut; ///<  size of image.
        tIm mImOut;      ///<  memory representation of the image.
        tDIm *mDImOut;   ///<  memory representation of the image.

        std::vector<cPt2dr> mVectorPts; // A vector containing a set of points.
        cTriangulation2D<tREAL8> mDelTri; // A Delaunay triangle.
        std::vector<cPt2dr> mFinalTranslatedPixelCoords; // Final translation coefficients.

        cResolSysNonLinear<tREAL8> *mSys; ///< Non Linear Sys for solving problem.
        cCalculator<double> *mEqHomTri;   ///< calculator giving access to values and derivatives.
    };

    cAppli_cTriangleDeformation::cAppli_cTriangleDeformation(const std::vector<std::string> &aVArgs,
                                                             const cSpecMMVII_Appli &aSpec) : cMMVII_Appli(aVArgs, aSpec),
                                                                                              mShow(true),
                                                                                              mGenerateDisplacementImage(false),
                                                                                              mSzImPre(cPt2di(1, 1)),
                                                                                              mImPre(mSzImPre),
                                                                                              mDImPre( nullptr),
                                                                                              mSzImPost(cPt2di(1, 1)),
                                                                                              mImPost(mSzImPost),
                                                                                              mDImPost(nullptr),
                                                                                              mSzImOut(cPt2di(1, 1)),
                                                                                              mImOut(mSzImOut),
                                                                                              mDImOut(nullptr),
                                                                                              mVectorPts({cPt2dr(0, 0)}),
                                                                                              mDelTri(mVectorPts),
                                                                                              mSys(nullptr),
                                                                                              mEqHomTri(nullptr)
    {
        cDenseVect<double> aVInit(6);
        /*
        aVInit(0) = 0.05;
        aVInit(1) = 0.01;
        aVInit(2) = 0.08;
        aVInit(3) = 0.06;
        aVInit(4) = 0.07;
        aVInit(5) = 0.04;
        */
        aVInit(0) = 0;
        aVInit(1) = 0;
        aVInit(2) = 0;
        aVInit(3) = 0;
        aVInit(4) = 0;
        aVInit(5) = 0;

        mSys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqDense, aVInit);
        mEqHomTri = EqDeformTriAffinity(true, 1); // true means with derivative,  1 is size of buffer
    }

    cAppli_cTriangleDeformation::~cAppli_cTriangleDeformation()
    {
        delete mSys;
        delete mEqHomTri;
    }

    cCollecSpecArg2007 &cAppli_cTriangleDeformation::ArgObl(cCollecSpecArg2007 &anArgObl)
    {
        return anArgObl
               << Arg2007(mNamePreImage, "Name of pre-image file.", {{eTA2007::FileImage}, {eTA2007::FileDirProj}})
               << Arg2007(mNamePostImage, "Name of post-image file.", {eTA2007::FileImage})
               << Arg2007(mNumberPointsToGenerate, "Number of points you want to generate for triangulation.")
			   << Arg2007(mRandomUniformLawUpperBound, "Maximum value that the uniform law can draw from.")
               << Arg2007(mNumberOfOptimisationIterations, "Total number of iterations to run in optimisation process.");
    }

    cCollecSpecArg2007 &cAppli_cTriangleDeformation::ArgOpt(cCollecSpecArg2007 &anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mShow, "Show", "Whether to print minimisation results.", {eTA2007::HDV})
               << AOpt2007(mGenerateDisplacementImage, "GenerateDisplacementImage", "Whether to generate and save an image having been translated.", {eTA2007::HDV});
    }

	void cAppli_cTriangleDeformation::ConstructUniformRandomVectorAndApplyDelaunay()
	{
        mVectorPts.pop_back(); // eliminate initialisation values.
		// Generate coordinates from drawing lines and columns of coordinates from a uniform distribution.
		for (int aNbPt = 0; aNbPt < mNumberPointsToGenerate; aNbPt++)
		{
			const double aUniformRandomLine = RandUnif_N(mRandomUniformLawUpperBound);
			const double aUniformRandomCol = RandUnif_N(mRandomUniformLawUpperBound);
			const cPt2dr aUniformRandomPt(aUniformRandomLine, aUniformRandomCol);
			mVectorPts.push_back(aUniformRandomPt);
		}
        mDelTri = mVectorPts;

		mDelTri.MakeDelaunay(); // Delaunay triangulate randomly generated points.
	}

	void cAppli_cTriangleDeformation::GeneratePointsForDelaunay()
	{
		const int aMinimumLinCol = std::min(mSzImPre.y(), mSzImPre.x());

		// make sure that values greater than image size can't be drawn from uniform law.
		while (mRandomUniformLawUpperBound >= aMinimumLinCol)
		{
			StdOut() << "Maximum value drawn from uniform law needs to be smaller than " << aMinimumLinCol << "." << std::endl;
			std::cin >> mRandomUniformLawUpperBound;
		}

        ConstructUniformRandomVectorAndApplyDelaunay();
	}

    cPt2dr cAppli_cTriangleDeformation::ApplyBarycenterTranslationFormulaToFilledPixel(cHomot2D<tREAL8> & aCurrentTranslationPointA, cHomot2D<tREAL8> & aCurrentTranslationPointB, 
                                                                                       cHomot2D<tREAL8> & aCurrentTranslationPointC, std::vector<double> & aVObs) // new addition
    {
        // auto aXTri = aXCoordinates + aAlphaCoordinate * aGeomTrXPointA + aBetaCoordinate * aGeomTrXPointB + aGammaCoordinate * aGeomTrXPointC;
        // auto aYTri = aYCoordinates + aAlphaCoordinate * aGeomTrYPointA + aBetaCoordinate * aGeomTrYPointB + aGammaCoordinate * aGeomTrYPointC;
        // size_t aIndObs = 0;
        // apply barycenter translation formula for x and y.
        auto aXTri = aVObs[0] + aVObs[2] * aCurrentTranslationPointA.Tr().x() + aVObs[3] * aCurrentTranslationPointB.Tr().x() + 
                     aVObs[4] * aCurrentTranslationPointC.Tr().x();
        auto aYTri = aVObs[1] + aVObs[2] * aCurrentTranslationPointA.Tr().y() + aVObs[3] * aCurrentTranslationPointB.Tr().y() + 
                     aVObs[4] * aCurrentTranslationPointC.Tr().y();

        cPt2dr aComputedTranslatedPixel = cPt2dr(aXTri, aYTri);

        return aComputedTranslatedPixel;
    }

    void cAppli_cTriangleDeformation::OneIterationFitModele(bool aIsLastIter)
    {
        //----------- index of unkown, basic here because the unknown are the same for each equation
        std::vector<int> aVecInd{0, 1, 2, 3, 4, 5};
        //----------- allocate vec of obs :
        std::vector<double> aVObs(12, 0.0); // 6 for ImagePre interpolation and 6 for ImagePost

        //----------- extract current parameters
        cDenseVect<double> aVCur = mSys->CurGlobSol();
        // new addition
        cHomot2D<tREAL8> aCurTrPointA(cPt2dr(aVCur(0), aVCur(1)), 0);    // current homothety translation point A
        cHomot2D<tREAL8> aCurTrPointB(cPt2dr(aVCur(2), aVCur(3)), 0);    // current homothety translation point B
        cHomot2D<tREAL8> aCurTrPointC(cPt2dr(aVCur(4), aVCur(5)), 0);    // current homothety translation point C
        // double aCurScR = aVCur(0);                                    // current scale on radiometry
        // double aCurTrR = aVCur(1);                                    // current translation on radiometry

        //----------- declaration of indicator of convergence
        double aSomDif = 0; // sum of difference between untranslated pixel and translated one.
        double aSomMod = 0; // sum of value of model, to normalize the difference
        double aNbOut = 0;  // number of translated pixels out of image

        // Parse all the point to add the observations on each point
        for (size_t aTr = 0; aTr < mDelTri.NbFace(); aTr++)
        {
            const cTriangle<tREAL8, 2> aTri = mDelTri.KthTri(aTr);

            const cTriangle2DCompiled aCompTri(aTri);

            std::vector<cPt2di> aVectorToFillWithInsidePixels;
		    aCompTri.PixelsInside(aVectorToFillWithInsidePixels);

            // long unsigned is necessary as there can be a lot of pixels in triangles.
            for (long unsigned int aFilledPixel=0; aFilledPixel < aVectorToFillWithInsidePixels.size(); aFilledPixel++) 
            {
                // prepare for barycenter translation formula by filling aVObs with different coordinates.
                FormalInterpBarycenter_SetObs(aVObs, 0, aCompTri, aVectorToFillWithInsidePixels, aFilledPixel, *mDImPre);

                // image of a point in triangle by current homothety
                cPt2dr aTranslatedFilledPoint = ApplyBarycenterTranslationFormulaToFilledPixel(aCurTrPointA, aCurTrPointB, aCurTrPointC, aVObs);
                if (mDImPost->InsideBL(aTranslatedFilledPoint)) // avoid error
                {
                    // prepare for application of bilinear formula
                    FormalBilinTri_SetObs(aVObs, TriangleDisplacement_NbObs, aTranslatedFilledPoint, *mDImPost);

                    // Now add observation
                    mSys->CalcAndAddObs(mEqHomTri, aVecInd, aVObs);

                    // compute indicators
                    double aDif = mDImPre->GetVBL(cPt2dr(aVObs[0], aVObs[1])) - mDImPost->GetVBL(aTranslatedFilledPoint); // residual
                    aSomMod += mDImPost->GetV(cPt2di(aVObs[0], aVObs[1]));
                    aSomDif += std::abs(aDif);

                    if (aIsLastIter && mGenerateDisplacementImage)
                        mFinalTranslatedPixelCoords.push_back(aTranslatedFilledPoint);
                }
                else
                    aNbOut++;
            }
        }

        // Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();

        if (mShow)
            StdOut() << aSomDif / aSomMod << ", " << aNbOut << std::endl;

        if (aIsLastIter)
        {
            cDenseVect<double> aVFinalSol = mSys->CurGlobSol();
            // new addition
            cHomot2D<tREAL8> aLastTrPointA(cPt2dr(aVFinalSol(0), aVFinalSol(1)), 0);    // final homothety translation point A
            cHomot2D<tREAL8> aLastTrPointB(cPt2dr(aVFinalSol(2), aVFinalSol(3)), 0);    // final homothety translation point B
            cHomot2D<tREAL8> aLastTrPointC(cPt2dr(aVFinalSol(4), aVFinalSol(5)), 0);    // final homothety translation point C

            StdOut() << "Final translation on x-axis of Point A : " << aLastTrPointA.Tr().x() << " and " 
                     << aLastTrPointA.Tr().y() << " for y-axis." << std::endl;
            StdOut() << "Final translation on x-axis of Point B : " << aLastTrPointB.Tr().x() << " and " 
                     << aLastTrPointB.Tr().y() << " for y-axis." << std::endl;
            StdOut() << "Final translation on x-axis of Point C : " << aLastTrPointC.Tr().x() << " and " 
                     << aLastTrPointC.Tr().y() << " for y-axis." << std::endl;
        }
    }

    int cAppli_cTriangleDeformation::Exe()
    {
        // read pre and post images and their sizes
        mImPre = tIm::FromFile(mNamePreImage);
        mImPost = tIm::FromFile(mNamePostImage);

        // cDataFileIm2D aDescFile = cDataFileIm2D::Create(mNamePreImage, false);

        mDImPre = &(mImPre.DIm());
        mSzImPre = (mDImPre->Sz());

        mDImPost = &(mImPost.DIm());
        mSzImPost = mDImPost->Sz();

        /*
        mDImOut = mImOut.DIm();
        mSzImOut = mDImPost.Sz();
        */
        if (mShow)
            StdOut() << "Diff," << "NbOut" << std::endl;

        GeneratePointsForDelaunay();

        // const int aNbIters = 20;
        for (int aIter = 0; aIter < mNumberOfOptimisationIterations; aIter++)
            OneIterationFitModele(aIter == (mNumberOfOptimisationIterations - 1));

        if (mGenerateDisplacementImage)
        {
            mImOut = tIm(mSzImPost);
            mDImOut = &mImOut.DIm();

            for (const auto &aNullPix : *mDImOut)
                mDImOut->SetV(aNullPix, 0);

            for (const auto &aDisplacedPix : *mDImOut)
            {
                /*
                tREAL8 aDx = aImDispx.DIm().GetV(aPix);
                tREAL8 aDy = aImDispy.DIm().GetV(aPix);
                cPt2dr aPixR = ToR(aPix) + cPt2dr(aDx, aDy);

                mDImOut->SetV(aPix, mDImIn->DefGetVBL(aPixR, 0));
                */
                cPt2dr aDisplacedPoint = cPt2dr(aDisplacedPix.x(), aDisplacedPix.y());
                mDImOut->SetV(cPt2di(aDisplacedPoint.x(), aDisplacedPoint.y()), mDImPost->DefGetVBL(aDisplacedPoint, 0));
            }
            // cDataFileIm2D aDescFile = cDataFileIm2D::Create(mNamePostImage, false);
            mDImOut->ToFile("DisplacedPixels.tif"); //, aDescFile.Type());
        }
        return EXIT_SUCCESS;
    }

    /********************************************/
    //              ::MMVII                     //
    /********************************************/

    tMMVII_UnikPApli Alloc_cTriangleDeformation(const std::vector<std::string> &aVArgs, const cSpecMMVII_Appli &aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_cTriangleDeformation(aVArgs, aSpec));
    }

    cSpecMMVII_Appli TheSpec_ComputeTriangleDeformation(
        "ComputeTriangleDeformation",
        Alloc_cTriangleDeformation,
        "Compute deformation of triangles between images",
        {eApF::ImProc}, // category
        {eApDT::Image}, // input
        {eApDT::Image}, // output
        __FILE__);

}; // namespace MMVII
