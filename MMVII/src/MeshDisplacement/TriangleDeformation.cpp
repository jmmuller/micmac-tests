#include "cMMVII_Appli.h"

#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"

#include "MMVII_TplSymbTriangle.h"

#include "TriangleDeformation.h"

/**
   \file TriangleDeformation.cpp

   \brief file for generating random points distributed uniformely
   and applying 2D Delaunay triangulation.
**/

namespace MMVII
{
    /****************************************/
    /*                                      */
    /*         cPtInsideTriangles           */
    /*                                      */
    /****************************************/

    cPtInsideTriangles::cPtInsideTriangles(const cTriangle2DCompiled<tREAL8> &aCompTri,              // a compiled triangle
                                           const std::vector<cPt2di> &aVectorFilledwithInsidePixels, // vector containing pixels insisde triangles
                                           const long unsigned int aFilledPixel,                     // a counter that is looping over pixels in triangles
                                           const cDataIm2D<tREAL8> &aDIm)                            // image
    {
        mFilledPoint = cPt2dr(aVectorFilledwithInsidePixels[aFilledPixel].x(), aVectorFilledwithInsidePixels[aFilledPixel].y());
        mBarycenterCoordinatesOfPixel = aCompTri.CoordBarry(mFilledPoint);
        mValueOfPixel = aDIm.GetV(cPt2di(aVectorFilledwithInsidePixels[aFilledPixel].x(),
                                         aVectorFilledwithInsidePixels[aFilledPixel].y()));
    }

    cPt3dr cPtInsideTriangles::GetBarycenterCoordinates() const {return mBarycenterCoordinatesOfPixel;} // Accessor
    cPt2dr cPtInsideTriangles::GetCartesianCoordinates() const {return mFilledPoint;}                   // Accessor
    tREAL8 cPtInsideTriangles::GetPixelValue() const {return mValueOfPixel;}                            // Accessor


    /********************************************/
    /*                                          */
    /*            cTriangleDeformation          */
    /*                                          */
    /********************************************/

    cAppli_cTriangleDeformation::cAppli_cTriangleDeformation(const std::vector<std::string> &aVArgs,
                                                             const cSpecMMVII_Appli &aSpec) : cMMVII_Appli(aVArgs, aSpec),
                                                                                              mRandomUniformLawUpperBound(20),
                                                                                              mShow(true),
                                                                                              mGenerateDisplacementImage(false),
                                                                                              mSzImPre(cPt2di(1, 1)),
                                                                                              mImPre(mSzImPre),
                                                                                              mDImPre(nullptr),
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
        mEqHomTri = EqDeformTriHomothety(true, 1); // true means with derivative, 1 is size of buffer
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
               << Arg2007(mNumberOfOptimisationIterations, "Total number of iterations to run in optimisation process.");
    }

    cCollecSpecArg2007 &cAppli_cTriangleDeformation::ArgOpt(cCollecSpecArg2007 &anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mRandomUniformLawUpperBound, "RandomUniformLawUpperBound", "Maximum value that the uniform law can draw from.", {eTA2007::HDV})
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
        mRandomUniformLawUpperBound = std::min(mSzImPre.y(), mSzImPre.x());

        ConstructUniformRandomVectorAndApplyDelaunay();
    }

    void cAppli_cTriangleDeformation::InitialisationAfterExe()
    {
        cDenseVect<double> aVInit(2 * mDelTri.NbPts(), eModeInitImage::eMIA_Null);

        mSys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqDense, aVInit);
    }

    cPt2dr cAppli_cTriangleDeformation::ApplyBarycenterTranslationFormulaToFilledPixel(tHomot2d &aCurrentTranslationPointA, tHomot2d &aCurrentTranslationPointB,
                                                                                       tHomot2d &aCurrentTranslationPointC, std::vector<double> &aVObs)
    {
        // auto aXTri = aXCoordinates + aAlphaCoordinate * aGeomTrXPointA + aBetaCoordinate * aGeomTrXPointB + aGammaCoordinate * aGeomTrXPointC;
        // auto aYTri = aYCoordinates + aAlphaCoordinate * aGeomTrYPointA + aBetaCoordinate * aGeomTrYPointB + aGammaCoordinate * aGeomTrYPointC;

        // apply barycenter translation formula for x and y on current observations.
        auto aXTri = aVObs[0] + aVObs[2] * aCurrentTranslationPointA.Tr().x() + aVObs[3] * aCurrentTranslationPointB.Tr().x() +
                     aVObs[4] * aCurrentTranslationPointC.Tr().x();
        auto aYTri = aVObs[1] + aVObs[2] * aCurrentTranslationPointA.Tr().y() + aVObs[3] * aCurrentTranslationPointB.Tr().y() +
                     aVObs[4] * aCurrentTranslationPointC.Tr().y();

        cPt2dr aComputedTranslatedPixel = cPt2dr(aXTri, aYTri);

        return aComputedTranslatedPixel;
    }

    void cAppli_cTriangleDeformation::DoOneIteration(bool aIsLastIter)
    {
        //----------- allocate vec of obs :
        std::vector<double> aVObs(12, 0.0); // 6 for ImagePre interpolation and 6 for ImagePost

        //----------- extract current parameters
        cDenseVect<double> aVCur = mSys->CurGlobSol();
        // double aCurScR = aVCur(0);                            // current scale on radiometry
        // double aCurTrR = aVCur(1);                            // current translation on radiometry

        //----------- declaration of indicator of convergence
        double aSomDif = 0; // sum of difference between untranslated pixel and translated one.
        double aNbOut = 0;  // number of translated pixels out of image
        
        size_t aTotalNumberOfInsidePixels=0;
        
        // Parse all the point to add the observations on each point
        for (size_t aTr = 0; aTr < mDelTri.NbFace(); aTr++)
        {
            // if (mDelTri.ValidFace())
            const tTri2dr aTri = mDelTri.KthTri(aTr);
            const cPt3di aIndicesOfTriKnots = mDelTri.KthFace(aTr);

            const cTriangle2DCompiled aCompTri(aTri);

            std::vector<cPt2di> aVectorToFillWithInsidePixels;
            aCompTri.PixelsInside(aVectorToFillWithInsidePixels);

            //----------- index of unkown, finds the pixels 
            std::vector<int> aVecInd = {2*aIndicesOfTriKnots.x(), 2*aIndicesOfTriKnots.x() + 1, 
                                        2*aIndicesOfTriKnots.y(), 2*aIndicesOfTriKnots.y() + 1,
                                        2*aIndicesOfTriKnots.z(), 2*aIndicesOfTriKnots.z() + 1};

            tHomot2d aCurTrPointA(cPt2dr(aVCur(aVecInd.at(0)), aVCur(aVecInd.at(1))), 0); // current homothety translation point A
            tHomot2d aCurTrPointB(cPt2dr(aVCur(aVecInd.at(2)), aVCur(aVecInd.at(3))), 0); // current homothety translation point B
            tHomot2d aCurTrPointC(cPt2dr(aVCur(aVecInd.at(4)), aVCur(aVecInd.at(5))), 0); // current homothety translation point C

            size_t aNumberOfInsidePixels = aVectorToFillWithInsidePixels.size();
            // size_t is necessary as there can be a lot of pixels in triangles.
            for (size_t aFilledPixel = 0; aFilledPixel < aNumberOfInsidePixels; aFilledPixel++)
            {
                const cPtInsideTriangles aPixInsideTriangle = cPtInsideTriangles(aCompTri, aVectorToFillWithInsidePixels, aFilledPixel, *mDImPre);
                // prepare for barycenter translation formula by filling aVObs with different coordinates.
                FormalInterpBarycenter_SetObs(aVObs, 0, aPixInsideTriangle);

                // image of a point in triangle by current homothety
                cPt2dr aTranslatedFilledPoint = ApplyBarycenterTranslationFormulaToFilledPixel(aCurTrPointA, aCurTrPointB, aCurTrPointC, aVObs);
                if (mDImPost->InsideBL(aTranslatedFilledPoint)) // avoid error
                {
                    // prepare for application of bilinear formula
                    FormalBilinTri_SetObs(aVObs, TriangleDisplacement_NbObs, aTranslatedFilledPoint, *mDImPost);

                    // Now add observation
                    mSys->CalcAndAddObs(mEqHomTri, aVecInd, aVObs);

                    // compute indicators
                    const double aDif = mDImPre->GetVBL(cPt2dr(aVObs[0], aVObs[1])) - mDImPost->GetVBL(aTranslatedFilledPoint); // residual
                    aSomDif += std::abs(aDif);

                    if (aIsLastIter && mGenerateDisplacementImage)
                        mFinalTranslatedPixelCoords.push_back(aTranslatedFilledPoint);
                    
                    aTotalNumberOfInsidePixels += aNumberOfInsidePixels;
                }
                else
                    aNbOut++;
            }
        }

        // Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();

        if (mShow)
            StdOut() << aSomDif / aTotalNumberOfInsidePixels << ", " << aNbOut << std::endl;

        if (aIsLastIter && mShow)
        {
            cDenseVect<double> aVFinalSol = mSys->CurGlobSol();
            // new addition
            tHomot2d aLastTrPointA(cPt2dr(aVFinalSol(0), aVFinalSol(1)), 0); // final homothety translation point 1
            tHomot2d aLastTrPointB(cPt2dr(aVFinalSol(2), aVFinalSol(3)), 0); // final homothety translation point 2
            tHomot2d aLastTrPointC(cPt2dr(aVFinalSol(4), aVFinalSol(5)), 0); // final homothety translation point 3

            StdOut() << "Final translation on x-axis of the 1st point : " << aLastTrPointA.Tr().x() << " and "
                     << aLastTrPointA.Tr().y() << " for y-axis." << std::endl;
            StdOut() << "Final translation on x-axis of the 2nd point : " << aLastTrPointB.Tr().x() << " and "
                     << aLastTrPointB.Tr().y() << " for y-axis." << std::endl;
            StdOut() << "Final translation on x-axis of the 3rd point : " << aLastTrPointC.Tr().x() << " and "
                     << aLastTrPointC.Tr().y() << " for y-axis." << std::endl;
        }

    }

    int cAppli_cTriangleDeformation::Exe()
    {
        // read pre and post images and their sizes
        mImPre = tIm::FromFile(mNamePreImage);
        mImPost = tIm::FromFile(mNamePostImage);

        mDImPre = &(mImPre.DIm());
        mSzImPre = (mDImPre->Sz());

        mDImPost = &(mImPost.DIm());
        mSzImPost = mDImPost->Sz();

        if (mShow)
            StdOut() << "Diff,"
                     << "NbOut" << std::endl;

        GeneratePointsForDelaunay();

        InitialisationAfterExe();

        for (int aIter = 0; aIter < mNumberOfOptimisationIterations; aIter++)
            DoOneIteration(aIter == (mNumberOfOptimisationIterations - 1));

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
