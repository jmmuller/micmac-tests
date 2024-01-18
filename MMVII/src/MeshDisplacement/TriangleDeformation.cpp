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

    cPt3dr cPtInsideTriangles::GetBarycenterCoordinates() const { return mBarycenterCoordinatesOfPixel; } // Accessor
    cPt2dr cPtInsideTriangles::GetCartesianCoordinates() const { return mFilledPoint; }                   // Accessor
    tREAL8 cPtInsideTriangles::GetPixelValue() const { return mValueOfPixel; }                            // Accessor

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
        mVectorPts.pop_back(); // eliminate initialisation values
        // Generate coordinates from drawing lines and columns of coordinates from a uniform distribution
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
        // If user hasn't defined another value than the default value, it is changed
        if (mRandomUniformLawUpperBound == 20)
            // Maximum value of coordinates are drawn from [0, MinimumSizeBetweenLinesAndColumnsOfPreImage[
            mRandomUniformLawUpperBound = std::min(mSzImPre.y(), mSzImPre.x());

        ConstructUniformRandomVectorAndApplyDelaunay();
    }

    void cAppli_cTriangleDeformation::InitialisationAfterExe()
    {
        tDensevect aVInit(2 * mDelTri.NbPts(), eModeInitImage::eMIA_Null);

        mSys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqDense, aVInit);
    }

    cPt2dr cAppli_cTriangleDeformation::ApplyBarycenterTranslationFormulaToFilledPixel(tHomot2d &aCurrentTranslationPointA, tHomot2d &aCurrentTranslationPointB,
                                                                                       tHomot2d &aCurrentTranslationPointC, std::vector<double> &aVObs)
    {
        // auto aXTri = aXCoordinates + aAlphaCoordinate * aGeomTrXPointA + aBetaCoordinate * aGeomTrXPointB + aGammaCoordinate * aGeomTrXPointC;
        // auto aYTri = aYCoordinates + aAlphaCoordinate * aGeomTrYPointA + aBetaCoordinate * aGeomTrYPointB + aGammaCoordinate * aGeomTrYPointC;

        // apply current barycenter translation formula for x and y on current observations.
        auto aXTri = aVObs[0] + aVObs[2] * aCurrentTranslationPointA.Tr().x() + aVObs[3] * aCurrentTranslationPointB.Tr().x() +
                     aVObs[4] * aCurrentTranslationPointC.Tr().x();
        auto aYTri = aVObs[1] + aVObs[2] * aCurrentTranslationPointA.Tr().y() + aVObs[3] * aCurrentTranslationPointB.Tr().y() +
                     aVObs[4] * aCurrentTranslationPointC.Tr().y();

        cPt2dr aComputedTranslatedPixel = cPt2dr(aXTri, aYTri);

        return aComputedTranslatedPixel;
    }

    void cAppli_cTriangleDeformation::LoopOverTrianglesAndUpdateParameters(const bool aIsLastIter)
    {
        //----------- allocate vec of obs :
        std::vector<double> aVObs(12, 0.0); // 6 for ImagePre interpolation and 6 for ImagePost

        //----------- extract current parameters
        tDensevect aVCur = mSys->CurGlobSol();           // Get current solution.
        // double aCurScR = aVCur(0);                            // current scale on radiometry
        // double aCurTrR = aVCur(1);                            // current translation on radiometry
        //----------- declaration of indicator of convergence
        double aSomDif = 0; // sum of difference between untranslated pixel and translated one.
        double aNbOut = 0;  // number of translated pixels out of image

        // Count number of pixels inside triangles for normalisation
        size_t aTotalNumberOfInsidePixels = 0;

        // Loop over all triangles to add the observations on each point
        for (size_t aTr = 0; aTr < mDelTri.NbFace(); aTr++)
        {
            // if (mDelTri.ValidFace())
            const tTri2dr aTri = mDelTri.KthTri(aTr);
            const cPt3di aIndicesOfTriKnots = mDelTri.KthFace(aTr);

            const cTriangle2DCompiled aCompTri(aTri);

            std::vector<cPt2di> aVectorToFillWithInsidePixels;
            aCompTri.PixelsInside(aVectorToFillWithInsidePixels);

            //----------- index of unkown, finds the pixels
            std::vector<int> aVecInd = {2 * aIndicesOfTriKnots.x(), 2 * aIndicesOfTriKnots.x() + 1,
                                        2 * aIndicesOfTriKnots.y(), 2 * aIndicesOfTriKnots.y() + 1,
                                        2 * aIndicesOfTriKnots.z(), 2 * aIndicesOfTriKnots.z() + 1};

            tHomot2d aCurTrPointA(cPt2dr(aVCur(aVecInd.at(0)), aVCur(aVecInd.at(1))), 0); // current homothety translation 1st point of triangle
            tHomot2d aCurTrPointB(cPt2dr(aVCur(aVecInd.at(2)), aVCur(aVecInd.at(3))), 0); // current homothety translation 2nd point of triangle
            tHomot2d aCurTrPointC(cPt2dr(aVCur(aVecInd.at(4)), aVCur(aVecInd.at(5))), 0); // current homothety translation 3rd point of triangle

            size_t aNumberOfInsidePixels = aVectorToFillWithInsidePixels.size();
            // Loop over all pixels inside triangle
            // size_t is necessary as there can be a lot of pixels in triangles
            for (size_t aFilledPixel = 0; aFilledPixel < aNumberOfInsidePixels; aFilledPixel++)
            {
                const cPtInsideTriangles aPixInsideTriangle = cPtInsideTriangles(aCompTri, aVectorToFillWithInsidePixels, aFilledPixel, *mDImPre);
                // prepare for barycenter translation formula by filling aVObs with different coordinates
                FormalInterpBarycenter_SetObs(aVObs, 0, aPixInsideTriangle);

                // image of a point in triangle by current homothety
                const cPt2dr aTranslatedFilledPoint = ApplyBarycenterTranslationFormulaToFilledPixel(aCurTrPointA, aCurTrPointB, aCurTrPointC, aVObs);
                if (mDImPost->InsideBL(aTranslatedFilledPoint)) // avoid errors
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
                    aNbOut++; // Count number of pixels translated outside post image
            }
        }

        // Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();

        if (mShow)
            StdOut() << aSomDif / aTotalNumberOfInsidePixels << ", " << aNbOut << std::endl;
    }

    void cAppli_cTriangleDeformation::DoOneIteration(bool aIsLastIter)
    {
        LoopOverTrianglesAndUpdateParameters(aIsLastIter); // Iterate over triangles and solve system

        // Show final translation results
        if (aIsLastIter && mShow)
        {
            tDensevect aVFinalSol = mSys->CurGlobSol();

            for (size_t aLastTrCoordinate=0; aLastTrCoordinate < mDelTri.NbPts(); aLastTrCoordinate++)
            {
                tHomot2d aLastTrPoint(cPt2dr(aVFinalSol(2*aLastTrCoordinate), aVFinalSol(2*aLastTrCoordinate + 1)), 0); // final homothety translation for points in triangulation

                StdOut() << "The un-translated point has the following coordinates : " << mDelTri.KthPts(aLastTrCoordinate) << ". The final translation on x-axis of this point is : " 
                         << aLastTrPoint.Tr().x() << " on the x-axis and " << aLastTrPoint.Tr().y() << " for y-axis." << std::endl;                         
            }
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
                     << " NbOut" << std::endl;

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

            std::vector<cPt2dr> AllPixCoordinatesInImOut;
            for (int aPixXCoord=0; aPixXCoord < mDImOut->Sz().x(); aPixXCoord++)
            {
                for (int aPixYCoord=0; aPixYCoord < mDImOut->Sz().y(); aPixYCoord++)
                    AllPixCoordinatesInImOut.push_back(cPt2dr(aPixXCoord, aPixYCoord));
            }

            for (const auto &aDisplacedPix : mFinalTranslatedPixelCoords)
            {
                cPt2dr aDisplacedPoint = cPt2dr(aDisplacedPix.x(), aDisplacedPix.y());
                if (std::find(AllPixCoordinatesInImOut.begin(), AllPixCoordinatesInImOut.end(), 
                              aDisplacedPoint) != AllPixCoordinatesInImOut.end())
                    mDImOut->SetV(cPt2di(aDisplacedPoint.x(), aDisplacedPoint.y()), 
                                  mDImPost->DefGetVBL(aDisplacedPoint, 0));
                else
                    mDImOut->SetV(cPt2di(aDisplacedPoint.x(), aDisplacedPoint.y()), mDImPost->DefGetVBL(aDisplacedPoint, 0));
            }

            mDImOut->ToFile("DisplacedPixels.tif");
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
