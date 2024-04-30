#ifndef _TRIANGLEDEFORMATIONTRANSLATION_H_
#define _TRIANGLEDEFORMATIONTRANSLATION_H_

#include "MMVII_Geom2D.h"
#include "MMVII_PhgrDist.h"

#include "TriangleDeformation.h"

using namespace NS_SymbolicDerivative;

namespace MMVII
{

    /******************************************/
    /*                                        */
    /*    cTriangleDeformationTranslation     */
    /*                                        */
    /******************************************/

    class cAppli_cTriangleDeformationTranslation : public cAppli_cTriangleDeformation
    {
    public:
        typedef cIm2D<tREAL8> tIm;
        typedef cDataIm2D<tREAL8> tDIm;
        typedef cTriangle<tREAL8, 2> tTri2dr;
        typedef cDenseVect<double> tDenseVect;
        typedef std::vector<int> tIntVect;
        typedef std::vector<double> tDoubleVect;

        cAppli_cTriangleDeformationTranslation(const std::vector<std::string> &aVArgs,
                                               const cSpecMMVII_Appli &aSpec);
        ~cAppli_cTriangleDeformationTranslation();

        int Exe() override;
        cCollecSpecArg2007 &ArgObl(cCollecSpecArg2007 &anArgObl) override;
        cCollecSpecArg2007 &ArgOpt(cCollecSpecArg2007 &anArgOpt) override;

        // Iterate of triangles and inside pixels
        void DoOneIterationTranslation(const int aIterNumber);
        // Loops over all triangles and solves system to update parameters at end of iteration
        void LoopOverTrianglesAndUpdateParametersTranslation(const int aIterNumber);
        // Generate displacement maps of last solution
        void GenerateDisplacementMaps(const tDenseVect &aVFinalSol, const int aIterNumber);
        // Generates Displacement maps and coordinates of points in triangulation at last iteration
        void GenerateDisplacementMapsAndLastTranslatedPoints(const int aIterNumber);
        // Initialise problem after user has input information
        void InitialisationAfterExeTranslation();
        // Fill displacement maps and output image
        void FillDisplacementMaps(const cPtInsideTriangles &aLastPixInsideTriangle,
                                  const cPt2dr &aLastTranslatedFilledPoint);
        // Apply  barycentration translation formula to last translation values in optimisation process
        cPt2dr ApplyLastBarycenterTranslationFormulaToInsidePixel(const cPt2dr &aLastTranslationPointA,
                                                                  const cPt2dr &aLastTranslationPointB,
                                                                  const cPt2dr &aLastTranslationPointC,
                                                                  const cPtInsideTriangles &aLastPixInsideTriangle);

    private:
        // ==  Mandatory args ====

        std::string mNamePreImage;   // Name of given pre-image
        std::string mNamePostImage;  // Name of given post-image
        int mNumberPointsToGenerate; // number of generated points
        int mNumberOfScales;         // number of iterations in optimisation process

        // ==  Optionnal args ====

        int mRandomUniformLawUpperBoundLines;         // Uniform law generates random coordinates in interval [0, mRandomUniformLawUpperBoundLines [
        int mRandomUniformLawUpperBoundCols;          // Uniform law generates random coordinates in interval [0, mRandomUniformLawUpperBoundCols [
        bool mShow;                                   // Print result, export image ...
        bool mComputeAvgMax;                          // Compute average and maximum pixel value of difference image between pre and post images
        bool mUseMultiScaleApproach;                  // Apply multi-scale approach or not
        int mSigmaGaussFilterStep;                    // Decreasing step of sigma value during iterations
        bool mGenerateDisplacementImage;              // Generate image with displaced pixels
        int mNumberOfIterGaussFilter;                 // Number of iterations to be done in Gauss filter algorithm
        int mNumberOfEndIterations;                   // Number of iterations to do while using original image in multi-scale approach
        bool mDisplayLastTranslatedPointsCoordinates; // Whether to display the final coordinates of the translated points

        // ==  Internal variables ====

        cPt2di mSzImPre; //  size of image
        tIm mImPre;      //  memory representation of the image
        tDIm *mDImPre;   //  memory representation of the image

        cPt2di mSzImPost; //  size of image
        tIm mImPost;      //  memory representation of the image
        tDIm *mDImPost;   //  memory representation of the image

        cPt2di mSzImOut; //  size of image
        tIm mImOut;      //  memory representation of the image
        tDIm *mDImOut;   //  memory representation of the image

        cPt2di mSzImDepX; //  size of image
        tIm mImDepX;      //  memory representation of the image
        tDIm *mDImDepX;   //  memory representation of the image

        cPt2di mSzImDepY; //  size of image
        tIm mImDepY;      //  memory representation of the image
        tDIm *mDImDepY;   //  memory representation of the image

        double mSigmaGaussFilter; // Value of sigma in gauss filter
        bool mIsLastIters;        // Determines whether optimisation process is at last iters to optimise on original image

        std::vector<cPt2dr> mVectorPts;   // A vector containing a set of points
        cTriangulation2D<tREAL8> mDelTri; // A Delaunay triangle

        cResolSysNonLinear<tREAL8> *mSys;       // Non Linear Sys for solving problem
        cCalculator<double> *mEqTranslationTri; // calculator giving access to values and derivatives
    };
}

#endif // _TRIANGLEDEFORMATIONTRANSLATION_H_