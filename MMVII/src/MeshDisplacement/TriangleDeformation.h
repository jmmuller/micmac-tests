#ifndef _TRIANGLEDEFORMATION_H_
#define _TRIANGLEDEFORMATION_H_

#include "MMVII_Geom2D.h"
#include "MMVII_PhgrDist.h"

using namespace NS_SymbolicDerivative;

namespace MMVII
{
    /****************************************/
    /*                                      */
    /*          cPtInsideTriangles          */
    /*                                      */
    /****************************************/
    class cPtInsideTriangles
    {
    public:
        cPtInsideTriangles(const cTriangle2DCompiled<tREAL8> &aCompTri,              // a compiled triangle.
                           const std::vector<cPt2di> &aVectorFilledwithInsidePixels, // vector containing pixels inside triangles.
                           const size_t aFilledPixel,                                // a counter that is looping over pixels in triangles.
                           const cDataIm2D<tREAL8> &aDIm);                           // image.
        cPt3dr GetBarycenterCoordinates() const;
        cPt2dr GetCartesianCoordinates() const;
        tREAL8 GetPixelValue() const;

    private:
        cPt3dr mBarycenterCoordinatesOfPixel; // Barycentric coordinates of pixel.
        cPt2dr mFilledIndices;                // 2D cartesian coordinates of pixel.
        tREAL8 mValueOfPixel;                 // Intensity in image at pixel.
    };

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
        typedef cTriangle<tREAL8, 2> tTri2dr;
        typedef cHomot2D<tREAL8> tHomot2d;
        typedef cDenseVect<double> tDensevect;

        cAppli_cTriangleDeformation(const std::vector<std::string> &aVArgs,
                                    const cSpecMMVII_Appli &aSpec);
        ~cAppli_cTriangleDeformation();

        int Exe() override;
        cCollecSpecArg2007 &ArgObl(cCollecSpecArg2007 &anArgObl) override;
        cCollecSpecArg2007 &ArgOpt(cCollecSpecArg2007 &anArgOpt) override;

        void ConstructUniformRandomVectorAndApplyDelaunay();
        void GeneratePointsForDelaunay();
        void SubtractPrePostImageAndComputeAvgAndMax();
        void DoOneIteration(const bool aIsLastIter);
        void LoopOverTrianglesAndUpdateParameters(const bool aIsLastIter);
        void InitialisationAfterExe();
        cPt2dr ApplyBarycenterTranslationFormulaToFilledPixel(const tHomot2d &aCurrentTranslationPointA, const tHomot2d &aCurrentTranslationPointB,
                                                              const tHomot2d &aCurrentTranslationPointC, std::vector<double> &aVObs);

    private:
        // ==  Mandatory args ====
        std::string mNamePreImage;           //< Name of given pre-image.
        std::string mNamePostImage;          //< Name of given post-image.
        int mNumberPointsToGenerate;         // number of generated points
        int mNumberOfOptimisationIterations; // number of iterations in optimisation process

        // ==  Optionnal args ====
        int mRandomUniformLawUpperBoundLines; // Uniform law generates random coordinates in interval [0, mRandomUniformLawUpperBoundLines [
        int mRandomUniformLawUpperBoundCols;  // Uniform law generates random coordinates in interval [0, mRandomUniformLawUpperBoundCols [
        bool mShow;                           // print result, export image ...
        bool mGenerateDisplacementImage;      // Generate image with displaced pixels

        // ==  Internal variables ====

        cPt2di mSzImPre; ///<  size of image
        tIm mImPre;      ///<  memory representation of the image
        tDIm *mDImPre;   ///<  memory representation of the image

        cPt2di mSzImPost; ///<  size of image
        tIm mImPost;      ///<  memory representation of the image
        tDIm *mDImPost;   ///<  memory representation of the image

        cPt2di mSzImOut; ///<  size of image
        tIm mImOut;      ///<  memory representation of the image
        tDIm *mDImOut;   ///<  memory representation of the image

        cPt2di mSzImDiff; ///<  size of image
        tIm mImDiff;      ///<  memory representation of the image
        tDIm *mDImDiff;   ///<  memory representation of the image

        cPt2di mSzImDepX; ///<  size of image
        tIm mImDepX;      ///<  memory representation of the image
        tDIm *mDImDepX;   ///<  memory representation of the image

        cPt2di mSzImDepY; ///<  size of image
        tIm mImDepY;      ///<  memory representation of the image
        tDIm *mDImDepY;   ///<  memory representation of the image

        std::vector<cPt2dr> mVectorPts;                  // A vector containing a set of points
        cTriangulation2D<tREAL8> mDelTri;                // A Delaunay triangle

        cResolSysNonLinear<tREAL8> *mSys; ///< Non Linear Sys for solving problem
        cCalculator<double> *mEqHomTri;   ///< calculator giving access to values and derivatives
    };
}

#endif // _TRIANGLEDEFORMATION_H_