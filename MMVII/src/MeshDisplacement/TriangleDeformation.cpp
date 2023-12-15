#include "cMMVII_Appli.h"
// #include "MMVII_enums.h"
// #include "MMVII_Image2D.h"
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

    
        void ComputeMinimumDistanceToCircle(const cTriangle<tREAL8, 2> & aTri);
        void BuildTrianglesAndApplyDelaunayTriangulation();
        void ConstructUniformRandomVector();
        void GeneratePointsForDelaunay();
        void OneIterationFitModele();

    private:
        // ==    Mandatory args ====
        std::string mNamePreImage; //< Name of given pre-image.
        std::string mNamePostImage; //< Name of given post-image.

        // ==    Optionnal args ====
        bool mShow; // print result, export image ...
        int mNumberPointsToGenerate;	 // number of generated points
		int mRandomUniformLawUpperBound; // Uniform law generate numbers from [0, mRandomUniformLawUpperBound [

        // ==    Internal variables ====
        // std::vector<cPt2dr> mVPtsMod;  ///<  points sampling the model.
        // std::vector<double> mValueMod; ///<  values of points in mVPtsMod.

        cPt2di mSzImPre; ///<  size of image.
        tIm mImPre;      ///<  memory representation of the image.
        tDIm &mDImPre;   ///<  memory representation of the image.

        cPt2di mSzImPost; ///<  size of image.
        tIm mImPost;      ///<  memory representation of the image.
        tDIm &mDImPost;   ///<  memory representation of the image.

        std::vector<cPt2dr> mVectorPts; // A vector containing a set of points
        cTriangulation2D<tREAL8> mDelTri; // A Delaunay triangle

        cResolSysNonLinear<tREAL8> *mSys; ///< Non Linear Sys for solving problem.
        cCalculator<double> *mEqHomTri;   ///< calculator giving access to values and derivatives.
    };

    cAppli_cTriangleDeformation::cAppli_cTriangleDeformation(const std::vector<std::string> &aVArgs,
                                                             const cSpecMMVII_Appli &aSpec) : cMMVII_Appli(aVArgs, aSpec),
                                                                                              mShow(true),
                                                                                              mSzImPre(cPt2di(1, 1)),
                                                                                              mImPre(mSzImPre),
                                                                                              mDImPre(mImPre.DIm()),
                                                                                              mSzImPost(cPt2di(1, 1)),
                                                                                              mImPost(mSzImPost),
                                                                                              mDImPost(mImPost.DIm()),
                                                                                              mVectorPts({cPt2dr(0, 0)}),
                                                                                              mDelTri(mVectorPts),
                                                                                              mSys(nullptr),
                                                                                              mEqHomTri(nullptr)
    {
        cDenseVect<double> aVInit(2);
        aVInit(0) = 0;
        aVInit(1) = 0;
        mSys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqDense, aVInit);
        mEqHomTri = EqDeformTriHomothety(true, 1); // true-> with derivative,  1=sz of buffer
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
               << Arg2007(mNamePostImage, "Name of post-image file.", {{eTA2007::FileImage}, {eTA2007::FileDirProj}})
               << Arg2007(mNumberPointsToGenerate, "Number of points you want to generate for triangulation.")
			   << Arg2007(mRandomUniformLawUpperBound, "Maximum value that the uniform law can draw from.");
    }

    cCollecSpecArg2007 &cAppli_cTriangleDeformation::ArgOpt(cCollecSpecArg2007 &anArgOpt)
    {

        return anArgOpt
               << AOpt2007(mShow, "Show", "Whether to print result", {eTA2007::HDV});
    }

	void cAppli_cTriangleDeformation::ComputeMinimumDistanceToCircle(const cTriangle<tREAL8, 2> & aTri)
	{
			// Compute center circle circum
			const cPt2dr aC = aTri.CenterInscribedCircle();
			// Compute min dist to this circle
			double aMinDist = 1e20;
			for (const auto &aPt : mVectorPts)
				aMinDist = std::min(aMinDist, Norm2(aC - aPt));
			// This  min dist must be (almost) equal to circum-radius
			const double aRadiusCircum = Norm2(aC - aTri.Pt(0));
			const double aDif = std::abs(aRadiusCircum - aMinDist);
			MMVII_INTERNAL_ASSERT_bench(aDif < 1e-5, "Inscribed circle property in Delaunay");
	}

	void cAppli_cTriangleDeformation::BuildTrianglesAndApplyDelaunayTriangulation()
	{
		cTriangulation2D<tREAL8> mDelTri(mVectorPts);

		mDelTri.MakeDelaunay();

		// Loop over all triangle
		for (size_t aKt = 0; aKt < mDelTri.NbFace(); aKt++)
		{
			const cTriangle<tREAL8, 2> aTri = mDelTri.KthTri(aKt);
			ComputeMinimumDistanceToCircle(aTri);
		}
	}

	void cAppli_cTriangleDeformation::ConstructUniformRandomVector()
	{
		// Generate coordinates from drawing lines and columns of coordinates from a uniform distribution
		for (int aNbPt = 0; aNbPt < mNumberPointsToGenerate; aNbPt++)
		{
			const double aUniformRandomLine = RandUnif_N(mRandomUniformLawUpperBound);
			const double aUniformRandomCol = RandUnif_N(mRandomUniformLawUpperBound);
			const cPt2dr aUniformRandomPt(aUniformRandomLine, aUniformRandomCol);
			mVectorPts.push_back(aUniformRandomPt);
		}
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
		ConstructUniformRandomVector();

		BuildTrianglesAndApplyDelaunayTriangulation(); // Apply Delaunay triangulation on generated points.
	}

    void cAppli_cTriangleDeformation::OneIterationFitModele()
    {
        //----------- index of unkown, basic here because the unknown are the same for each equation
        std::vector<int> aVecInd{0, 1, 2, 3, 4, 5};
        //----------- allocate vec of obs :
        std::vector<double> aVObs(5, 0.0);

        //----------- extract current parameters
        cDenseVect<double> aVCur = mSys->CurGlobSol();
        cHomot2D<tREAL8> aCurHomM2I(cPt2dr(aVCur(3), aVCur(4)), aVCur(2));    // current homothety
        // double aCurScR = aVCur(0);                                         // current scale on radiometry
        // double aCurTrR = aVCur(0);                                         // current translation on radiometry

        //----------- declaration of indicator of convergence
        // double aSomDif = 0; // sum of difference between model and image
        // double aSomMod = 0; // sum of value of model, to normalize the difference
        double aNbOut = 0;  // number of points out of image

        // Parse all the point to add the observations on each point
        for (size_t aTr = 0; aTr < mDelTri.NbFace(); aTr++)
        {
            const cTriangle<tREAL8, 2> aTri = mDelTri.KthTri(aTr);
            // cPt2dr aPMod = mVPtsMod[aKPt];         // point of model
            // cPt2dr aPIm = aCurHomM2I.Value(aPMod); // image of aPMod by current homotethy
            // put observations in vectors
            //  observations on image and point-image
            // Set_FormalBilinIm2D_Obs(aVObs, 0, aPIm, mDIm);
            const cTriangle2DCompiled aCompTri(aTri);
            std::vector<cPt2di> aVectorToFillWithInsidePixels;
		    aCompTri.PixelsInside(aVectorToFillWithInsidePixels);

            for (int aKnot=0; aKnot < 3; aKnot++)
            {
                if (mDImPre.InsideBL(aTri.Pt(aKnot)))             // avoid error
                {
                    FormalBilinTri_SetObs(aVObs, 0, aTri.Pt(aKnot), mDImPre);

                    // Now add observation
                    mSys->CalcAndAddObs(mEqHomTri, aVecInd, aVObs);
/*
                    double aDif = mDImPre.GetVBL(aTri.Pt(aKnot)) - (mValueMod[aKPt]); // residual
                    aSomMod += mValueMod[aKPt];
                    aSomDif += std::abs(aDif);
                    */
                }
                else
                    aNbOut++;
                // Set_TriangleDeform_Obs(aVObs, 0, aPIm, mDImPre, mDImPost);

                //  observation point model and value model
                /*
                aVObs[6] = aPMod.x();
                aVObs[7] = aPMod.y();
                aVObs[8] = mValueMod[aKPt];
                */

                // compute indicator
                // double aDif = mDImIn.GetVBL(aPIm) - (aCurTrR + mValueMod[aKPt]); // residual
            }
            
        }

        // Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();
        /*
        if (mShow)
            StdOut() << " Dif= " << aSomDif / aSomMod << " NbOut=" << aNbOut << std::endl;
            */
    }

    int cAppli_cTriangleDeformation::Exe()
    {
        mImPre = tIm::FromFile(mNamePreImage);
        mImPost = tIm::FromFile(mNamePostImage);

        // mDImIn = &mImIn.DIm();
        // mSz = mDImIn->Sz();

        const int aNbIters = 8;
        for (int aIter = 0; aIter < aNbIters; aIter++)
            OneIterationFitModele();

        return EXIT_SUCCESS;
    }

    /********************************************/
    /*              ::MMVII                     */
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