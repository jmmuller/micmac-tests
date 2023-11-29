#include <string>
#include <vector>

#include "cMMVII_Appli.h"
#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"
#include "MMVII_nums.h"
#include "MMVII_enums.h"

/**
   \file ApplyDelaunayOnRandomGeneratedPoints.cpp

   \brief file for generating random points distributed uniformely
   and applying 2D Delaunay triangulation.

**/

namespace MMVII
{

	/* ==================================================== */
	/*                                                      */
	/*          cAppli_RandomGeneratedDelaunay              */
	/*                                                      */
	/* ==================================================== */

	class cAppli_RandomGeneratedDelaunay : public cMMVII_Appli
	{
	public:

		typedef cIm2D<tREAL4> tIm;
        typedef cDataIm2D<tREAL4> tDIm;

		cAppli_RandomGeneratedDelaunay(const std::vector<std::string> &aVArgs,
									   const cSpecMMVII_Appli &aSpec);

		int Exe() override;
		cCollecSpecArg2007 &ArgObl(cCollecSpecArg2007 &anArgObl) override;
		cCollecSpecArg2007 &ArgOpt(cCollecSpecArg2007 &anArgOpt) override;

		void ApplyAndSaveDelaunayTriangulationOnPoints(const std::vector<cPt2dr> &aVP);
		void GeneratePointsForDelaunay(const int aNbCol, const int aNbLin);
		void ConstructUniformRandomVector(std::vector<cPt2dr> & aVP);

	private:
		// ==   Mandatory args ====
		std::string mNamePlyFile;
		int mNumberPointsToGenerate; // number of generated points
		int mRandomUniformLawUpperBound; // Uniform law generate numbers from [0, mRandomUniformLawUpperBound [
		std::string mNameInputImage;		

		// ==   Optionnal args ====

		// long unsigned int mLimitNumberOfPoints;
		bool mPlyFileisBinary;

		// ==    Internal variables ====
		tIm mImIn;    ///<  memory representation of the image
        tDIm *mDImIn; ///<  memory representation of the image
        cPt2di mSz;
        // tIm mImOut;    ///<  memory representation of the image
        // tDIm *mDImOut; ///<  memory representation of the image
	};

	cAppli_RandomGeneratedDelaunay::cAppli_RandomGeneratedDelaunay(const std::vector<std::string> &aVArgs,
        const cSpecMMVII_Appli &aSpec) :
             cMMVII_Appli(aVArgs, aSpec),
			 mPlyFileisBinary(false),
			 mImIn(cPt2di(1, 1)),
             mDImIn(nullptr)
	{
	}

	cCollecSpecArg2007 &cAppli_RandomGeneratedDelaunay::ArgObl(cCollecSpecArg2007 &anArgObl)
	{
		return anArgObl
				<< Arg2007(mNameInputImage, "Name of input image file.", {{eTA2007::FileImage}, {eTA2007::FileDirProj}})		
                << Arg2007(mNamePlyFile, "Name of file to save in .ply format.", {{eTA2007::FileCloud}})
			   	<< Arg2007(mNumberPointsToGenerate, "Number of points you want to generate for triangulation.")
			   	<< Arg2007(mRandomUniformLawUpperBound, "Maximum value that the uniform law can draw from.");
	}

	cCollecSpecArg2007 &cAppli_RandomGeneratedDelaunay::ArgOpt(cCollecSpecArg2007 &anArgOpt)
	{

		return anArgOpt
				// << AOpt2007(mLimitNumberOfPoints, "PointsLimit", "Maximum possible number of points to generate for triangulation.", {eTA2007::HDV})
				<< AOpt2007(mPlyFileisBinary, "PlyFileIsBinary", "Whether to save the .ply file binarised or not.", {eTA2007::HDV});

	}

	//=========================================================

	void cAppli_RandomGeneratedDelaunay::ApplyAndSaveDelaunayTriangulationOnPoints(const std::vector<cPt2dr> &aVP)
	{
		cTriangulation2D<tREAL8> aDelTri(aVP);

		aDelTri.MakeDelaunay();

		// Parse all triangle
		for (size_t aKt = 0; aKt < aDelTri.NbFace(); aKt++)
		{
			cTriangle<tREAL8, 2> aTri = aDelTri.KthTri(aKt);
			// Compute center circle circum
			cPt2dr aC = aTri.CenterInscribedCircle();
			// Compute min dist to this circle
			double aMinDist = 1e20;
			for (const auto &aPt : aVP)
				aMinDist = std::min(aMinDist, Norm2(aC - aPt));
			// This  min dist must be (almost) equal to circum-radius
			double aRadiusCircum = Norm2(aC - aTri.Pt(0));
			double aDif = std::abs(aRadiusCircum - aMinDist);
			MMVII_INTERNAL_ASSERT_bench(aDif < 1e-5, "Inscribed circle property in delaunay");
		}

		aDelTri.WriteFile(mNamePlyFile, mPlyFileisBinary);
	}

	void cAppli_RandomGeneratedDelaunay::ConstructUniformRandomVector(std::vector<cPt2dr> & aVPts)
	{
		for (int aNbPt = 0; aNbPt < mNumberPointsToGenerate; aNbPt++)
		{
			double uniform_random_line = RandUnif_N(mRandomUniformLawUpperBound);
			double uniform_random_col = RandUnif_N(mRandomUniformLawUpperBound);
			cPt2dr aUniformRandomPt(uniform_random_line, uniform_random_col);
			aVPts.push_back(aUniformRandomPt);
		}
	}

	void cAppli_RandomGeneratedDelaunay::GeneratePointsForDelaunay(const int aNbLin, const int aNbCol)
	{
		std::vector<cPt2dr> aVPts;

		const int minimum_lin_col = std::min(aNbLin, aNbCol);

		while (mRandomUniformLawUpperBound >= minimum_lin_col)
		{
			StdOut() << "Maximum value drawn from uniform law needs to be smaller than " << minimum_lin_col << "." << std::endl;
			std::cin >> mRandomUniformLawUpperBound;
		}
		ConstructUniformRandomVector(aVPts);

		ApplyAndSaveDelaunayTriangulationOnPoints(aVPts); // Apply Delaunay triangulation on generated points.
	}

	//----------------------------------------

	int cAppli_RandomGeneratedDelaunay::Exe()
	{
		mImIn = tIm::FromFile(mNameInputImage);
        // cDataFileIm2D aDescFile = cDataFileIm2D::Create(mNameInputImage, false);

        mDImIn = &mImIn.DIm();
        mSz = mDImIn->Sz();

		GeneratePointsForDelaunay(mSz.y(), mSz.x());

		// StdOut() << "hello , size of image = " << mImIn.DIm().Sz().x() << " and " << mDImIn->Sz().y() << "\n";
		return EXIT_SUCCESS;
	}

	/* ==================================================== */
	/*                                                      */
	/*               MMVII                                  */
	/*                                                      */
	/* ==================================================== */

	tMMVII_UnikPApli Alloc_RandomGeneratedDelaunay(const std::vector<std::string> &aVArgs, const cSpecMMVII_Appli &aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_RandomGeneratedDelaunay(aVArgs, aSpec));
    }

	cSpecMMVII_Appli TheSpec_RandomGeneratedDelaunay(
		"RandomGeneratedDelaunay",
		Alloc_RandomGeneratedDelaunay,
		"Generate random points and apply Delaunay triangulation",
		{eApF::ImProc}, // category
		{eApDT::Image},	   // input
		{eApDT::Ply},	   // output
		__FILE__);


}; // MMVII
