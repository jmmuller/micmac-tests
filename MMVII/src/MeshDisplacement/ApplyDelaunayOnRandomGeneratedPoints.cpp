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

		void ApplyAndSaveDelaunayTriangulationOnPoints(const std::vector<cPt2dr> &aVP, const bool shift_points_of_triangles);
		void GeneratePointsForDelaunay(const int aNbCol, const int aNbLin, const bool shift_points_of_triangles);
		void ConstructUniformRandomVector(std::vector<cPt2dr> &aVP);
		// void GetCoordinatesInsideTriangles(const cTriangulation2D<tREAL8> & aDelaunayTriangle);
		void ComputeMinimumDistanceToCircle(const std::vector<cPt2dr> & aVP, const cTriangle<tREAL8, 2> & aTri);
		void ComputeShiftedTrianglesAndPoints(cTriangle<tREAL8, 2> & aTri, std::vector<cPt2dr> & ShiftedTriangleCoordinates,
											  const double uniform_random_red_channel, const double uniform_random_green_channel,
											  const double uniform_random_blue_channel);

	private:
		// ==   Mandatory args ====
		std::string mNamePlyFile;
		int mNumberPointsToGenerate;	 // number of generated points
		int mRandomUniformLawUpperBound; // Uniform law generate numbers from [0, mRandomUniformLawUpperBound [
		std::string mNameInputImage;

		// ==   Optionnal args ====

		bool mPlyFileisBinary;
		bool mShiftTriangles;
		std::string mNameModifedTrianglePlyFile;

		// ==    Internal variables ====
		tIm mImIn;    ///<  memory representation of the image
        tDIm *mDImIn; ///<  memory representation of the image
        cPt2di mSz;
	};

	cAppli_RandomGeneratedDelaunay::cAppli_RandomGeneratedDelaunay(const std::vector<std::string> &aVArgs,
																   const cSpecMMVII_Appli &aSpec) : cMMVII_Appli(aVArgs, aSpec),
																									mPlyFileisBinary(false),
																									mShiftTriangles(true),
																									mImIn(cPt2di(1, 1)),
																									mDImIn(nullptr)
	{
	}

	cCollecSpecArg2007 &cAppli_RandomGeneratedDelaunay::ArgObl(cCollecSpecArg2007 &anArgObl)
	{
		return anArgObl
			   << Arg2007(mNameInputImage, "Name of input image file.", {{eTA2007::FileImage}, {eTA2007::FileDirProj}})
			   << Arg2007(mNamePlyFile, "Name of main triangulation file to save in .ply format.", {{eTA2007::FileCloud}})
			   << Arg2007(mNumberPointsToGenerate, "Number of points you want to generate for triangulation.")
			   << Arg2007(mRandomUniformLawUpperBound, "Maximum value that the uniform law can draw from.");
	}

	cCollecSpecArg2007 &cAppli_RandomGeneratedDelaunay::ArgOpt(cCollecSpecArg2007 &anArgOpt)
	{

		return anArgOpt
				<< AOpt2007(mShiftTriangles, "ShiftTriangles", "Whether to shift points of triangles after application of Delaunay triangulation.", {eTA2007::HDV})
				<< AOpt2007(mPlyFileisBinary, "PlyFileIsBinary", "Whether to save the .ply file binarised or not.", {eTA2007::HDV})
				<< AOpt2007(mNameModifedTrianglePlyFile, "NamePlyFileShiftedTriangles", "Name of .ply file for shifted triangles.", {eTA2007::FileCloud});
	}

	//=========================================================

	/*
	void cAppli_RandomGeneratedDelaunay::GetCoordinatesInsideTriangles(const cTriangulation2D<tREAL8> & aDelaunayTriangle)
	{
		std::vector<cPt2dr> vector_of_barycenters;
		for (size_t aFace = 0; aFace < aDelaunayTriangle.NbFace(); aFace++)
		{
			cTriangle<tREAL8, 2> aTriangulatedTri = aDelaunayTriangle.KthTri(aFace);
			cTriangle2DCompiled<tREAL8> aCompiledTri(aTriangulatedTri);
			std::vector<cPt2di> aVectorToFillwithInsidePixels;
			aCompiledTri.PixelsInside(aVectorToFillwithInsidePixels);
			vector_of_barycenters.push_back(aTriangulatedTri.Barry());
		}
		StdOut() << vector_of_barycenters << std::endl;
	}
	*/

	void cAppli_RandomGeneratedDelaunay::ComputeShiftedTrianglesAndPoints(cTriangle<tREAL8, 2> & aTri,
																		  std::vector<cPt2dr> & ShiftedTriangleCoordinates,
																		  const double uniform_random_red_channel, 
																		  const double uniform_random_green_channel,
																		  const double uniform_random_blue_channel)
	{
		cTriangle2DCompiled aCompTri(aTri);

		// Compute shifts depending on point of triangle
		cPt2dr PercentDiffA = 0.01 * (aTri.Pt(0) - aTri.Pt(2));
		cPt2dr PercentDiffB = 0.015 * (aTri.Pt(1) - aTri.Pt(0));
		cPt2dr PercentDiffC = 0.02 * (aTri.Pt(2) - aTri.Pt(1));

		ShiftedTriangleCoordinates.push_back(aTri.Pt(0) + PercentDiffA);
		ShiftedTriangleCoordinates.push_back(aTri.Pt(1) + PercentDiffB);
		ShiftedTriangleCoordinates.push_back(aTri.Pt(2) + PercentDiffC);
		
		// Get pixels inside each triangle and shift them
		std::vector<cPt2di> aVectorToFillwithInsidePixels;
		aCompTri.PixelsInside(aVectorToFillwithInsidePixels);
		for (long unsigned int filledPixel=0; filledPixel < aVectorToFillwithInsidePixels.size(); filledPixel++)
		{					
			if (filledPixel % 10 == 0)
			{
				cPt2dr aFilledPoint(aVectorToFillwithInsidePixels[filledPixel].x(), aVectorToFillwithInsidePixels[filledPixel].y());
				cPt3dr barycenter_coordinates = aCompTri.CoordBarry(aFilledPoint);
				cPt2dr ShiftedInsidePixels = cPt2dr(aFilledPoint.x() + barycenter_coordinates.x() * PercentDiffA.x() +
													barycenter_coordinates.y() * PercentDiffB.x() + barycenter_coordinates.z() * PercentDiffC.x(), 
													aFilledPoint.y() + barycenter_coordinates.x() * PercentDiffA.y() +
													barycenter_coordinates.y() * PercentDiffB.y() + barycenter_coordinates.z() * PercentDiffC.y());

				StdOut() << aFilledPoint.x() << " " << aFilledPoint.y() << " " << uniform_random_red_channel 
							<< " " << uniform_random_green_channel << " " << uniform_random_blue_channel << std::endl;
				StdOut() << ShiftedInsidePixels.x() << " " << ShiftedInsidePixels.y() << " " 
							<< uniform_random_red_channel << " " << uniform_random_green_channel << " " << uniform_random_blue_channel << std::endl;
			}
		}
	}

	void cAppli_RandomGeneratedDelaunay::ComputeMinimumDistanceToCircle(const std::vector<cPt2dr> & aVP, const cTriangle<tREAL8, 2> & aTri)
	{
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

	void cAppli_RandomGeneratedDelaunay::ApplyAndSaveDelaunayTriangulationOnPoints(const std::vector<cPt2dr> & aVP, const bool shift_points_of_triangles)
	{
		cTriangulation2D<tREAL8> aDelTri(aVP);

		aDelTri.MakeDelaunay();

		std::vector<cPt2dr> ShiftedTriangleCoordinates;

		// Loop over all triangle
		for (size_t aKt = 0; aKt < aDelTri.NbFace(); aKt++)
		{
			cTriangle<tREAL8, 2> aTri = aDelTri.KthTri(aKt);
			ComputeMinimumDistanceToCircle(aVP, aTri);

			// for colouring points in representation
			double uniform_random_red_channel = RandUnif_N(256);
			double uniform_random_green_channel = RandUnif_N(256);
			double uniform_random_blue_channel = RandUnif_N(256);

			if (shift_points_of_triangles)
				ComputeShiftedTrianglesAndPoints(aTri, ShiftedTriangleCoordinates, uniform_random_red_channel, 
												 uniform_random_green_channel, uniform_random_blue_channel);
		}

		cTriangulation2D<tREAL8> aModifiedDelTri(ShiftedTriangleCoordinates);

		aModifiedDelTri.MakeDelaunay();
		// GetCoordinatesInsideTriangles(aDelTri); // Other development

		// Save files to .ply format
		aDelTri.WriteFile(mNamePlyFile, mPlyFileisBinary);
		aModifiedDelTri.WriteFile(mNameModifedTrianglePlyFile, mPlyFileisBinary);
	}

	void cAppli_RandomGeneratedDelaunay::ConstructUniformRandomVector(std::vector<cPt2dr> & aVPts)
	{
		// Generate coordinates from drawing lines and columns of coordinates from a uniform distribution
		for (int aNbPt = 0; aNbPt < mNumberPointsToGenerate; aNbPt++)
		{
			double uniform_random_line = RandUnif_N(mRandomUniformLawUpperBound);
			double uniform_random_col = RandUnif_N(mRandomUniformLawUpperBound);
			cPt2dr aUniformRandomPt(uniform_random_line, uniform_random_col);
			aVPts.push_back(aUniformRandomPt);
		}
	}

	void cAppli_RandomGeneratedDelaunay::GeneratePointsForDelaunay(const int aNbLin, const int aNbCol, const bool shift_points_of_triangles)
	{
		std::vector<cPt2dr> aVPts;

		const int minimum_lin_col = std::min(aNbLin, aNbCol);

		// make sure that values greater than image size can't be drawn from uniform law.
		while (mRandomUniformLawUpperBound >= minimum_lin_col)
		{
			StdOut() << "Maximum value drawn from uniform law needs to be smaller than " << minimum_lin_col << "." << std::endl;
			std::cin >> mRandomUniformLawUpperBound;
		}
		ConstructUniformRandomVector(aVPts);

		ApplyAndSaveDelaunayTriangulationOnPoints(aVPts, shift_points_of_triangles); // Apply Delaunay triangulation on generated points.
	}

	//----------------------------------------

	int cAppli_RandomGeneratedDelaunay::Exe()
	{
		/* 
		MMVII RandomGeneratedDelaunay pair18_im1_720.png OriginalTriangles.ply 20 50 NamePlyFileShiftedTriangles=ShiftedTriangles.ply > OriginalAndShiftedPoints.xyz
		awk NR%2==1 < OriginalAndShiftedPoints.xyz > OriginalPoints.xyz
		awk NR%2==0 < OriginalAndShiftedPoints.xyz > ShiftedPoints.xyz
		*/

		mImIn = tIm::FromFile(mNameInputImage);

		mDImIn = &mImIn.DIm();
		mSz = mDImIn->Sz();

		GeneratePointsForDelaunay(mSz.y(), mSz.x(), mShiftTriangles);

		return EXIT_SUCCESS;
	}

	/* ================================== */
	/*                                    */
	/*               MMVII                */
	/*                                    */
	/* ================================== */

	tMMVII_UnikPApli Alloc_RandomGeneratedDelaunay(const std::vector<std::string> &aVArgs, const cSpecMMVII_Appli &aSpec)
	{
		return tMMVII_UnikPApli(new cAppli_RandomGeneratedDelaunay(aVArgs, aSpec));
	}

	cSpecMMVII_Appli TheSpec_RandomGeneratedDelaunay(
		"RandomGeneratedDelaunay",
		Alloc_RandomGeneratedDelaunay,
		"Generate random points and apply Delaunay triangulation",
		{eApF::ImProc}, // category
		{eApDT::Image}, // input
		{eApDT::Ply},	// output
		__FILE__);

}; // MMVII
