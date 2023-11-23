#include "MMVII_Geom2D.h"
#include "MMVII_nums.h"

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
	public :
	
    cAppli_RandomGeneratedDelaunay();

    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

	void  ApplyAndSaveDelaunayTriangulationOnPoints();
	void GeneratePointsForDelaunay;

    private :

	// ==   Mandatory args ====
	std::string mNamePlyFile;
	int mNumberPointsToGenerate;
	int mRandomUniformLawUpperBound;

	// ==   Optionnal args ====

	int mLimitNumberOfPoints;

	// ==    Internal variables ====

};


cAppli_RandomGeneratedDelaunay::cAppli_RandomGeneratedDelaunay (const std::vector<std::string> & aVArgs,
																const cSpecMMVII_Appli & aSpec)
{
	cMMVII_Appli (aVArgs, aSpec);
	int mRandomUniformLawUpperBound(1000);
}


cCollecSpecArg2007 & cAppli_RandomGeneratedDelaunay::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return anArgObl
             << Arg2007(mNamePlyFile, "Name of file to save in .ply format.",{{eTA2007::FileName},{eTA2007::FileDirProj}})
			 << Arg2007(mNumberPointsToGenerate, "Number of points you want to generate for triangulation.")
			 << Arg2007(mRandomUniformLawUpperBound, "Maximum value that the uniform law can draw from.")
           ;
}


cCollecSpecArg2007 & cAppli_RandomGeneratedDelaunay::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return    anArgOpt
	   << AOpt2007(mLimitNumberOfPoints, "PointsLimit", "Maximum possible number of points to generate for triangulation.", {eTA2007::HDV})
    ;
}


//=========================================================


void ApplyAndSaveDelaunayTriangulationOnPoints(const std::vector<cPt2dr> & aVP, const int aVectorSizeLimit, 
											   const std::string & aFileName, const bool aisBinary)
{
	cTriangulation2D<tREAL8> aDelTri(aVP);
	aDelTri.MakeDelaunay();
   
	if (aVP.size()<=aVectorSizeLimit)
    {
		// Parse all triangle
		for (size_t aKt=0 ; aKt < aDelTri.NbFace() ; aKt++)
		{
			cTriangle<tREAL8,2> aTri = aDelTri.KthTri(aKt);
			// Compute center circle circum
			cPt2dr aC = aTri.CenterInscribedCircle() ;
			// Compute min dist to this circle
			double aMinDist = 1e20;
			for (const auto & aPt : aVP)
				aMinDist = std::min(aMinDist, Norm2(aC-aPt));
			// This  min dist must be (almost) equal to circum-radius
			double aRadiusCircum = Norm2(aC-aTri.Pt(0));
			double aDif = std::abs(aRadiusCircum -aMinDist);
			MMVII_INTERNAL_ASSERT_bench(aDif<1e-5, "Inscribed circle property in delaunay");
	   }
	}
	aDelTri.WriteFile(aFileName, aisBinary);
}


void GeneratePointsForDelaunay(const int aNbPts, const int aUpperBoundUniform)
{
   std::vector<cPt2dr> aVP;
   for (int aK=0 ; aK<aNbPts ; aK++)
   {
        double coords = RandUnif_N(aUpperBoundUniform);
		aVP.push_back(coords);
   } 
   ApplyAndSaveDelaunayTriangulationOnPoints(aVP); // Apply Delaunay on generated points.
}

//----------------------------------------


int cAppli_RandomGeneratedDelaunay::Exe()
{
	GeneratePointsForDelaunay();
}



/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */


cSpecMMVII_Appli  TheSpec_RandomGeneratedDelaunay
(
     "RandomGeneratedDelaunay",
      RandomGeneratedDelaunay,
      "Generate random points and apply Delaunay triangulation",
      {eApF::ComputationalGeom}, // category
      {eApDT::VectorPoint}, // input 
      {eApDT::VectorPoint}, // output
      __FILE__
);

}; //MMVII