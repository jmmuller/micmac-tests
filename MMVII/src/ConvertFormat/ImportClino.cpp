#include "MMVII_PCSens.h"
#include "MMVII_MMV1Compat.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_BundleAdj.h"
#include "MMVII_Stringifier.h"

/**
   \file cConvCalib.cpp  testgit

   \brief file for conversion between calibration (change format, change model) and tests
*/


namespace MMVII
{

   /* ********************************************************** */
   /*                                                            */
   /*                 cAppli_CERN_ImportClino                    */
   /*                                                            */
   /* ********************************************************** */

class cAppli_CERN_ImportClino : public cMMVII_Appli
{
     public :
        cAppli_CERN_ImportClino(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
     private :

        void MakeOneDir(const std::string & aDir,cMMVII_Ofs &) const;
	cPhotogrammetricProject  mPhProj;

	// Mandatory Arg
	std::string               mNameRes;
	std::vector<std::string>  Samples() const override;

	std::vector<std::string>  mNamesClino;
};

cAppli_CERN_ImportClino::cAppli_CERN_ImportClino(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
   cMMVII_Appli  (aVArgs,aSpec),
   mPhProj       (*this),
   mNamesClino   {"A1","B1","B2","A2"}
{
}

cCollecSpecArg2007 & cAppli_CERN_ImportClino::ArgObl(cCollecSpecArg2007 & anArgObl) 
{
    return anArgObl
	      <<  Arg2007(mNameRes ,"Folder of external project where the calib is to be imported")
           ;
}

cCollecSpecArg2007 & cAppli_CERN_ImportClino::ArgOpt(cCollecSpecArg2007 & anArgFac) 
{
    
    return anArgFac
       //  << AOpt2007(mNameGCP,"NameGCP","Name of GCP set")
       //  << AOpt2007(mNbDigName,"NbDigName","Number of digit for name, if fixed size required (only if int)")
       //  << AOpt2007(mL0,"NumL0","Num of first line to read",{eTA2007::HDV})
       //  << AOpt2007(mLLast,"NumLast","Num of last line to read (-1 if at end of file)",{eTA2007::HDV})
       //  << AOpt2007(mPatternTransfo,"PatName","Pattern for transforming name (first sub-expr)")
    ;
}

//template<class Type> inline Type GetV(std::istringstream & iss);


void cAppli_CERN_ImportClino::MakeOneDir(const std::string & aDir,cMMVII_Ofs & anOFS) const
{
    std::string aNameF = aDir + StringDirSeparator() +  "ClinoValue.json";
    std::ifstream infile(aNameF);

    std::vector<std::string>  aVFileIm =  GetFilesFromDir(aDir+StringDirSeparator(),AllocRegex("043.*"));

    MMVII_INTERNAL_ASSERT_tiny(aVFileIm.size()==1,"cAppli_CERN_ImportClino : bad size for image pattern match");
    anOFS.Ofs() << aVFileIm.at(0) ;
    StdOut() << "DDD " << aDir << " " << aVFileIm << "\n";

    std::string line;
    // mNbLineRead = 0;
    // int aNumL = 0;
    cCarLookUpTable aLUT;
    aLUT.InitIdGlob();
    aLUT.Init("[],",' ');
    while (std::getline(infile, line))
    {
        line = aLUT.Translate(line);
        //  StdOut() << "DDD=" << aDir << " [" << line << "]\n";
	std::istringstream iss(line);

	for (const auto & aNameClino : mNamesClino)
	{
		tREAL8 aAvg = GetV<tREAL8>(iss);
		tREAL8 aStdDev = GetV<tREAL8>(iss);
		StdOut() << "   * " << aNameClino  << " " << aAvg << " " << aStdDev << "\n";
		anOFS.Ofs() << " " << aNameClino  << " " << aAvg << " " << aStdDev ;
	}
    }
    anOFS.Ofs() << std::endl;
}


int cAppli_CERN_ImportClino::Exe()
{
    mPhProj.FinishInit();

    tNameSelector   aSelec = AllocRegex("Calibration_Clino_.*");
    std::vector<std::string>   aLD = GetSubDirFromDir("./",aSelec);
    std::sort(aLD.begin(),aLD.end());

    cMMVII_Ofs anOFS(mNameRes,eFileModeOut::CreateText);
    for (const auto & aDir : aLD)
        MakeOneDir(aDir,anOFS);


    return EXIT_SUCCESS;
}


std::vector<std::string>  cAppli_CERN_ImportClino::Samples() const
{
	return {"MMVII V2ImportCalib ../../Pannel/ BA_725 CalibInit725"};
}


tMMVII_UnikPApli Alloc_CERN_ImportClino(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CERN_ImportClino(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_CERN_ImportClino
(
     "CERN_ImportClino",
      Alloc_CERN_ImportClino,
      "A temporary command to arrange clino format",
      {eApF::Ori},
      {eApDT::Ori},
      {eApDT::Ori},
      __FILE__
);
/*
*/


}; // MMVII

