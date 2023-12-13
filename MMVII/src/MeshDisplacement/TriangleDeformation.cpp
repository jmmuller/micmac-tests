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

        void OneIterationFitModele();

    private:
        // ==    Mandatory args ====
        std::string mNamePreImage; //< Name of given pre-image.
        std::string mNamePostImage; //< Name of given post-image.

        // ==    Optionnal args ====
        bool mShow; // print result, export image ...

        // ==    Internal variables ====
        std::vector<cPt2dr> mVPtsMod;  ///<  points sampling the model.
        std::vector<double> mValueMod; ///<  values of points in mVPtsMod.

        cPt2di mSzImPre; ///<  size of image.
        tIm mImPre;      ///<  memory representation of the image.
        tDIm &mDImPre;   ///<  memory representation of the image.

        cPt2di mSzImPost; ///<  size of image.
        tIm mImPost;      ///<  memory representation of the image.
        tDIm &mDImPost;   ///<  memory representation of the image.

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
               << Arg2007(mNamePostImage, "Name of post-image file.", {{eTA2007::FileImage}, {eTA2007::FileDirProj}});
    }

    cCollecSpecArg2007 &cAppli_cTriangleDeformation::ArgOpt(cCollecSpecArg2007 &anArgOpt)
    {

        return anArgOpt
               << AOpt2007(mShow, "Show", "Whether to print result", {eTA2007::HDV});
    }

    void cAppli_cTriangleDeformation::OneIterationFitModele()
    {
        //----------- index of unkown, basic here because the unknown are the same for each equation
        std::vector<int> aVecInd{0, 1, 2};
        //----------- allocate vec of obs :
        std::vector<double> aVObs(4, 0.0);

        //----------- extract current parameters
        cDenseVect<double> aVCur = mSys->CurGlobSol();
        cHomot2D<tREAL8> aCurHomM2I(cPt2dr(aVCur(3), aVCur(4)), aVCur(2));    // current homothety
        // double aCurScR = aVCur(0);                                         // current scale on radiometry
        // double aCurTrR = aVCur(0);                                         // current translation on radiometry

        //----------- declaration of indicator of convergence
        double aSomDif = 0; // sum of difference between model and image
        double aSomMod = 0; // sum of value of model, to normalize the difference
        double aNbOut = 0;  // number of points out of image

        // Parse all the point to add the observations on each point
        for (size_t aKPt = 0; aKPt < mVPtsMod.size(); aKPt++)
        {
            cPt2dr aPMod = mVPtsMod[aKPt];         // point of model
            cPt2dr aPIm = aCurHomM2I.Value(aPMod); // image of aPMod by current homotethy
            if (mDImPre.InsideBL(aPIm))             // avoid error
            {
                // put observations in vectors
                //  observations on image and point-image
                // Set_FormalBilinIm2D_Obs(aVObs, 0, aPIm, mDIm);
                Set_TriangleDeform_Obs(aVObs, 0, aPIm, mDImPre, mDImPost);

                //  observation point model and value model
                aVObs[6] = aPMod.x();
                aVObs[7] = aPMod.y();
                aVObs[8] = mValueMod[aKPt];

                // Now add observation
                mSys->CalcAndAddObs(mEqHomTri, aVecInd, aVObs);

                // compute indicator
                // double aDif = mDImIn.GetVBL(aPIm) - (aCurTrR + mValueMod[aKPt]); // residual
                double aDif = mDImPre.GetVBL(aPIm) - (mValueMod[aKPt]); // residual
                aSomMod += mValueMod[aKPt];
                aSomDif += std::abs(aDif);
            }
            else
                aNbOut++;
        }

        // Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();

        if (mShow)
            StdOut() << " Dif= " << aSomDif / aSomMod << " NbOut=" << aNbOut << std::endl;
        /*
        //  If we are at end, check that the model is equal (up to numerical accuracy)  to the target
        if (IsLast)
        {
            MMVII_INTERNAL_ASSERT_bench(aNbOut == 0, "Gauss compos");
            MMVII_INTERNAL_ASSERT_bench((aSomDif / aSomMod) < 1e-10, "Gauss compos");
        }
        */
    }

    int cAppli_cTriangleDeformation::Exe()
    {
        mImPre = tIm::FromFile(mNamePreImage);
        mImPost = tIm::FromFile(mNamePostImage);

        // mDImIn = &mImIn.DIm();
        // mSz = mDImIn->Sz();

        // cAppli_cTriangleDeformation aTriangleDeform(); // not thinking object this isn't great.

        const int aNbIters = 8;
        for (int aIter = 0; aIter < aNbIters; aIter++)
            OneIterationFitModele();

        return EXIT_SUCCESS;
    }

    /********************************************/
    /*              ::MMVII                     */
    /********************************************/
    /*
    void ComputeTriangleDeformation()
    {
        cPt2di SizeOfImage = {500, 300};
        cAppli_cTriangleDeformation aTriangleDeform(false, SizeOfImage);
        int aNbIters = 8;
        for (int aIter = 0; aIter < aNbIters; aIter++)
            aTriangleDeform.OneIterationFitModele(aIter == (aNbIters - 1));
    }
    */

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