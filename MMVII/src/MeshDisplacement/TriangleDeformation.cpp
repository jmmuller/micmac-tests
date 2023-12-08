#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"
#include "MMVII_PhgrDist.h"
#include "../TutoBenchTrianguRSNL/TrianguRSNL.h"

#include "MMVII_TplSymbImage.h"

using namespace NS_SymbolicDerivative;

namespace MMVII
{
    /********************************************/
    /*                                          */
    /*            cTriangleDeformation          */
    /*                                          */
    /********************************************/
    class cTriangleDeformation
    {
    public:
        cTriangleDeformation(bool Show, cPt2di ImageSize);
        ~cTriangleDeformation();
        void OneIterationFitModele(bool aIsLast);

    private:
        bool mShow; // print result, export image ...
        cPt2di             mSzIm;   ///<  size of image
        cIm2D<tREAL8>      mIm;     ///<  image matched
        cDataIm2D<tREAL8> &mDIm;          ///< data-image
        cResolSysNonLinear<tREAL8> *mSys; ///< Non Linear Sys for solving problem
        cCalculator<double> *mEqHomTri;    ///< calculator giving access to values and derivatives
        std::vector<cPt2dr> mVPtsMod;     ///<  points sampling the model
        std::vector<double> mValueMod;    ///<  values of points in mVPtsMod
    };

    cTriangleDeformation::cTriangleDeformation(bool Show, cPt2di ImageSize) : mShow(Show),
                                                                              mSzIm(ImageSize),
                                                                              mIm(mSzIm),
                                                                              mDIm(mIm.DIm()),
                                                                              mSys(nullptr),
                                                                              mEqHomTri(nullptr)
    {
    }

    cTriangleDeformation::~cTriangleDeformation() 
    {
    delete mSys;
    delete mEqHomTri;
    }

    void cTriangleDeformation::OneIterationFitModele(bool IsLast)
    {
        //----------- index of unkown, basic here because the unknown are the same for each equation
        std::vector<int> aVecInd{0, 1, 2};
        //----------- allocate vec of obs : 6 for image, 3 for model
        std::vector<double> aVObs(4, 0.0);

        //----------- extract current parameters
        cDenseVect<double> aVCur = mSys->CurGlobSol();
        cHomot2D<tREAL8> aCurHomM2I(cPt2dr(aVCur(3), aVCur(4)), aVCur(2)); // current homothety
        // double aCurScR = aVCur(0);                                         // current scale on radiometry
        double aCurTrR = aVCur(0);                                            // current translation on radiometry

        //----------- declaration of indicator of convergence
        double aSomDif = 0; // sum of difference between model and image
        double aSomMod = 0; // sum of value of model, to normalize the difference
        double aNbOut = 0;  //  number of points out of image

        // Parse all the point to add the observations on each point
        for (size_t aKPt = 0; aKPt < mVPtsMod.size(); aKPt++)
        {
            cPt2dr aPMod = mVPtsMod[aKPt];         // point of model
            cPt2dr aPIm = aCurHomM2I.Value(aPMod); // image of aPMod by current homotethy
            if (mDIm.InsideBL(aPIm))               // avoid error
            {
                // put observations in vectors
                //  observations on image and point-image
                // Set_FormalBilinIm2D_Obs(aVObs, 0, aPIm, mDIm);
                FormalBilinIm2D_SetObs(aVObs, 0, aPIm, mDIm);

                //  observation point model and value model
                aVObs[6] = aPMod.x();
                aVObs[7] = aPMod.y();
                aVObs[8] = mValueMod[aKPt];

                // Now add observation
                mSys->CalcAndAddObs(mEqHomTri, aVecInd, aVObs);

                // compute indicator
                double aDif = mDIm.GetVBL(aPIm) - (aCurTrR + mValueMod[aKPt]); // residual
                aSomMod += mValueMod[aKPt];
                aSomDif += std::abs(aDif);
            }
            else
                aNbOut++;
        }
        //   Update all parameter taking into account previous observation
        mSys->SolveUpdateReset();

        if (mShow)
            StdOut() << " Dif= " << aSomDif / aSomMod << " NbOut=" << aNbOut << std::endl;
        //  If we are at end, check that the model is equal (up to numerical accuracy)  to the target
        if (IsLast)
        {
            MMVII_INTERNAL_ASSERT_bench(aNbOut == 0, "Gauss compos");
            MMVII_INTERNAL_ASSERT_bench((aSomDif / aSomMod) < 1e-10, "Gauss compos");
        }
    }

    /********************************************/
    /*              ::MMVII                     */
    /********************************************/

    void ComputeTriangleDeformation()
    {
        cPt2di SizeOfImage = {500, 300};
        cTriangleDeformation aTriangleDeform(false, SizeOfImage);
        int aNbIters = 8;
        for (int aIter = 0; aIter < aNbIters; aIter++)
            aTriangleDeform.OneIterationFitModele(aIter == (aNbIters - 1));
    }

}; // namespace MMVII