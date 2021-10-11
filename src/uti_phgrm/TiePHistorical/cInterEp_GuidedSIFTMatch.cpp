/*Header-MicMac-eLiSe-25/06/2007

    MicMac : Multi Image Correspondances par Methodes Automatiques de Correlation
    eLiSe  : ELements of an Image Software Environnement

    www.micmac.ign.fr


    Copyright : Institut Geographique National
    Author : Marc Pierrot Deseilligny
    Contributors : Gregoire Maillet, Didier Boldo.

[1] M. Pierrot-Deseilligny, N. Paparoditis.
    "A multiresolution and optimization-based image matching approach:
    An application to surface reconstruction from SPOT5-HRS stereo imagery."
    In IAPRS vol XXXVI-1/W41 in ISPRS Workshop On Topographic Mapping From Space
    (With Special Emphasis on Small Satellites), Ankara, Turquie, 02-2006.

[2] M. Pierrot-Deseilligny, "MicMac, un lociel de mise en correspondance
    d'images, adapte au contexte geograhique" to appears in
    Bulletin d'information de l'Institut Geographique National, 2007.

Francais :

   MicMac est un logiciel de mise en correspondance d'image adapte
   au contexte de recherche en information geographique. Il s'appuie sur
   la bibliotheque de manipulation d'image eLiSe. Il est distibue sous la
   licences Cecill-B.  Voir en bas de fichier et  http://www.cecill.info.


English :

    MicMac is an open source software specialized in image matching
    for research in geographic information. MicMac is built on the
    eLiSe image library. MicMac is governed by the  "Cecill-B licence".
    See below and http://www.cecill.info.

Header-MicMac-eLiSe-25/06/2007*/

#include "TiePHistorical.h"



/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant à la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est régi par la licence CeCILL-B soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffusée par le CEA, le CNRS et l'INRIA
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant
donné sa spécificité de logiciel libre, qui peut le rendre complexe à
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement,
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
pris connaissance de la licence CeCILL-B, et que vous en avez accepté les
termes.
aooter-MicMac-eLiSe-25/06/2007*/

extern ElSimilitude SimilRobustInit(const ElPackHomologue & aPackFull,double aPropRan,int aNbTir);

void FilterKeyPt(std::vector<Siftator::SiftPoint> aVSIFTPt, std::vector<Siftator::SiftPoint>& aVSIFTPtNew, double dMinScale, double dMaxScale)
{
    int nSizeL = aVSIFTPt.size();

    for(int i=0; i<nSizeL; i++)
    {
        if(aVSIFTPt[i].scale < dMaxScale && aVSIFTPt[i].scale > dMinScale)
            aVSIFTPtNew.push_back(aVSIFTPt[i]);
    }
}

void GetMinMaxScale(std::vector<Siftator::SiftPoint> aVSIFTPt, double & dMinScale, double & dMaxScale)
{
    int nSizeL = aVSIFTPt.size();

    if(nSizeL == 0)
        return;

    dMinScale = aVSIFTPt[0].scale;
    dMaxScale = aVSIFTPt[0].scale;

    for(int i=0; i<nSizeL; i++)
    {
        if(aVSIFTPt[i].scale > dMaxScale)
            dMaxScale = aVSIFTPt[i].scale;
        if(aVSIFTPt[i].scale < dMinScale)
            dMinScale = aVSIFTPt[i].scale;
    }
}

//use co-registered orientation or DSM to predict key points in another image
void PredictKeyPt(std::string aImg, std::vector<Pt2dr>& aVPredL, std::vector<Siftator::SiftPoint> aVSiftL, std::string aDSMFileL, std::string aDSMDirL, std::string aNameOriL, std::string aNameOriR, cTransform3DHelmert aTrans3DH, bool bPrint)
{
    //bool bDSM = false;

    cGet3Dcoor a3DL(aNameOriL);
    cDSMInfo aDSMInfoL = a3DL.SetDSMInfo(aDSMFileL, aDSMDirL);
    /*
    TIm2D<float,double> aTImProfPxL(a3DL.GetDSMSz(aDSMFileL, aDSMDirL));
    if(aDSMDirL.length() > 0)
    {
        bDSM = true;
        aTImProfPxL = a3DL.SetDSMInfo(aDSMFileL, aDSMDirL);
    }
    */
    cGet3Dcoor a3DR(aNameOriR);

    int nSizeL = aVSiftL.size();

    double dGSD1 = a3DL.GetGSD();
    printf("GSD of %s: %lf\n", aImg.c_str(), dGSD1);

    //printf("---------------------\n");
    for(int i=0; i<nSizeL; i++)
    {
        Pt2dr aPL = Pt2dr(aVSiftL[i].x, aVSiftL[i].y);
        Pt3dr aPTer1;
        /*if(bDSM == true)
        {*/
            bool bPreciseL;
            aPTer1 = a3DL.Get3Dcoor(aPL, aDSMInfoL, bPreciseL, bPrint);//, dGSD1);
        /*}
        else
        {
            aPTer1 = a3DL.GetRough3Dcoor(aPL);
        }*/

        aPTer1 = aTrans3DH.Transform3Dcoor(aPTer1);
        Pt2dr aPLPred = a3DR.Get2Dcoor(aPTer1);

        aVPredL.push_back(aPLPred);
    }
}

void SetAngleToValidRange(double& dAngle, double d2PI)
{
    while(dAngle > d2PI)
        dAngle = dAngle - d2PI;
    while(dAngle < 0)
        dAngle = dAngle + d2PI;
}
//for the key points in one image (master or secondary image), find their nearest neighbor in another image and record it
void MatchOneWay(std::vector<int>& matchIDL, std::vector<Siftator::SiftPoint> aVSiftL, std::vector<Siftator::SiftPoint> aVSiftR, std::vector<Pt2dr> aVPredL, Pt2di ImgSzR, double dScale, double dAngle, double threshScale, double threshAngle, bool bCheckScale, bool bCheckAngle, double dSearchSpace, bool bPredict, bool bRatioT)
{
    int nSIFT_DESCRIPTOR_SIZE = 128;
    const double d2PI = 3.1415926*2;

    int nSizeL = aVSiftL.size();
    int nSizeR = aVSiftR.size();

    int nStartL = 0;
    int nStartR = 0;
    int nEndL = nSizeL;
    int nEndR = nSizeR;

    std::time_t t1 = std::time(nullptr);
    std::cout << std::put_time(std::localtime(&t1), "%Y-%m-%d %H:%M:%S") << std::endl;

    long nSkiped = 0;
    float alpha = 2;

    int i, j, k;
    int nProgress = nSizeL/10;
    for(i=nStartL; i<nEndL; i++)
    {
        if(i%nProgress == 0)
        {
            printf("%.2lf%%\n", i*100.0/nSizeL);
        }

        double dBoxLeft  = 0;
        double dBoxRight = 0;
        double dBoxUpper = 0;
        double dBoxLower = 0;

        double dEuDisMin = DBL_MAX;
        double dEuDisSndMin = DBL_MAX;
        int nMatch = -1;

        double x, y;
        if(bPredict == true)
        {
            x = aVPredL[i].x;
            y = aVPredL[i].y;

            //if predicted point is out of the border of the other image, skip searching
            if(x<0 || x> ImgSzR.x || y<0 || y>ImgSzR.y)
            {
                matchIDL.push_back(-1);
                continue;
            }

            dBoxLeft  = x - dSearchSpace;
            dBoxRight = x + dSearchSpace;
            dBoxUpper = y - dSearchSpace;
            dBoxLower = y + dSearchSpace;
        }

        for(j=nStartR; j<nEndR; j++)
        {
            if(bPredict == true)
            {
                if(aVSiftR[j].x<=dBoxLeft || aVSiftR[j].x>= dBoxRight || aVSiftR[j].y<= dBoxUpper || aVSiftR[j].y>= dBoxLower)
                {
                    nSkiped++;
                    continue;
                }
            }
            if(bCheckScale == true)
            {
                /*
                double dScaleDif = fabs(aVSiftR[j].scale/aVSiftL[i].scale - dScale);
                if(dScaleDif > threshScale)
                    continue;
                    */
                double dScaleRatio = aVSiftR[j].scale/aVSiftL[i].scale;
                if((dScaleRatio < dScale*(1-threshScale)) || (dScaleRatio > dScale*(1+threshScale)))
                    continue;
                //printf("%.2lf ", dScaleRatio);
            }
            if(bCheckAngle == true)
            {
                /*
                double dAngleDif = fabs(aVSiftR[j].angle-aVSiftL[i].angle - dAngle);
                if((dAngleDif > threshAngle) && (dAngleDif < d2PI-threshAngle))
                    continue;
                    */
                double dAngleDif = aVSiftR[j].angle-aVSiftL[i].angle;
                SetAngleToValidRange(dAngleDif, d2PI);
                if((dAngleDif < dAngle-threshAngle) || (dAngleDif > dAngle+threshAngle))
                    continue;
                //printf("[%.2lf] ", dAngleDif);
            }

            double dDis = 0;
            for(k=0; k<nSIFT_DESCRIPTOR_SIZE; k++)
            {
                double dDif = aVSiftL[i].descriptor[k] - aVSiftR[j].descriptor[k];
                dDis += pow(dDif, alpha);
            }
            dDis = pow(dDis, 1.0/alpha);

            //save master and secondary nearest neigbor
            if(dDis < dEuDisMin)
            {
                dEuDisMin = dDis;
                nMatch = j;
                if(dEuDisMin > dEuDisSndMin)
                {
                    dEuDisSndMin = dEuDisMin;
                }
            }
            else if(dDis < dEuDisSndMin)
            {
                dEuDisSndMin = dDis;
            }
        }

        if(bRatioT == true && dEuDisMin/dEuDisSndMin > 0.8)
            nMatch = -1;
        matchIDL.push_back(nMatch);
    }
    std::time_t t2 = std::time(nullptr);
    std::cout << std::put_time(std::localtime(&t2), "%Y-%m-%d %H:%M:%S") << std::endl;
}

void MutualNearestNeighbor(bool bMutualNN, std::vector<int> matchIDL, std::vector<int> matchIDR, std::vector<Pt2di> & match)
{
    int nStartL = 0;
    int nStartR = 0;
    int nEndL = matchIDL.size();
    int nEndR = matchIDR.size();

    int i, j;
    if (bMutualNN == true){
        printf("Mutual nearest neighbor applied.\n");
        for(i=nStartL; i<nEndL; i++)
        {
            j = matchIDL[i-nStartL];

            if(j-nStartR < 0 || j-nStartR >= nEndR)
                 continue;
            if(matchIDR[j-nStartR] == i)
            {
                    Pt2di mPair = Pt2di(i, j);
                    match.push_back(mPair);
            }
        }
    }
    else
    {
        printf("Mutual nearest neighbor NOT applied.\n");
        for(i=nStartL; i<nEndL; i++)
        {
            j = matchIDL[i-nStartL];
            if(j-nStartR < 0 || j-nStartR >= nEndR)
                 continue;
            Pt2di mPair = Pt2di(i, j);
            match.push_back(mPair);

            //if the current pair is not mutual, save the other pair
            int nMatch4j = matchIDR[j-nStartR];
            if(nMatch4j != i && nMatch4j >= nStartL && nMatch4j-nStartL<nEndL)
            {
                Pt2di mPair = Pt2di(i, j);
                match.push_back(mPair);
            }
        }
    }
}

//transform the descriptor to rootSIFT if neccessary
void AmendSIFTKeyPt(std::vector<Siftator::SiftPoint>& aVSiftL, bool aRootSift)
{
    int nSizeL = aVSiftL.size();

    int nSIFT_DESCRIPTOR_SIZE = 128;
    for(int i=0; i<nSizeL; i++)
    {
        if(aRootSift == true)
        {
            double dSum = 0;
            for(int j=0; j<nSIFT_DESCRIPTOR_SIZE; j++)
                dSum += aVSiftL[i].descriptor[j];
            for(int j=0; j<nSIFT_DESCRIPTOR_SIZE; j++)
            {
                aVSiftL[i].descriptor[j] = sqrt(aVSiftL[i].descriptor[j]/dSum);
            }
        }
    }
}

Pt2dr GetScaleRotate(std::string aImg1, std::string aImg2, std::string aDSMFileL, std::string aDSMFileR, std::string aDSMDirL, std::string aDSMDirR, std::string aOri1, std::string aOri2, cInterfChantierNameManipulateur * aICNM, cTransform3DHelmert aTrans3DHL, cTransform3DHelmert aTrans3DHR, bool bPrint)
{
    if (ELISE_fp::exist_file(aImg1) == false || ELISE_fp::exist_file(aImg2) == false)
    {
        cout<<aImg1<<" or "<<aImg2<<" didn't exist, hence set ScaleRef=1 and RotateRef=0 instead"<<endl;
        return Pt2dr(1,0);
    }

    Tiff_Im aRGBIm1(aImg1.c_str());
    Pt2di ImgSzL = aRGBIm1.sz();
    Tiff_Im aRGBIm2(aImg2.c_str());
    Pt2di ImgSzR = aRGBIm2.sz();

    std::string aNameOriL = aICNM->StdNameCamGenOfNames(aOri1, aImg1);
    std::string aNameOriR = aICNM->StdNameCamGenOfNames(aOri2, aImg2);
    cGet3Dcoor a3DL(aNameOriL);
    cDSMInfo aDSMInfoL = a3DL.SetDSMInfo(aDSMFileL, aDSMDirL);
    /*
    TIm2D<float,double> aTImProfPxL(a3DL.GetDSMSz(aDSMFileL, aDSMDirL));
    bool bDSML = false;
    if(aDSMDirL.length() > 0)
    {
        bDSML = true;
        aTImProfPxL = a3DL.SetDSMInfo(aDSMFileL, aDSMDirL);
    }
    */
    cGet3Dcoor a3DR(aNameOriR);
    cDSMInfo aDSMInfoR = a3DR.SetDSMInfo(aDSMFileR, aDSMDirR);
    /*
    TIm2D<float,double> aTImProfPxR(a3DR.GetDSMSz(aDSMFileR, aDSMDirR));
    bool bDSMR = false;
    if(aDSMDirR.length() > 0)
    {
        bDSMR = true;
        aTImProfPxR = a3DR.SetDSMInfo(aDSMFileR, aDSMDirR);
    }
    */

    Pt2dr aPCornerL[4];
    Pt2dr origin = Pt2dr(0, 0);
    aPCornerL[0] = origin;
    aPCornerL[1] = Pt2dr(origin.x+ImgSzL.x, origin.y);
    aPCornerL[2] = Pt2dr(origin.x+ImgSzL.x, origin.y+ImgSzL.y);
    aPCornerL[3] = Pt2dr(origin.x, origin.y+ImgSzL.y);

    //for the 4 corners in master images, get corresponding points in secondary image, then get overlapping area of the 2 images in the frame of secondary image
    Pt2dr aPLPredinR[4];
    for(int i=0; i<4; i++)
    {
        Pt2dr aPL = aPCornerL[i];
        Pt3dr aPTer1;
        /*if(bDSML == true)
        {*/
            bool bPreciseL;
            aPTer1 = a3DL.Get3Dcoor(aPL, aDSMInfoL, bPreciseL, bPrint);//, a3DL.GetGSD());
        /*}
        else
        {
            aPTer1 = a3DL.GetRough3Dcoor(aPL);
        }*/

        aPTer1 = aTrans3DHL.Transform3Dcoor(aPTer1);
        Pt2dr ptPred = a3DR.Get2Dcoor(aPTer1);
        aPLPredinR[i] = ptPred;

        CheckRange(0, ImgSzR.x, aPLPredinR[i].x);
        CheckRange(0, ImgSzR.y, aPLPredinR[i].y);

        if(bPrint)
        {
            printf("%dth: CornerL: [%.2lf\t%.2lf], ImEtProf2Terrain: [%.2lf\t%.2lf\t%.2lf], CornerR: [%.2lf\t%.2lf], CornerRNew: [%.2lf\t%.2lf]\n", i, aPCornerL[i].x, aPCornerL[i].y, aPTer1.x, aPTer1.y, aPTer1.z, ptPred.x, ptPred.y, aPLPredinR[i].x, aPLPredinR[i].y);
        }
    }

    double dMaxX = max(max(aPLPredinR[0].x, aPLPredinR[1].x), max(aPLPredinR[2].x, aPLPredinR[3].x));
    double dMinX = min(min(aPLPredinR[0].x, aPLPredinR[1].x), min(aPLPredinR[2].x, aPLPredinR[3].x));
    double dMaxY = max(max(aPLPredinR[0].y, aPLPredinR[1].y), max(aPLPredinR[2].y, aPLPredinR[3].y));
    double dMinY = min(min(aPLPredinR[0].y, aPLPredinR[1].y), min(aPLPredinR[2].y, aPLPredinR[3].y));
    //printf("dMaxX: %lf\tdMaxX: %lf\tdMaxX: %lf\tdMaxX: %lf\n", dMaxX, dMinX, dMaxY, dMinY);

    Pt2dr aPtL[2];
    Pt2dr aPtR[2];
    aPtR[0] = Pt2dr(dMinX, dMinY);
    aPtR[1] = Pt2dr(dMaxX, dMaxY);

    //for the 2 ends of the diagonal of the overlapping area, get the corresponding points in the frame of master image
    ElPackHomologue aPackSeed;
    for(int i=0; i<2; i++)
    {
        Pt2dr aPR = aPtR[i];
        Pt3dr aPTer1;
        /*if(bDSMR == true)
        {*/
            bool bPreciseR;
            aPTer1 = a3DR.Get3Dcoor(aPR, aDSMInfoR, bPreciseR, bPrint);//, a3DL.GetGSD());
        /*}
        else
        {
            aPTer1 = a3DR.GetRough3Dcoor(aPR);
        }*/

        aPTer1 = aTrans3DHR.Transform3Dcoor(aPTer1);
        Pt2dr ptPred = a3DL.Get2Dcoor(aPTer1);
        aPtL[i] = ptPred;

        aPackSeed.Cple_Add(ElCplePtsHomologues(aPtL[i],aPtR[i]));
        printf("%dth end of the diagonal of the overlapping area of the image pair: %.2lf %.2lf %.2lf %.2lf\n",i, aPtL[i].x,aPtL[i].y,aPtR[i].x,aPtR[i].y);
    }

    //based on the 2 correspoinding points of the overlapping area, calculate the scale and rotate of the master and secondary images
    double aPropRan = 0.8;
    ElSimilitude aSimCur = SimilRobustInit(aPackSeed,aPropRan,1);

    Pt2dr sc = aSimCur.sc();
    Pt2dr aRes;
    //scale
    aRes.x = pow(sc.x*sc.x + sc.y*sc.y, 0.5);
    int nSign = (asin(sc.y/aRes.x) > 0) ? 1 : -1;
    //angle
    aRes.y = nSign*acos(sc.x/aRes.x);

    return aRes;
}

void GuidedSIFTMatch(std::string aDir,std::string aImg1, std::string aImg2, std::string outSH, std::string aDSMFileL, std::string aDSMFileR, std::string aDSMDirL, std::string aDSMDirR, std::string aOri1, std::string aOri2, cInterfChantierNameManipulateur * aICNM, bool bRootSift, double dSearchSpace, bool bPredict, bool bRatioT, bool bMutualNN, cTransform3DHelmert aTrans3DHL, cTransform3DHelmert aTrans3DHR, bool bCheckScale, bool bCheckAngle, bool bPrint, bool bCheckFile, double threshScale, double threshAngle)//, double dScale=1, double dAngle=0)
{
    if(bRatioT == true)
        printf("Ratio test applied.\n");
    else
        printf("Ratio test NOT applied.\n");

    if (ELISE_fp::exist_file(aImg1) == false || ELISE_fp::exist_file(aImg2) == false)
    {
        cout<<aImg1<<" or "<<aImg2<<" didn't exist, hence skipped"<<endl;
        return;
    }

    std::string aSHDir = aDir + "/Homol" + outSH + "/";
    ELISE_fp::MkDir(aSHDir);
    std::string aNewDir = aSHDir + "Pastis" + aImg1;
    ELISE_fp::MkDir(aNewDir);
    std::string aNameFile1 = aNewDir + "/"+aImg2+".txt";

    aNewDir = aSHDir + "Pastis" + aImg2;
    ELISE_fp::MkDir(aNewDir);
    std::string aNameFile2 = aNewDir + "/"+aImg1+".txt";
    if (bCheckFile == true && ELISE_fp::exist_file(aNameFile1) == true && ELISE_fp::exist_file(aNameFile2) == true)
    {
        cout<<aNameFile1<<" already exist, hence skipped"<<endl;
        return;
    }

    Tiff_Im aRGBIm1(aImg1.c_str());
    Pt2di ImgSzL = aRGBIm1.sz();
    Tiff_Im aRGBIm2(aImg2.c_str());
    Pt2di ImgSzR = aRGBIm2.sz();

    std::string aImg1Key = aImg1.substr(0, aImg1.rfind(".")) + ".key";
    std::string aImg2Key = aImg2.substr(0, aImg2.rfind(".")) + ".key";

    //*********** 0. calculate ScaleTh and AngleTh for CheckScale and CheckAngle
    double d2PI = 3.1415926*2;
    Pt2dr ScaleRotateL = Pt2dr(1, 0);
    Pt2dr ScaleRotateR = Pt2dr(1, 0);
    if(bCheckScale == true || bCheckAngle == true)
    {
        ScaleRotateL = GetScaleRotate(aImg1, aImg2, aDSMFileL, aDSMFileR, aDSMDirL, aDSMDirR, aOri1, aOri2, aICNM, aTrans3DHL, aTrans3DHR, bPrint);
        SetAngleToValidRange(ScaleRotateL.y, d2PI);
        ScaleRotateR.x = 1.0/ScaleRotateL.x;
        ScaleRotateR.y = d2PI-ScaleRotateL.y;
        SetAngleToValidRange(ScaleRotateR.y, d2PI);
        printf("From master image to secondary image (for CheckScale and CheckAngle): \n");
        printf("    Reference of scale: [%.2lf], ScaleTh: [%.2lf]; Range of valid scale ratio: [%.2lf, %.2lf]\n", ScaleRotateL.x, threshScale, ScaleRotateL.x*(1-threshScale), ScaleRotateL.x*(1+threshScale));
        printf("    Reference of angle: [%.2lf], AngleTh: [%.2lf]; Range of valid angle difference: [%.2lf, %.2lf]\n", ScaleRotateL.y, threshAngle, ScaleRotateL.y-threshAngle, ScaleRotateL.y+threshAngle);

        printf("From secondary image to master image (for CheckScale and CheckAngle): \n");
        printf("    Reference of scale: [%.2lf], ScaleTh: [%.2lf]; Range of valid scale ratio: [%.2lf, %.2lf]\n", ScaleRotateR.x, threshScale, ScaleRotateR.x*(1-threshScale), ScaleRotateR.x*(1+threshScale));
        printf("    Reference of angle: [%.2lf], AngleTh: [%.2lf]; Range of valid angle difference: [%.2lf, %.2lf]\n", ScaleRotateR.y, threshAngle, ScaleRotateR.y-threshAngle, ScaleRotateR.y+threshAngle);
    }
    if(bCheckScale == false && bCheckAngle == false)
        printf("won't check scale and rotate\n");

    //*********** 1. read SIFT key-pts
    std::vector<Siftator::SiftPoint> aVSiftOriL;
    std::vector<Siftator::SiftPoint> aVSiftOriR;
    if(read_siftPoint_list(aImg1Key,aVSiftOriL) == false || read_siftPoint_list(aImg2Key,aVSiftOriR) == false)
    {
        cout<<"Read SIFT of "<<aImg1Key<<" or "<<aImg2Key<<" went wrong."<<endl;
        return;
    }

    double dMinScaleL, dMaxScaleL;
    GetMinMaxScale(aVSiftOriL, dMinScaleL, dMaxScaleL);
    double dMinScaleR, dMaxScaleR;
    GetMinMaxScale(aVSiftOriR, dMinScaleR, dMaxScaleR);


    std::vector<Siftator::SiftPoint> aVSiftL;
    FilterKeyPt(aVSiftOriL, aVSiftL, dMinScaleR*ScaleRotateR.x*(1-threshScale), dMaxScaleR*ScaleRotateR.x*(1+threshScale));
    std::vector<Siftator::SiftPoint> aVSiftR;
    FilterKeyPt(aVSiftOriR, aVSiftR, dMinScaleL*ScaleRotateL.x*(1-threshScale), dMaxScaleL*ScaleRotateL.x*(1+threshScale));

    printf("Original key point number of master image: %d. (With scales between [%.2lf, %.2lf].)\nOriginal key point number of secondary image: %d. (With scales between [%.2lf, %.2lf].)\nFiltered key point number of master image: %d. (Only key points with scale between [%.2lf, %.2lf] are kept.)\niltered key point number of secondary image: %d. (Only key points with scale between [%.2lf, %.2lf] are kept.)\n", int(aVSiftOriL.size()), dMinScaleL, dMaxScaleL, int(aVSiftOriR.size()), dMinScaleR, dMaxScaleR, int(aVSiftL.size()), dMinScaleR*ScaleRotateR.x*(1-threshScale), dMaxScaleR*ScaleRotateR.x*(1+threshScale), int(aVSiftR.size()), dMinScaleL*ScaleRotateL.x*(1-threshScale), dMaxScaleL*ScaleRotateL.x*(1+threshScale));
/*
    printf("Original key point number of master image: %d. (With scales between [%.2lf, %.2lf].)\n", int(aVSiftOriL.size()), dMinScaleL, dMaxScaleL);
    printf("Original key point number of secondary image: %d. (With scales between [%.2lf, %.2lf].)\n", int(aVSiftOriR.size()), dMinScaleR, dMaxScaleR);
    printf("Filtered key point number of master image: %d. (Only key points with scale between [%.2lf, %.2lf] are kept.)\n", int(aVSiftL.size()), dMinScaleR*ScaleRotateR.x*(1-threshScale), dMaxScaleR*ScaleRotateR.x*(1+threshScale));
    printf("Filtered key point number of secondary image: %d. (Only key points with scale between [%.2lf, %.2lf] are kept.)\n", int(aVSiftR.size()), dMinScaleL*ScaleRotateL.x*(1-threshScale), dMaxScaleL*ScaleRotateL.x*(1+threshScale));
*/
    //transform the descriptor to rootSIFT if neccessary
    if(bRootSift == true){
        printf("Use RootSIFT as descriptor.\n");
        AmendSIFTKeyPt(aVSiftL, bRootSift);
        AmendSIFTKeyPt(aVSiftR, bRootSift);
    }

    //*********** 2. Predict SIFT key-pts
    std::vector<Pt2dr> aVPredL;
    std::vector<Pt2dr> aVPredR;
    std::string aNameOriL = aICNM->StdNameCamGenOfNames(aOri1, aImg1);
    std::string aNameOriR = aICNM->StdNameCamGenOfNames(aOri2, aImg2);
    PredictKeyPt(aImg1, aVPredL, aVSiftL, aDSMFileL, aDSMDirL, aNameOriL, aNameOriR, aTrans3DHL, bPrint);
    PredictKeyPt(aImg2, aVPredR, aVSiftR, aDSMFileR, aDSMDirR, aNameOriR, aNameOriL, aTrans3DHR, bPrint);

    //*********** 3. match SIFT key-pts
    threshAngle = 3.14*threshAngle/180;
    /*
    double threshScale, threshAngle;
    threshScale = 0.2;//0.05;
    threshAngle = 3.14*30/180;//3.14*10/180;
    */
    /*
    if(bCheckScale == true)
        printf("Check scale. dScale: %lf; threshScale: %lf\n", dScale, threshScale);
    else
        printf("won't check scale\n");
    if(bCheckAngle == true)
        printf("Check angle. dAngle: %lf; threshAngle: %lf\n", dAngle, threshAngle);
    else
        printf("won't check angle\n");
    */

    std::vector<int> matchIDL;
    std::vector<int> matchIDR;
    printf("*****processing Left*****\n");
    if(bPredict)
        printf("SearchSpace = %.2lf when searching %s\n", dSearchSpace*ScaleRotateL.x, aImg2.c_str());
    MatchOneWay(matchIDL, aVSiftL, aVSiftR, aVPredL, ImgSzR, ScaleRotateL.x, ScaleRotateL.y, threshScale, threshAngle, bCheckScale, bCheckAngle, dSearchSpace*ScaleRotateL.x, bPredict, bRatioT);
    printf("*****processing Right*****\n");
    if(bPredict)
        printf("SearchSpace = %.2lf when searching %s\n", dSearchSpace, aImg1.c_str());
    MatchOneWay(matchIDR, aVSiftR, aVSiftL, aVPredR, ImgSzL, ScaleRotateR.x, ScaleRotateR.y, threshScale, threshAngle, bCheckScale, bCheckAngle, dSearchSpace, bPredict, bRatioT);

    //cout<<matchIDL.size()<<",,,,"<<matchIDR.size()<<endl;

    std::vector<Pt2di> match;
    MutualNearestNeighbor(bMutualNN, matchIDL, matchIDR, match);

    //*********** 4. Save tie pt
    FILE * fpTiePt1 = fopen(aNameFile1.c_str(), "w");
    FILE * fpTiePt2 = fopen(aNameFile2.c_str(), "w");

    int nTiePtNum = match.size();
    for(int i = 0; i < nTiePtNum; i++)
    {
        int idxL = match[i].x;
        int idxR = match[i].y;
        Pt2dr aP1, aP2;
        aP1 = Pt2dr(aVSiftL[idxL].x, aVSiftL[idxL].y);
        aP2 = Pt2dr(aVSiftR[idxR].x, aVSiftR[idxR].y);
        fprintf(fpTiePt1, "%lf %lf %lf %lf\n", aP1.x, aP1.y, aP2.x, aP2.y);
        fprintf(fpTiePt2, "%lf %lf %lf %lf\n", aP2.x, aP2.y, aP1.x, aP1.y);

        if(bPrint){
            printf("%dth match, master and secondary scale and angle: %lf %lf %lf %lf\n", i, aVSiftL[idxL].scale, aVSiftL[idxL].angle, aVSiftR[idxR].scale, aVSiftR[idxR].angle);
            printf("scaleDif, angleDif: %lf %lf\n", aVSiftR[idxR].scale/aVSiftL[idxL].scale,aVSiftR[idxR].angle-aVSiftL[idxL].angle);
        }
    }
    fclose(fpTiePt1);
    fclose(fpTiePt2);

    cout<<"Extracted tie point number: "<<match.size()<<endl;

    std::string aCom = "mm3d SEL" + BLANK + aDir + BLANK + aImg1 + BLANK + aImg2 + BLANK + "KH=NT SzW=[600,600] SH="+outSH;
    std::string aComInv = "mm3d SEL" + BLANK + aDir + BLANK + aImg2 + BLANK + aImg1 + BLANK + "KH=NT SzW=[600,600] SH="+outSH;
    printf("%s\n%s\n", aCom.c_str(), aComInv.c_str());
}

void ExtractSIFT(std::string aFullName, std::string aDir)
{
    cInterfChantierNameManipulateur::BasicAlloc(DirOfFile(aFullName));
    cout<<aFullName<<endl;

    //Tiff_Im::StdConvGen(aFullName,1,true,true);
    Tiff_Im::StdConvGen(aFullName,1,false,true);

    std::string aGrayImgName = aFullName + "_Ch1.tif";

    //if RGB image
    if( ELISE_fp::exist_file(aDir + "/Tmp-MM-Dir/" + aGrayImgName) == true)
    {
        std::string aComm;
        aComm = "mv " + aDir + "/Tmp-MM-Dir/" + aGrayImgName + " " + aGrayImgName;
        cout<<aComm<<endl;
        System(aComm);

        aComm = MMBinFile(MM3DStr) + "SIFT " + aGrayImgName;
        cout<<aComm<<endl;
        System(aComm);

        aComm = "mv " + StdPrefix(aGrayImgName)+".key" + " "+StdPrefix(aFullName)+".key";
        cout<<aComm<<endl;
        System(aComm);

        aComm = "rm " + aGrayImgName;
        cout<<aComm<<endl;
        System(aComm);
    }
    //gray image
    else
    {
        std::string aCom = MMBinFile(MM3DStr) + "SIFT " + aFullName;
        cout<<aCom<<endl;
        System(aCom);
    }
}

int GuidedSIFTMatch_main(int argc,char ** argv)
{
   cCommonAppliTiepHistorical aCAS3D;

   std::string aImg1;
   std::string aImg2;
   std::string aOri1;
   std::string aOri2;

   std::string aDSMDirL = "";
   std::string aDSMDirR = "";
   std::string aDSMFileL;
   std::string aDSMFileR;

   aDSMFileL = "MMLastNuage.xml";
   aDSMFileR = "MMLastNuage.xml";

   std::string aPara3DHL = "";
   std::string aPara3DHR = "";

   bool bCheckFile = false;

   ElInitArgMain
    (
        argc,argv,
        LArgMain()   << EAMC(aImg1,"Master image name")
               << EAMC(aImg2,"Secondary image name")
               << EAMC(aOri1,"Orientation of master image")
               << EAMC(aOri2,"Orientation of secondary image"),
        LArgMain()
                    << aCAS3D.ArgBasic()
                    << aCAS3D.ArgGuidedSIFT()
               << EAM(aDSMDirL, "DSMDirL", true, "DSM of master image (for improving the reprojecting accuracy), Def=none")
               << EAM(aDSMDirR, "DSMDirR", true, "DSM of secondary image (for improving the reprojecting accuracy), Def=none")
               << EAM(aDSMFileL, "DSMFileL", true, "DSM File of master image, Def=MMLastNuage.xml")
               << EAM(aDSMFileR, "DSMFileR", true, "DSM File of secondary image, Def=MMLastNuage.xml")
               << EAM(aPara3DHL, "Para3DHL", false, "Input xml file that recorded the paremeter of the 3D Helmert transformation from orientation of master image to secondary image, Def=none")
               << EAM(aPara3DHR, "Para3DHR", false, "Input xml file that recorded the paremeter of the 3D Helmert transformation from orientation of secondary image to master image, Def=none")
               << EAM(bCheckFile, "CheckFile", true, "Check if the result files of inter-epoch correspondences exist (if so, skip to avoid repetition), Def=false")

    );
    StdCorrecNameOrient(aOri1,"./",true);
    StdCorrecNameOrient(aOri2,"./",true);
   //if SIFT key point is not extracted before (aCAS3D.mSkipSIFT = true), extract SIFT first
   if(aCAS3D.mSkipSIFT == false)
   {
       ExtractSIFT(aImg1, aCAS3D.mDir);
       ExtractSIFT(aImg2, aCAS3D.mDir);

   }
   cTransform3DHelmert aTrans3DHL(aPara3DHL);
   cTransform3DHelmert aTrans3DHR(aPara3DHR);

   GuidedSIFTMatch( aCAS3D.mDir, aImg1,  aImg2,  aCAS3D.mGuidedSIFTOutSH, aDSMFileL, aDSMFileR, aDSMDirL, aDSMDirR,  aOri1, aOri2, aCAS3D.mICNM, aCAS3D.mRootSift, aCAS3D.mSearchSpace, aCAS3D.mPredict, aCAS3D.mRatioT, aCAS3D.mMutualNN, aTrans3DHL, aTrans3DHR, aCAS3D.mCheckScale, aCAS3D.mCheckAngle, aCAS3D.mPrint, bCheckFile, aCAS3D.mThreshScale, aCAS3D.mThreshAngle);//, aCAS3D.mScale, aCAS3D.mAngle);

   return EXIT_SUCCESS;
}
