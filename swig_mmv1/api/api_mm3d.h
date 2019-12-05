#include "StdAfx.h"

/**
@file
@brief New methods for python API and existing classes
**/

class CamStenope;

//-------------------- Nouvelles methodes ---------------------

//! Create CamStenope form a XML file
CamStenope *  CamOrientFromFile(std::string filename);

//-------------------- classes MM a exporter ------------------
#ifdef SWIG

#include "general/CMake_defines.h"
#include "general/sys_dep.h"

typedef  Pt2d<INT> Pt2di;
typedef  Pt2d<double> Pt2dr;
typedef  Pt2d<long double> Pt2dlr;
typedef  Pt2d<float> Pt2df;
typedef  Pt2d<U_INT2> Pt2dUi2;

template <class Type> class Pt2d : public  ElStdTypeScal<Type>
{
   public :

     typedef typename TCompl<Type>::TypeCompl  tCompl;
     typedef Pt2d<REAL>  TypeReel;
     typedef Type        TypeScal;
     typedef Pt2d<Type>  TypeEff;
     static Pt2d  El0 () {return Pt2d(0,0);}

     typedef Pt2d<typename ElStdTypeScal<Type>::TypeVarProvReel> TypeProvPtScalR;

     typedef Type (& t2)[2] ;
     Type   x;
     Type   y;


  // Constructeur

     Pt2d<Type>()  : x (0), y (0) {}
     Pt2d<Type>(Type X,Type Y) : x (X), y (Y) {}

     Pt2d<Type>(const Pt2d<Type>& p) : x (p.x), y (p.y) {}
     explicit Pt2d<Type>(const Pt2d<tCompl>& p) :
              x( TCompl<Type>::FromC( p.x)),
              y( TCompl<Type>::FromC( p.y))
     {
     }

     static  Pt2d<Type> IP2ToThisT(const Pt2d<int> & aP){return Pt2d<Type>(Type(aP.x),Type(aP.y));}
     static  Pt2d<Type> RP2ToThisT(const Pt2d<double> & aP){return Pt2d<Type>(Type(aP.x),Type(aP.y));}
     static  Pt2d<Type> FP2ToThisT(const Pt2d<float> & aP){return Pt2d<Type>(Type(aP.x),Type(aP.y));}

/*

     Pt2d<Type>(const Pt2d<INT>& p) : x (p.x), y (p.y) {};
     Pt2d<Type>(const Pt2d<REAL>& p): x (Pt2d<Type>::RtoT(p.x)), y (Pt2d<Type>::RtoT(p.y)) {};
*/


     static  Pt2d<Type>  FromPolar(REAL rho,REAL teta)
     {
        return   Pt2d<Type>(ElStdTypeScal<Type>::RtoT(cos(teta)*rho),ElStdTypeScal<Type>::RtoT(sin(teta)*rho));
     }

     static Pt2d<double> polar(const Pt2d<double> & p,REAL AngDef);

 // Operateurs

         // unaires,  Pt => Pt

     TypeProvPtScalR  ToPtProvR() const
     {
           return TypeProvPtScalR (this->T2R(x),this->T2R(y));
     }

     Pt2d<Type> operator - () const { return Pt2d<Type>(-x,-y); }
     Pt2d<Type> yx() const { return Pt2d(y,x);}
     Pt2d<Type> conj() const { return Pt2d(x,-y);}
     Pt2d<typename ElStdTypeScal<Type>::TypeScalReel> inv() const
     {
         typename ElStdTypeScal<Type>::TypeVarProvReel  n= this->T2R(x)*x+y*y;
         return Pt2d<typename ElStdTypeScal<Type>::TypeScalReel>(x/n,-y/n);
     };
      Pt2d<Type> Square() const;
      Type XtY() const {return x * y;}


         // binaires,  PtxPt => Pt

     Pt2d<Type> operator + (const Pt2d<Type> & p2) const
                {return Pt2d<Type>(x+p2.x,y+p2.y);}
     Pt2d<Type> operator * (const Pt2d<Type> & p2) const
                {return Pt2d<Type>(x*p2.x-y*p2.y,x*p2.y+y*p2.x);}

     // TCompl
     Pt2d<Type> operator / (const Pt2d<Type> & p2) const
     {
            TypeProvPtScalR aRes = this->ToPtProvR() * p2.inv().ToPtProvR();
             return Pt2d<Type> ((Type)aRes.x,(Type)aRes.y);
     }

     Pt2d<Type> operator - (const Pt2d<Type> & p2) const
                {return Pt2d<Type>(x-p2.x,y-p2.y);}
     Pt2d<Type> mcbyc(const Pt2d<Type> & p2) const
                {return Pt2d(x*p2.x,y*p2.y);}
     Pt2d<Type> dcbyc(const Pt2d<Type> & p2) const
                {return Pt2d(x/p2.x,y/p2.y);}


     void SetSup(const Pt2d<Type> & p){ElSetMax(x,p.x); ElSetMax(y,p.y);}
     void SetInf(const Pt2d<Type> & p){ElSetMin(x,p.x); ElSetMin(y,p.y);}

    // RatioMin :  return Min ( x/(TypeScalReel)p.x, y/(TypeScalReel)p.y);
     typename Pt2d<Type>::TypeScalReel RatioMin(const Pt2d<Type> & p) const;

         // binnaire, affectation composee

     Pt2d<Type> & operator += (const Pt2d<Type> & p2)
                { x += p2.x; y += p2.y; return * this;}
     Pt2d<Type> & operator -= (const Pt2d<Type> & p2)
                { x -= p2.x; y -= p2.y; return * this;}

     Pt2d<Type>  &  operator = (const Pt2d<Type> & p2)
     {
            x = p2.x;
            y = p2.y;
            return * this;
     }

         // binaire,  PtxPt => bool
     typename ElStdTypeScal<Type>::TypeBool  operator == (const Pt2d<Type> & p2) const {return (x==p2.x) && (y==p2.y);}
     typename ElStdTypeScal<Type>::TypeBool  operator != (const Pt2d<Type> & p2) const {return (x!=p2.x) || (y!=p2.y);}
     // p1 < p2 , utile par ex ds les map<Pt2di,Machin>
     typename ElStdTypeScal<Type>::TypeBool  operator <  (const Pt2d<Type> & p2) const {return (x<p2.x) || ((x==p2.x)&&(y<p2.y));}

     typename ElStdTypeScal<Type>::TypeBool   xety_inf_ou_egal (const Pt2d<Type> & p2) const
            {return (x<=p2.x) && (y<=p2.y);}

         // binaires,  PtxScalaire => Pt

     Pt2d<Type> operator * (INT  lambda) const { return Pt2d<Type>(x*lambda,y*lambda);}


     Pt2d<typename ElStdTypeScal<Type>::TypeScalReel> operator * (REAL lambda) const { return Pt2d<typename ElStdTypeScal<Type>::TypeScalReel>(x*lambda,y*lambda);}

     Pt2d<Type> operator / (INT  lambda) const { return Pt2d<Type>(x/lambda,y/lambda);}
     Pt2d<typename ElStdTypeScal<Type>::TypeScalReel> operator / (REAL lambda) const { return Pt2d<typename ElStdTypeScal<Type>::TypeScalReel>(x/lambda,y/lambda);}


      // operator * est deja surcharge
      Pt2d<Type> mul (const Type & aL) const { return Pt2d<Type>(x*aL,y*aL);}
      Pt2d<Type> div (const Type & aL) const { return Pt2d<Type>(x/aL,y/aL);}


         // binaires,  PtxPt => scalaire

     Type  operator ^ (const Pt2d<Type> & p2) const{return x*p2.y-y*p2.x;}


          // lies a une distance
    //friend Type  dist4(const Pt2d<Type> & p){return ElAbs(p.x)+  ElAbs(p.y);}
    //friend Type  dist8(const Pt2d<Type> & p){return ElMax(ElAbs(p.x),ElAbs(p.y));}

     typename ElStdTypeScal<Type>::TypeBool in_box(const Pt2d<Type> & p0, const Pt2d<Type> & p1)
     {
         return (x>=p0.x)&&(y>=p0.y)&&(x< p1.x)&&(y<p1.y);
     }

     friend void pt_set_min_max<>(Pt2d<Type> & p0,Pt2d<Type> & p1);

            // tertiaire

     // in_sect_angulaire :  est que le pt est dans le secteur partant de p1
     // et defini par un parcourt trigo jusqu'a p2
     bool in_sect_angulaire(const Pt2d<Type> & p1,const Pt2d<Type> & p2) const;

     // Ceux-ci n'ont aucun interet a etre iniline


     void to_tab(Type (& t)[2] ) const;
     static Pt2d<Type> FromTab(const Type *);
     static Pt2d<Type> FromTab(const std::vector<Type> &);
     std::vector<Type> ToTab() const;
     Output sigma();
     Output VMax();
     Output VMin();
     Output WhichMax();
     Output WhichMin();

     Pt2d<Type> AbsP() const {return Pt2d<Type>(ElAbs(x),ElAbs(y));}

     private :
          void Verif_adr_xy();

};


//original
//#define Pt3di  Pt3d<INT>
//#define Pt3dr  Pt3d<REAL>
//#define Pt3df  Pt3d<float>
typedef  Pt3d<INT> Pt3di;
typedef  Pt3d<double> Pt3dr;
typedef  Pt3d<float> Pt3df;

template <class Type> class Pt3d : public  ElStdTypeScal<Type>
{
   public :
     typedef typename TCompl<Type>::TypeCompl  tCompl;
     Type   x;
     Type   y;
     Type   z;

     Pt3d();

     Pt3d<Type>(const Pt3d<Type>& ); // to please visual
     explicit Pt3d<Type>(const Pt3d<tCompl>& p);

     static  Pt3d<Type> P3ToThisT(const Pt3d<int> & aP){return Pt3d<Type>(Type(aP.x),Type(aP.y),Type(aP.z));}
     static  Pt3d<Type> P3ToThisT(const Pt3d<double> & aP){return Pt3d<Type>(Type(aP.x),Type(aP.y),Type(aP.z));}
     static  Pt3d<Type> P3ToThisT(const Pt3d<float> & aP){return Pt3d<Type>(Type(aP.x),Type(aP.y),Type(aP.z));}


     Pt3d<Type>(const Pt2d<Type>&,Type z); // to please visual



     Pt3d(Type X,Type Y,Type Z);
     Pt3d<Type> operator + (const Pt3d<Type> & p2) const;

     Pt3d<Type> operator * (Type) const;
     Pt3d<Type> operator / (Type) const;

     Pt3d<Type> operator - (const Pt3d & p2) const;
     Pt3d<Type> operator - () const;

     typename ElStdTypeScal<Type>::TypeBool  operator == (const Pt3d<Type> & p2) const {return (x==p2.x) && (y==p2.y) && (z==p2.z);}

     // multiplication coordinate by coordinate

     // friend Type  scal<Type> (const Pt3d<Type> & p1,const Pt3d<Type> & p2);
/*
     friend Type Det(const Pt3d<Type> & p1,const Pt3d<Type> & p2,const Pt3d<Type> & p3)
     {
         return scal(p1 ,p2^p3);
     }
*/

     Pt3d<Type>  operator ^ (const Pt3d<Type> & p2) const;
     Pt3d<Type>  &  operator = (const Pt3d<Type> & p2) ;


     void to_tab(Type (& t)[3] ) const;
     static Pt3d<Type> FromTab(const Type *);
     std::vector<Type> ToTab() const;
     static Pt3d<Type> FromTab(const std::vector<Type> &);

     Pt3d<Type> AbsP() const {return Pt3d<Type>(ElAbs(x),ElAbs(y),ElAbs(z));}
     /*
     friend Pt3d<Type> Sup (const Pt3d<Type> & p1,const Pt3d<Type> & p2)
           { return Pt3d<Type>(ElMax(p1.x,p2.x),ElMax(p1.y,p2.y),ElMax(p1.z,p2.z));}
     friend Pt3d<Type> Inf (const Pt3d<Type> & p1,const Pt3d<Type> & p2)
            { return Pt3d<Type>(ElMin(p1.x,p2.x),ElMin(p1.y,p2.y),ElMin(p1.z,p2.z));}
      */
     Output sigma();
     Output VMax();
     Output VMin();
     Output WhichMax();
     Output WhichMin();

     static Type instantiate();



     // ! Convention Phi = 0 "a l'equateur"
     static  Pt3d<Type>  TyFromSpherique(Type Rho,Type Teta,Type Phi)
     {
        return   Pt3d<Type>
         (
          ElStdTypeScal<Type>::RTtoT(cos(Phi)*cos(Teta)*Rho),
          ElStdTypeScal<Type>::RTtoT(cos(Phi)*sin(Teta)*Rho),
          ElStdTypeScal<Type>::RTtoT(sin(Phi)*Rho)
         );
     }

     private :
          void Verif_adr_xy();
};

class ElCamera : public cCapture3D
{
     public :

         const bool &   IsScanned() const;
         void  SetScanned(bool mIsSC);

         Pt3dr DirRayonR3(const Pt2dr & aPIm) const;
         Pt2di    SzBasicCapt3D() const; 
         double GetVeryRoughInterProf() const;
         bool  CaptHasData(const Pt2dr &) const ;
         Pt2dr    Ter2Capteur   (const Pt3dr & aP) const;
         bool     PIsVisibleInImage   (const Pt3dr & aP,cArgOptionalPIsVisibleInImage * =0) const ;
         ElSeg3D  Capteur2RayTer(const Pt2dr & aP) const;
         double ResolImRefFromCapteur() const ;
         bool  HasRoughCapteur2Terrain() const ;
         Pt2dr ImRef2Capteur   (const Pt2dr & aP) const ;
         bool  HasPreciseCapteur2Terrain() const ;
         Pt3dr RoughCapteur2Terrain   (const Pt2dr & aP) const ;
         Pt3dr PreciseCapteur2Terrain   (const Pt2dr & aP) const ;
         double ResolSolOfPt(const Pt3dr &) const;
         double ResolSolGlob() const;

         double  ScaleAfnt() const;

         Pt3dr Vitesse() const;
         void  SetVitesse(const Pt3dr &);
         bool  VitesseIsInit() const;
         Pt3dr IncCentre() const;
         void  SetIncCentre(const Pt3dr &);

         void TestCam(const std::string & aMes);
         const double & GetTime() const;
         void   SetTime(const double &);
     // ProfIsZ si true, ZProf est l'altisol habituel, sinon c'est une profondeur de champ
         cOrientationConique ExportCalibGlob(Pt2di aSzIm,double AltiSol,double Prof,int AddVerif,bool ModMatr,const char * aNameAux,const Pt3di * aNbVeridDet=0) const;

         cCalibrationInternConique ExportCalibInterne2XmlStruct(Pt2di aSzIm) const;
         // cCalibrationInternConique ExportCalibInterne2XmlStruct(Pt2di aSzIm) const;
         cVerifOrient MakeVerif( int aNbVerif,double aProf,const char *,const Pt3di  * aNbDeterm=0) const;
         cOrientationConique  StdExportCalibGlob(bool Matr) const;
         cOrientationConique  StdExportCalibGlob() const;
         std::string StdExport2File(cInterfChantierNameManipulateur *,const std::string & aDirOri,const std::string & aNameIm);  // Test -> Ori-R

      virtual  Pt3dr ImEtProf2Terrain(const Pt2dr & aP,double aZ) const = 0;
      virtual  Pt3dr NoDistImEtProf2Terrain(const Pt2dr & aP,double aZ) const = 0;
          void SetAltiSol(double );
          void SetProfondeur(double );

           // void ChangeSys(const cSysCoord & a1Source,const cSysCoord & a2Cible,const Pt3dr & aP);
         static void ChangeSys
                     (
                            const std::vector<ElCamera *>& ,
                            const cTransfo3D & aTr3D,
                            bool ForceRot,
                            bool AtGroundLevel
                     );

          // Pour compatibilite stricte avec ce qui etait fait avant
         // dans cDistStdFromCam::Diff
          virtual double SzDiffFinie() const = 0;
          Pt3dr DirVisee() const;
          double ProfondeurDeChamps(const Pt3dr & aP) const;

          virtual double ResolutionSol() const = 0;
          virtual double ResolutionSol(const Pt3dr &) const = 0;
          double GetAltiSol() const;
          bool AltisSolIsDef() const;
          void UndefAltisSol() ;


          double GetProfondeur() const;
          virtual double GetRoughProfondeur() const; // Tente Prof puis Alti
          bool   ProfIsDef() const;
          eTypeProj GetTypeProj() const;
          CamStenope * CS();
          const CamStenope * CS() const;


          virtual cCamStenopeBilin * CSBil_SVP();
          cCamStenopeBilin * CSBil();

         double  RatioInterSol(const ElCamera &) const;

         double   EcartAngulaire(Pt2dr aPF2A, const ElCamera & CamB, Pt2dr aPF2B) const;
         double   SomEcartAngulaire(const ElPackHomologue &, const ElCamera & CamB, double & aSomP) const;
         double   EcartAngulaire(const Appar23 &) const;
         double   SomEcartAngulaire(const std::vector<Appar23> & aVApp) const;
    // L'identifiant ou le nom d'un camera, est qq chose d'optionnel , rajoute a posteriori
    // pour la debugage-tracabilite

         const std::string &  IdentCam() const;
         void SetIdentCam(const std::string & aName);

         const std::string &  NameIm() const;
         void SetNameIm(const std::string & aName);

     // ========================

          ElRotation3D & Orient();
          const ElRotation3D & Orient() const;

          void SetOrientation(const ElRotation3D &);
          void AddToCenterOptical(const Pt3dr & aOffsetC);
          void MultiToRotation(const ElMatrix<double> & aOffsetR);

          Pt3dr  PseudoInter(Pt2dr aPF2A,const ElCamera & CamB,Pt2dr aPF2B,double * aD=0) const;
          // Idem PseudoInter mais la precision est celle de reprojection
          Pt3dr  PseudoInterPixPrec(Pt2dr aPF2A,const ElCamera & CamB,Pt2dr aPF2B,double & aD) const;
          Pt3dr  CdgPseudoInter(const ElPackHomologue &,const ElCamera & CamB,double & aD) const;

          REAL EcartProj(Pt2dr aPF2A,const ElCamera & CamB,Pt2dr aPF2B) const;

          REAL EcartProj(Pt2dr aPF2A,Pt3dr aPR3,Pt3dr aDirR3) const;


          double  ScaleCamNorm() const;
          Pt2dr   TrCamNorm() const;

        //   R3 : "reel" coordonnee initiale
        //   L3 : "Locale", apres rotation
        //   C2 :  camera, avant distortion
        //   F2 : finale apres Distortion
        //
        //       Orientation      Projection      Distortion
        //   R3 -------------> L3------------>C2------------->F2

          Pt2dr R3toF2(Pt3dr) const;
          Pt2dr R3toC2(Pt3dr) const;

          virtual Pt3dr R3toL3(Pt3dr) const;
          virtual Pt3dr L3toR3(Pt3dr) const;

          // Direction en terrain de l'axe camera
          Pt3dr  DirK() const; // OO

          // A la orilib
          Pt3dr F2AndZtoR3(const Pt2dr & aPIm,double aZ) const;

      Pt2dr F2toC2(Pt2dr) const;
      void F2toRayonL3(Pt2dr,Pt3dr &aP0,Pt3dr & aP1) const;
      void F2toRayonR3(Pt2dr,Pt3dr &aP0,Pt3dr & aP1) const;

          Pt3dr PtFromPlanAndIm(const cElPlan3D  & aPlan,const Pt2dr& aP) const;


          ElSeg3D F2toRayonR3(Pt2dr) const;
      Pt3dr   F2toDirRayonL3(Pt2dr) const;
      Pt3dr   F2toDirRayonR3(Pt2dr) const;
      Pt3dr   C2toDirRayonR3(Pt2dr) const;
      Pt2dr   F2toPtDirRayonL3(Pt2dr) const;  // Meme chose, enleve la z a 1
      Pt2dr   L3toF2(Pt3dr) const;
      Pt2dr   PtDirRayonL3toF2(Pt2dr) const;

      Pt2dr Pixel2Radian(const Pt2dr & aP) const;
      Pt2dr Radian2Pixel(const Pt2dr & aP) const;

      Pt3dr   C2toDirRayonL3(Pt2dr) const;
      Pt2dr   L3toC2(Pt3dr) const;

          // Transforme en points photogrammetriques
      ElPackHomologue F2toPtDirRayonL3(const ElPackHomologue &,ElCamera * aCam2=0);  // Def = this
      ElCplePtsHomologues F2toPtDirRayonL3(const ElCplePtsHomologues &,ElCamera * aCam2=0); // Def = this

         Appar23   F2toPtDirRayonL3(const Appar23 &);
     std::list<Appar23>  F2toPtDirRayonL3(const std::list<Appar23>&);

          // Renvoie la somme des ecarts entre la projection des points
          // 3D et les points 2D

          bool Devant(const Pt3dr &) const;
          bool TousDevant(const std::list<Pt3dr> &) const;
          REAL EcProj(const ElSTDNS list<Pt3dr> & PR3 ,
                      const ElSTDNS list<Pt2dr> & PF2) const;

          REAL EcProj ( const ElSTDNS list<Appar23> & P23);

          // Differentielle de l'application globale
                // par rapport a un point
          void  DiffR3F2(ElMatrix<REAL> &,Pt3dr) const;
          ElMatrix<REAL>  DiffR3F2(Pt3dr) const;
                // par rapport aux params
          void  DiffR3F2Param(ElMatrix<REAL> &,Pt3dr) const;
          ElMatrix<REAL>  DiffR3F2Param(Pt3dr) const;

      // void SetDistInverse();
      // void SetDistDirecte();

          bool DistIsDirecte() const;
          bool DistIsC2M() const;
      Pt2dr DistDirecte(Pt2dr aP) const;
      Pt2dr DistInverse(Pt2dr aP) const;
      Pt2dr DistDirecteSsComplem(Pt2dr aP) const;
      Pt2dr DistInverseSsComplem(Pt2dr aP) const;


       // Les tailles representent des capteurs avant Clip et Reech
      const  Pt2di & Sz() const;
          Pt2dr  SzPixel() const;
          Pt2dr  SzPixelBasik() const;
          void  SetSzPixel(const Pt2dr &) ;

      void  SetSz(const Pt2di &aSz,bool AcceptInitMult=false);
          bool SzIsInit() const;

         void SetParamGrid(const cParamForGrid &);
              // AVANT REECH etc... , sz soit etre connu
      void  SetRayonUtile(double aRay,int aNbDisc);

        // La Box utile tient compte d'une eventuelle  affinite
        // elle peut tres bien avoir des coord negatives
           Box2dr BoxUtile() const;

          void HeritComplAndSz(const ElCamera &);
          void CamHeritGen(const ElCamera &,bool WithCompl,bool WithOrientInterne=true);

          void AddCorrecRefrac(cCorrRefracAPost *);
      void AddDistCompl(bool isDirect,ElDistortion22_Gen *);
      void AddDistCompl
           (
               const std::vector<bool> &  isDirect,
               const std::vector<ElDistortion22_Gen *> &
           );
      Pt2dr DComplC2M(Pt2dr,bool UseTrScN = true  ) const;
      Pt2dr DComplM2C(Pt2dr,bool UseTrScN = true  ) const;
          Pt2dr NormC2M(Pt2dr aP) const;
          Pt2dr NormM2C(Pt2dr aP) const;

          ElDistortion22_Gen   &  Get_dist()        ;
          const ElDistortion22_Gen   &  Get_dist() const  ;
      const std::vector<ElDistortion22_Gen *> & DistCompl() const;
      const std::vector<bool> & DistComplIsDir() const;


          // Ajoute une transfo finale pour aller vers la
          // camera, typiquement pour un crop/scale





          // const ElSimilitude & SimM2C();
          static const Pt2di   TheSzUndef ;
          const std::vector<Pt2dr> &  ContourUtile() ;
          bool  HasRayonUtile() const;
          bool IsInZoneUtile(const Pt2dr & aP,bool Pixel=false) const;
          bool     GetZoneUtilInPixel() const;

          double  RayonUtile() const;
     // A priori lie a HasRayonUtile, mais eventuellement
     // autre chose
          bool    HasDomaineSpecial() const;

         virtual ElDistortion22_Gen   *  DistPreCond() const ;
         ElDistortion22_Gen   *  StaticDistPreCond() const ;
  // Eventuellement a redef; now : DistPreCond != 0
         bool IsForteDist() const;


         virtual bool IsGrid() const;
         virtual ~ElCamera();
   // Coincide avec le centre optique pour les camera stenope, est la position
   // du centre origine pour les camera ortho (utilise pour la geom faisceau)
         virtual Pt3dr OrigineProf() const;
         virtual bool  HasOrigineProf() const;
         const cElPolygone &  EmpriseSol() const;
         const Box2dr &  BoxSol() const;

         const tOrIntIma & IntrOrImaC2M() const;


         Pt2dr ResiduMond2Cam(const Pt2dr & aRes)const;
         tOrIntIma  InhibeScaneOri();
         void RestoreScaneOri(const tOrIntIma &);
    protected :







  // Translation et scale de Normalisation
         Pt2dr                          mTrN;
         double                         mScN;


     std::vector<ElDistortion22_Gen *> mDistCompl;
     std::vector<bool>                 mDComplIsDirect;
         cCorrRefracAPost *                mCRAP;

         ElCamera(bool isDistC2M,eTypeProj);
         ElRotation3D     _orient;

         virtual       ElProj32 &        Proj()       = 0;
         virtual const ElProj32       &  Proj() const = 0;
     Pt2di    mSz;
         Pt2dr    mSzPixel;



     // Une distorsion de "pre-conditionnement" est une fonction "simple"
     // qui approxime la partie non lineaire de la distorsion, si !=0 elle
     // est exprimee dans le sens M->C , 0 signifie identite
     //
     // Elle est utilisee notamment parce que les distorsions "compliquees"
     // peuvent etre exprimees comme la composition d'une distorsion
     // grille a faible distorsion de la distorsion de "pre-conditionnement"
   protected :
     bool             mDIsDirect;
   public :
         virtual       ElDistortion22_Gen   &  Dist()        = 0;
   protected :
         virtual void InstanceModifParam(cCalibrationInternConique &) const  =0;
         virtual const ElDistortion22_Gen   &  Dist() const  = 0;

         void AssertSolInit() const;


         eTypeProj   mTypeProj;
   protected :
         bool        mAltisSolIsDef;
         double      mAltiSol;
         bool        mProfondeurIsDef;
         double      mProfondeur;

   private :

         std::string  mIdentCam;
         std::string  mNameIm;

         //double      mPrecisionEmpriseSol;
         cElPolygone mEmpriseSol;
         Box2dr      mBoxSol;

         double              mRayonUtile;
         bool                mHasDomaineSpecial;
         bool                mDoneScanContU;
         std::vector<Pt2dr>  mContourUtile;

   protected :
         bool                 mParamGridIsInit;
         Pt2dr                mStepGrid;
         double               mRayonInvGrid;
         double               mTime;
         bool                 mScanned;

   private :
         Pt3dr  mVitesse;
         bool   mVitesseIsInit;
         Pt3dr  mIncCentre;

         mutable ElDistortion22_Gen *mStatDPC;
         mutable bool                mStatDPCDone;
};




class CamStenope : public ElCamera
{
      public :
         CamStenope * DownCastCS() ;
         virtual std::string Save2XmlStdMMName(  cInterfChantierNameManipulateur * anICNM,
                                        const std::string & aOriOut,
                                        const std::string & aNameImClip,
                                        const ElAffin2D & anOrIntInit2Cur
                    ) const;

         double GetRoughProfondeur() const; // Tente Prof puis Alti
         const tParamAFocal   & ParamAF() const;

         void StdNormalise(bool doScale,bool  doTr);
         void StdNormalise(bool doScale,bool  doTr,double aS,Pt2dr  aTr);
         void UnNormalize();
         // .xml ou .ori
         static CamStenope * StdCamFromFile(bool UseGr,const std::string &,cInterfChantierNameManipulateur * anICNM);

         virtual const cCamStenopeDistRadPol * Debug_CSDRP() const;



          // renvoit la distance de p1 a la projection de la droite
          // Inhibee car non testee

          // La methode a ete definie dans la mere, il n'y a aucun interet
          // apparement a la specialiser

          // REAL EcartProj(Pt2dr aPF2A,const ElCamera & CamB,Pt2dr aPF2B);

         // Helas, le SzIm n'est pas integre dans mes CamStenope ...





         CamStenope(bool isDistC2M,REAL Focale,Pt2dr centre,const std::vector<double>  & AFocalParam);
         CamStenope(const CamStenope &,const ElRotation3D &);

         // Par defaut true, mais peut redefini, par exemple pour
         // un fish-eye
         virtual bool CanExportDistAsGrid() const;

         void OrientFromPtsAppui
              (
                 ElSTDNS list<ElRotation3D> &,
                 Pt3dr R3A, Pt3dr R3B, Pt3dr R3C,
                 Pt2dr F2A, Pt2dr F2B, Pt2dr F2C
              );
         void OrientFromPtsAppui
              (
                 ElSTDNS list<ElRotation3D> & Res,
                 const ElSTDNS list<Pt3dr> & PR3 ,
                 const ElSTDNS list<Pt2dr> & PF2
              );
         void OrientFromPtsAppui
              (
                 ElSTDNS list<ElRotation3D>  & Res,
                 const ElSTDNS list<Appar23> & P32
              );

        // Si  NbSol ==  0 et resultat vide => Erreur
        // Sinon *NbSol Contient  le nombre de solution

         ElRotation3D  OrientFromPtsAppui
              (
                 bool TousDevant,
                 const ElSTDNS list<Pt3dr> & PR3 ,
                 const ElSTDNS list<Pt2dr> & PF2 ,
                 REAL * Ecart = 0,
                 INT  * NbSol    = 0
              );

         ElRotation3D  OrientFromPtsAppui
              (
                                bool TousDevant,
                 const ElSTDNS list<Appar23> & P32 ,
                 REAL * Ecart = 0,
                 INT  * NbSol    = 0
              );
     ElRotation3D  CombinatoireOFPAGen
               (
                                bool TousDevant,
                INT  NbTest,
                const ElSTDNS list<Pt3dr> & PR3 ,
                const ElSTDNS list<Pt2dr> & PF2,
                REAL * Res_Dmin,
                bool   ModeRansac,
                                Pt3dr * aDirApprox = 0
                       );

     ElRotation3D  CombinatoireOFPA
               (
                                bool TousDevant,
                INT  NbTest,
                const ElSTDNS list<Pt3dr> & PR3 ,
                const ElSTDNS list<Pt2dr> & PF2,
                REAL * Res_Dmin,
                                Pt3dr * aDirApprox = 0
               );

     ElRotation3D  RansacOFPA
               (
                                bool TousDevant,
                INT  NbTest,
                const ElSTDNS list<Appar23> & P23 ,
                REAL * Res_Dmin,
                                Pt3dr * aDirApprox = 0
               );



     ElRotation3D  CombinatoireOFPA
               (
                                bool TousDevant,
                INT  NbTest,
                                const ElSTDNS list<Appar23> & P32 ,
                REAL * Res_Dmin,
                                Pt3dr * aDirApprox = 0
               );


         // Orientations avec "GPS", i.e. avec centre fixe

         void Set_GPS_Orientation_From_Appuis
                      (
                           const Pt3dr & aGPS,
                           const std::vector<Appar23> & aVApp,
                           int  aNbRansac
                      );

         // Pour compatibilite temporaire avec la proj carto d'orilib
         virtual Ori3D_Std * CastOliLib();  // OO  Def return 0
         Ori3D_Std * NN_CastOliLib();  //OO   Erreur si 0
         double ResolutionPDVVerticale();  //OO   OriLib::resolution, assume implicitement une
                                           // PDV sub verticale
         double ResolutionAngulaire() const;  // OO
         double ResolutionSol() const ;
         double ResolutionSol(const Pt3dr &) const ;
         // Pour l'instant bovin, passe par le xml
         virtual CamStenope * Dupl() const;   // OO


     REAL Focale() const ;
     Pt2dr PP() const ;
     Pt3dr VraiOpticalCenter() const;
     Pt3dr PseudoOpticalCenter() const;
     Pt3dr    OpticalCenterOfPixel(const Pt2dr & aP) const ; 
     Pt3dr OpticalVarCenterIm(const Pt2dr &) const;
     Pt3dr OpticalVarCenterTer(const Pt3dr &) const;
     Pt3dr ImEtProf2Terrain(const Pt2dr & aP,double aZ) const;
     Pt3dr NoDistImEtProf2Terrain(const Pt2dr & aP,double aZ) const;
     Pt3dr ImEtZ2Terrain(const Pt2dr & aP,double aZ) const;
     void  Coins(Pt3dr &aP1, Pt3dr &aP2, Pt3dr &aP3, Pt3dr &aP4, double aZ) const;
     void  CoinsProjZ(Pt3dr &aP1, Pt3dr &aP2, Pt3dr &aP3, Pt3dr &aP4, double aZ) const;
     Box2dr BoxTer(double aZ) const;

         Pt3dr  ImEtProfSpherik2Terrain(const Pt2dr & aPIm,const REAL & aProf) const; //OO
         Pt3dr  ImDirEtProf2Terrain(const Pt2dr & aPIm,const REAL & aProf,const Pt3dr & aNormPl) const; //OO
         Pt3dr Im1DirEtProf2_To_Terrain  //OO
               (Pt2dr p1,const CamStenope &  ph2,double prof2,const Pt3dr & aDir) const;
         Pt3dr Im1EtProfSpherik2_To_Terrain (Pt2dr p1,const CamStenope &  ph2,double prof2) const;
    void ExpImp2Bundle(const Pt2di aGridSz, const std::string aName) const;

     double ProfInDir(const Pt3dr & aP,const Pt3dr &) const; // OO


         // Sert pour un clonage, par defaut null
         virtual ElProj32             &  Proj();
         virtual const ElProj32       &  Proj() const;
         virtual ElDistortion22_Gen   &  Dist();
         virtual const ElDistortion22_Gen   &  Dist() const;

// Def  : erreur fatale
         virtual cParamIntrinsequeFormel * AllocParamInc(bool isDC2M,cSetEqFormelles &);


         cCamStenopeDistRadPol *Change2Format_DRP
                            (
                      bool C2M,
                      int  aDegreOut,
                      bool CDistPPLie,
                      double Resol,
                      Pt2dr  Origine
                );




         void InstanceModifParam(cCalibrationInternConique &)  const;
         Pt3dr OrigineProf() const;
         bool  HasOrigineProf() const;
         bool  UseAFocal() const;
      private :
         CamStenope(const CamStenope &); // N.I.

      protected :
         ElProjStenope  _PrSten;
         bool                 mUseAF;
         ElDistortion22_Gen * mDist;

         double SzDiffFinie() const;
};





#endif
