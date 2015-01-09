#ifndef __EXITFBMEASLAYER__
#define __EXITFBMEASLAYER__
//*************************************************************************
//* ===================
//*  EXITFBMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXITHit.
//* (Requires)
//* (Provides)
//*     class EXITFBMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/07/25  Kim, Youngim      Forward & Backward version.
//*
//*   2011/06/30  D.Kamai       Modified to handle turbine-blade-like FTD.
//*************************************************************************
//
#include "TVector3.h"
#include "TPlane.h"
#include "TKalMatrix.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"

class TVTrackHit;

class EXITFBMeasLayer : public EXVMeasLayer, public TPlane {
public:
   // Ctors and Dtor

   EXITFBMeasLayer(TMaterial &min,
                   TMaterial &mout,
                   const TVector3  &xc,
                   const TVector3  &normal,
                   Double_t   SortingPolicy,
                   Double_t   rin,
                   Double_t   rout,
                   Double_t   dxMax,
                   Double_t   dxMin,
                   Double_t   sigmax,
                   Double_t   sigmar,
                   Bool_t     type = EXVMeasLayer::kActive,
                   Int_t      mode = 0,
                   const Char_t    *name = "FTDML");

   virtual ~EXITFBMeasLayer();

   // Parrent's pure virtuals that must be implemented
   
   virtual TKalMatrix XvToMv           (const TVTrackHit &ht,
                                        const TVector3   &xv) const;
   virtual TKalMatrix XvToMv           (const TVector3   &xv) const;
   virtual TVector3   HitToXv          (const TVTrackHit &ht) const;
   virtual void       CalcDhDa         (const TVTrackHit &ht,
                                        const TVector3   &xv,
                                        const TKalMatrix &dxphiada,
                                              TKalMatrix &H)  const;

   
   virtual Bool_t     IsOnSurface       (const TVector3 &xx) const;

   virtual void       ProcessHit        (const TVector3  &xx,
                                         TObjArray &hits);

   Double_t GetSortingPolicy() const {return fSortingPolicy;}
   Double_t GetRin() const {return fRin;}
   Double_t GetRout() const {return fRout;}
   Double_t GetdxMax() const {return fdxMax;}
   Double_t GetdxMin() const {return fdxMin;}
   Double_t GetSigmaX() const { return fSigmaX; } // by yikim 2005/07/28
   Double_t GetSigmaY() const { return fSigmaY; }
   Int_t    GetMode() const { return fMode; }
   
   inline Double_t Cosphi() const {return GetXc().Y()/GetXc().Perp(); }
   inline Double_t Sinphi() const {return GetXc().X()/GetXc().Perp(); }
   inline Double_t Cosalpha() const {return TMath::Abs(GetNormal().Z()/GetNormal().Mag()); }
   inline Double_t Sinalpha() const {return -TMath::Abs(GetNormal().X()/Cosphi()); }
   
   using TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt);

private:
   Double_t fSortingPolicy;
   Double_t fRin;     // inner radius
   Double_t fRout;    // outer radius
   Double_t fdxMax;   // length of trapezoid bottom base
   Double_t fdxMin;   // length of trapozoid upper base
   Double_t fSigmaX;  // sigma_x
   Double_t fSigmaY;  // sigma_y
   Int_t    fMode;    // Right-side, Left-side, or Both
   
   ClassDef(EXITFBMeasLayer,1)  // Sample measurement layer class
};

#endif
