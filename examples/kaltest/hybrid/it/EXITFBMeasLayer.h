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
//*************************************************************************
//
#include "TVector3.h"
#include "TPlane.h"
#include "TKalMatrix.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXITFBMeasLayer : public EXVMeasLayer, public TPlane {
public:
   // Ctors and Dtor

   EXITFBMeasLayer(TMaterial &min,
                  TMaterial  &mout,
                  TVector3   &xc,
                  Double_t    rin,
                  Double_t    rout,
                  Double_t    sigmax,
                  Double_t    sigmar,
                  Bool_t      type = EXVMeasLayer::kActive);
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
   virtual Int_t      CalcXingPointWith(const TVTrack  &hel,
                                              TVector3 &xx,
                                              Double_t &phi,
                                              Double_t  eps = 1.e-8) const;

   inline virtual Double_t GetSortingPolicy ()            const;
   inline virtual Bool_t IsOnSurface (const TVector3 &xx) const;

   inline Double_t GetSigmaX() const { return fSigmaX; } // by yikim 2005/07/28
   inline Double_t GetSigmaY() const { return fSigmaY; }

private:
   Double_t fRin;
   Double_t fRout;
   Double_t fSigmaX;
   Double_t fSigmaY;

   ClassDef(EXITFBMeasLayer,1) 	// Sample measurement layer class
};

Double_t EXITFBMeasLayer::GetSortingPolicy() const
{
   return IsActive() ? fRout : fRout + 1.e-10;
} 

Bool_t EXITFBMeasLayer::IsOnSurface(const TVector3 &xx) const
{
   return xx.Pt() < fRout && xx.Pt() > fRin ? kTRUE : kFALSE; 
}

#endif
