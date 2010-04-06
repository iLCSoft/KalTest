#ifndef __EXTPCMEASLAYER__
#define __EXTPCMEASLAYER__
//*************************************************************************
//* ===================
//*  EXTPCMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXTPCHit.
//* (Requires)
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXTPCMeasLayer : public EXVMeasLayer, public TCylinder {
public:
   // Ctors and Dtor

   EXTPCMeasLayer(TMaterial &min,
                  TMaterial &mout,
                  Double_t   r0,
                  Double_t   lhalf,
                  Double_t   sigmax0,
                  Double_t   sigmax1,
                  Double_t   sigmaz,
                  Bool_t     type = EXVMeasLayer::kActive,
            const Char_t     *name = "TPCML");
   virtual ~EXTPCMeasLayer();

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv)   const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv,
                                       Int_t       side) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)    const;
   virtual void       ProcessHit(const TVector3   &xx,
                                       TObjArray  &hits);

   Double_t GetSigmaX(Double_t z) const;
   Double_t GetSigmaZ()           const { return fSigmaZ; }

private:
   Double_t fSigmaX0;   // xy resolution
   Double_t fSigmaX1;   // xy resolution
   Double_t fSigmaZ;    // z  resolution

   ClassDef(EXTPCMeasLayer,1) 	// Sample measurement layer class
};

#endif
