#ifndef __EXMEASLAYER__
#define __EXMEASLAYER__
//*************************************************************************
//* ===================
//*  EXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXHit.
//* (Requires)
//* (Provides)
//*     class EXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "TVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXMeasLayer : public TVMeasLayer, public TCylinder {
public:
   static Bool_t kActive;
   static Bool_t kDummy;

   // Ctors and Dtor

   EXMeasLayer(TMaterial &min,
               TMaterial &mout,
               Double_t   r0,
               Double_t   lhalf,
               Bool_t     type = EXMeasLayer::kActive);
   virtual ~EXMeasLayer();

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;

   // Methods for MC event generation
   virtual void       ProcessHit(const TVector3    &xx,
                                       TObjArray   &hits);

   Double_t GetSigmaX() const { return fgSigmaX; }
   Double_t GetSigmaZ() const { return fgSigmaZ; }

private:
   static Double_t fgSigmaX;   // xy resolution
   static Double_t fgSigmaZ;   // z  resolution

   ClassDef(EXMeasLayer,1) 	// Sample measurement layer class
};

#endif
