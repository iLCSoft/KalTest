#ifndef EXITMEASLAYER_H
#define EXITMEASLAYER_H
//*************************************************************************
//* ===================
//*  EXITMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXITHit.
//* (Requires)
//* (Provides)
//*     class EXITMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/07/25  Kim, Youngim 
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXITMeasLayer : public EXVMeasLayer , public TCylinder {
public:
   // Ctors and Dtor

   EXITMeasLayer(TMaterial &min,
                 TMaterial &mout,
                 Double_t   r0,
                 Double_t   lhalf,
                 Double_t   sigmax,
                 Double_t   sigmaz,
                 Bool_t     type = EXVMeasLayer::kActive,
           const Char_t    *name = "ITML");
   virtual ~EXITMeasLayer();

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;
   virtual void       ProcessHit(const TVector3   &xx,
                                       TObjArray  &hits);

   Double_t GetSigmaX() const { return fSigmaX; }
   Double_t GetSigmaZ() const { return fSigmaZ; }

   using TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt);

private:
   Double_t fSigmaX;      // sigma_x
   Double_t fSigmaZ;      // sigma_z

   ClassDef(EXITMeasLayer,1)    // Sample measurement layer class
};

#endif
