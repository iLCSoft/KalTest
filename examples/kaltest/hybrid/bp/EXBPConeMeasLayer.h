#ifndef EXBPCONEMEASLAYER_H
#define EXBPCONEMEASLAYER_H
//*************************************************************************
//* ===================
//*  EXBPConeMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXBPConeHit.
//* (Requires)
//* (Provides)
//*     class EXBPConeMeasLayer
//* (Update Recored)
//*   2012/01/19  K.Fujii       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCutCone.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXBPConeMeasLayer : public EXVMeasLayer, public TCutCone {
public:
   // Ctors and Dtor

   EXBPConeMeasLayer(TMaterial &min,
                     TMaterial &mout,
                     Double_t   z1,
                     Double_t   r1,
                     Double_t   z2,
                     Double_t   r2,
                     Double_t   sigmax,
                     Double_t   sigmaz,
                     Bool_t     type = EXVMeasLayer::kDummy,
               const Char_t    *name = "BPCONEML");
   virtual ~EXBPConeMeasLayer();

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;

   // Parent's virtuals overridden
   virtual Bool_t     IsOnSurface(const TVector3 &xx) const;


   virtual void       ProcessHit(const TVector3   &xx,
                                       TObjArray  &hits);


   Double_t GetSigmaX() const { return fSigmaX; }
   Double_t GetSigmaZ() const { return fSigmaZ; }

   using TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt);

private:
   Double_t fZ1;      // z of front face
   Double_t fR1;      // r of front face
   Double_t fZ2;      // z of back end
   Double_t fR2;      // r of back end
   Double_t fSigmaX;  // sigma_x
   Double_t fSigmaZ;  // sigma_z

   ClassDef(EXBPConeMeasLayer,1)   // Sample measurement layer class
};

#endif
