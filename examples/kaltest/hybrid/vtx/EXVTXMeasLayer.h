#ifndef __EXVTXMEASLAYER__
#define __EXVTXMEASLAYER__
//*************************************************************************
//* ===================
//*  EXVTXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXVTXHit.
//* (Requires)
//* (Provides)
//*     class EXVTXMeasLayer
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

class EXVTXMeasLayer : public EXVMeasLayer, public TCylinder {
public:
   // Ctors and Dtor

   EXVTXMeasLayer(TMaterial &min,
                  TMaterial &mout,
                  Double_t   r0,
                  Double_t   lhalf,
                  Double_t   sigmax,
                  Double_t   sigmaz,
                  Bool_t     type = EXVMeasLayer::kActive,
            const Char_t    *name = "VTXML");
   virtual ~EXVTXMeasLayer();

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
   Double_t fSigmaX;  // sigma_x
   Double_t fSigmaZ;  // sigma_z

   ClassDef(EXVTXMeasLayer,1) 	// Sample measurement layer class
};

#endif
