#ifndef __TVMEASLAYER__
#define __TVMEASLAYER__
//*************************************************************************
//* ====================
//*  TVMeasLayer Class
//* ====================
//*
//* (Description)
//*   Measurement layer interface class.
//* (Requires)
//* (Provides)
//*     class TVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added new data members, fFwdX0Inv,
//*                                 fBwdX0Inv and fIndex, and their
//*                                 corresponding getters and setters.
//*                                 Added a new method, GetX0Inv().
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TAttElement.h"
#include "TKalMatrix.h"
#include "KalTrackDim.h"
#include "TMaterial.h"

class TVTrackHit;

class TVMeasLayer : public TAttElement {
public:
   // Ctors and Dtor

   TVMeasLayer(TMaterial &matIn, TMaterial &matOut);
   virtual ~TVMeasLayer() {}

   // Utiliy Methods

   virtual TKalMatrix XvToMv   (const TVTrackHit &ht,
                                const TVector3   &xv) const = 0;
   virtual TVector3   HitToXv  (const TVTrackHit &ht) const = 0;
   virtual void       CalcDhDa (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                      TKalMatrix &H)  const = 0;

   inline virtual TMaterial &GetMaterial(Bool_t isfwd) const
   { return isfwd ? *fMaterialOutPtr : *fMaterialInPtr; }
     
   inline Int_t GetIndex() const  { return fIndex; }
   inline void  SetIndex(Int_t i) { fIndex = i;    } 

private:
   TMaterial     *fMaterialInPtr;   // pointer of inner Material
   TMaterial     *fMaterialOutPtr;  // pointer of outer Material
   Int_t          fIndex;           // index in TKalDetCradle

   ClassDef(TVMeasLayer,1) 	// Measurement layer interface class
};

#endif
