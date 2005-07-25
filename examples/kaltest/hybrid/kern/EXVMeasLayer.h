#ifndef __EXVMEASLAYER__
#define __EXVMEASLAYER__
//*************************************************************************
//* ===================
//*  EXVMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by TVTrackHit.
//* (Requires)
//* (Provides)
//*     class EXVMeasLayer
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

class EXVMeasLayer : public TVMeasLayer {
public:
   static Bool_t kActive;
   static Bool_t kDummy;

   // Ctors and Dtor

   EXVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Bool_t     type = EXVMeasLayer::kActive);
   virtual ~EXVMeasLayer();

   // Parrent's pure virtuals that must be implemented

   Bool_t IsActive () const { return fType;   }

private:
   Bool_t   fType;      // (true, false) = (active layer, dummy layer)

   ClassDef(EXVMeasLayer,1) 	// Sample measurement layer class
};

#endif
