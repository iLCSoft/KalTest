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

#include "EXVMeasLayer.h"

Bool_t   EXVMeasLayer::kActive = kTRUE;
Bool_t   EXVMeasLayer::kDummy = kFALSE;

ClassImp(EXVMeasLayer)
                                                                                
EXVMeasLayer::EXVMeasLayer(TMaterial &min,
                           TMaterial &mout,
                           Bool_t     isactive,
                     const Char_t    *name)  
            : TVMeasLayer(min, mout, isactive),
	      fName(name),
	      fNodePtr(0)
{
}

EXVMeasLayer::~EXVMeasLayer()
{
}
