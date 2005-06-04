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
                           Double_t   r0,
                           Double_t   lhalf,
                           Bool_t     type)
            : TVMeasLayer(min, mout),
              TCylinder(r0, lhalf),
              fType(type)
{
}

EXVMeasLayer::~EXVMeasLayer()
{
}
