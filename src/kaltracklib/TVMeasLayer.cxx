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
//*   2005/0X/XX  A.Yamaguchi       Replaced fFwdX0Inv, fBwdX0Inv, and
//*                                 their getters and setters by
//*                                 fMaterialOutPtr, fMaterialInPtr, and
//*                                 their getters and setters.
//*   2005/08/15  K.Fujii           Added fIsActive and IsActive().
//*                                 Moved GetEnergyLoss() and CalcQms()
//*                                 from TVKalDetector.
//*   2011/12/03  S.Aplin           Added new member: name 
//*                                 default value set to "TVMeasLayer"
//*                                 and corresponding member function
//*                                 TString GetName()
//*   2012/11/29  K.Fujii           Removed fgCurInstancePtr and its
//*                                 getter & setter.
//*                                 Set parent pointer in Add() instead.
//*
//*************************************************************************

#include "TVMeasLayer.h"  // from KalTrackLib
#include "TKalTrack.h"    // from KalTrackLib
#include "TVTrack.h"      // from KalTrackLib

ClassImp(TVMeasLayer)

//_________________________________________________________________________
// ----------------------------------
//  Ctor
// ----------------------------------
TVMeasLayer::TVMeasLayer(TMaterial     &matIn,
                         TMaterial     &matOut,
                         Bool_t         isactive,
                         const Char_t    *name)
           : fMaterialInPtr(&matIn),
             fMaterialOutPtr(&matOut),
             fIndex(0),
             fIsActive(isactive),
             fname(name)
{
}
