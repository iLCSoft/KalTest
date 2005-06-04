#ifndef __TVKALDETECTOR__
#define __TVKALDETECTOR__
//*************************************************************************
//* =====================
//*  TVKalDetector Class
//* =====================
//*
//* (Description)
//*   Base class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalDetector
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2005/02/23  A.Yamaguchi	Moved most methods to TKalDetCradle.
//*
//*************************************************************************

#include "TObjArray.h"
#include "TAttElement.h"
#include "TVMeasLayer.h"
#include "TKalTrackState.h"

//_____________________________________________________________________
//  ------------------------------
//  Detector system class
//  ------------------------------
//

class TVKalDetector : public TObjArray, public TAttElement {
public:
   // Ctors and Dtor
   TVKalDetector(Int_t n = 1) : TObjArray(n) {}
   virtual ~TVKalDetector() {}

   virtual void       CalcEnergyLoss(const TMaterial      &mat,
                                           Double_t        df,
                                     const TKalTrackState &afrom,
                                           TKalTrackState &ato)  const;
   virtual TKalMatrix CalcSigmaMS0  (const TMaterial      &mat, 
                                           Double_t        df,
                                     const TKalTrackState &from) const;
   
   ClassDef(TVKalDetector,1)  // Base class for detector system
};

#endif
