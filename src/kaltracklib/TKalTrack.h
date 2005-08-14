#ifndef __TKALTRACK__
#define __TKALTRACK__
//*************************************************************************
//* =================
//*  TKalTrack Class
//* =================
//*
//* (Description)
//*   Track class for Kalman filter
//* (Requires)
//*     TVKalSystem
//* (Provides)
//*     class TKalTrack
//* (Update Recored)
//*   2003/09/30  Y.Nakashima	Original version.
//*   2005/02/23  A.Yamaguchi	Added a new data member fMass and its
//*                             getter and setter.
//*
//*************************************************************************
                                                                                
#include "TVKalSystem.h"       // from KalLib
#include "TKalTrackState.h"    // from KalTrackLib

//_________________________________________________________________________
//  ------------------------------
//   TKalTrack: Kalman Track class
//  ------------------------------
                                                                                
class TKalTrack : public TVKalSystem {
public:
   TKalTrack(Int_t n = 1);
   ~TKalTrack() {} 

   inline virtual void      SetFitDirection(Bool_t isfwd) { fDir = isfwd; }
   inline virtual void      SetMass(Double_t m)           { fMass = m;    }
   inline virtual Bool_t    GetFitDirection()     const   { return fDir;  }
   inline virtual Double_t  GetMass()             const   { return fMass; }

   Double_t FitToHelix(TKalTrackState &a, TKalMatrix &C, Int_t &ndf);

private:
   Double_t     fMass;        // mass [GeV]
   Bool_t       fDir;         // fitting direction

   ClassDef(TKalTrack,1)  // Base class for Kalman Filter
};

#endif
