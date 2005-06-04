#ifndef __TKALTRACKSTATE__
#define __TKALTRACKSTATE__
//*************************************************************************
//* ======================
//*  TKalTrackState Class
//* ======================
//*
//* (Description)
//*   Track state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class TKalTrackState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added CalcDapDa method.
//*
//*************************************************************************
//
#include "TVKalState.h"
#include "THelicalTrack.h"
#include "TStraightTrack.h"
#include "KalTrackDim.h"

class TKalTrackSite;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//

class TKalTrackState : public TVKalState {
                                                                                
public:
                                                                                
   // Ctors and Dtor
                                                                                
   TKalTrackState(Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c, 
                        Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TVKalSite &site, 
                        Int_t type = 0, Int_t p = kSdim);
   TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                  const TVKalSite &site, Int_t type = 0, Int_t p = kSdim);
   virtual ~TKalTrackState() {}
                                                                                
   // Implementation of paraent class pure virtuals
                                                                                
   TKalTrackState * MoveTo(const TVKalSite  &to, 
                                 TKalMatrix &F, 
                                 TKalMatrix *QPtr = 0) const;
   TKalTrackState & MoveTo(const TVKalSite  &to, 
                                 TKalMatrix &F, 
                                 TKalMatrix &Q) const;
   void         DebugPrint() const;

   // Derived class methods

   THelicalTrack   GetHelix() const;
   TStraightTrack  GetLine () const;
   TVTrack        &CreateTrack() const;

private:
   TKalMatrix      CalcProcessNoise(const TKalTrackSite  &to,
                                          TKalTrackState &ato,
                                    const TVTrack        &tto,
                                          Double_t        dfi)   const;

private:

   TVector3 fX0;		// pivot

   ClassDef(TKalTrackState,1)      // sample state vector class
};
                                                                                
//=======================================================
// inline functions
//=======================================================
                                                                                
#endif
