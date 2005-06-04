#ifndef __EXKALSITE__
#define __EXKALSITE__
//*************************************************************************
//* ===================
//*  EXKalSite Class
//* ===================
//*
//* (Description)
//*   Sample measurement site class used by Kalman filter.
//* (Requires)
//*     EXKalState
//* (Provides)
//*     class EXKalSite
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVKalSite.h"
class EXHit;
class TVKalState;
class EXKalState;

//_____________________________________________________________________
//  ---------------------------------
//  Class for Kalman measurement site
//  ---------------------------------
//
class EXKalSite : public TVKalSite {
public:
   EXKalSite(Int_t m = 1, Int_t p = 2);
   EXKalSite(const EXHit &ht, Int_t m = 1, Int_t p = 2);
   ~EXKalSite() {}

   Int_t   CalcExpectedMeasVec(const TVKalState &a, TKalMatrix &h);
   Int_t   CalcMeasVecDerivative(const TVKalState &a, TKalMatrix &H);
   Bool_t  IsAccepted();

   void  DebugPrint() const;

private:
   TVKalState & CreateState(const TKalMatrix &sv, Int_t type = 0);
   TVKalState & CreateState(const TKalMatrix &sv, const TKalMatrix &C,
                            Int_t type = 0);

private:
   EXHit *fHitPtr;   // pointer to corresponding hit

   ClassDef(EXKalSite,1)  // sample measurement site class
};
                                                                                
//=======================================================
// inline functions
//=======================================================
                                                                                
#endif
