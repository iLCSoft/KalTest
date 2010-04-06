#ifndef __EXKALSTATE__
#define __EXKALSTATE__
//*************************************************************************
//* ====================
//*  EXKalState Class
//* ====================
//*
//* (Description)
//*   Sample state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class EXKalState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVKalState.h"

class EXKalSite;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//

class EXKalState : public TVKalState {
                                                                                
public:
                                                                                
   // Ctors and Dtor
                                                                                
   EXKalState(Int_t p = 2);
   EXKalState(const TKalMatrix &sv, Int_t type = 0, Int_t p = 2);
   EXKalState(const TKalMatrix &sv, const TKalMatrix &c, 
              Int_t type = 0, Int_t p = 2);
   EXKalState(const TKalMatrix &sv, const TVKalSite &site, 
              Int_t type = 0, Int_t p = 2);
   EXKalState(const TKalMatrix &sv, const TKalMatrix &c,
              const TVKalSite &site, Int_t type = 0, Int_t p = 2);
   virtual ~EXKalState() {}
                                                                                
   // Implementation of paraent class pure virtuals
                                                                                
   EXKalState * MoveTo(TVKalSite  &to, 
                       TKalMatrix &F, 
                       TKalMatrix *QPtr = 0) const;
   EXKalState & MoveTo(TVKalSite  &to, 
                       TKalMatrix &F, 
                       TKalMatrix &Q) const;
   void         DebugPrint() const;

   // Derived class methods

   TKalMatrix CalcXAt(Double_t t);
   TKalMatrix CalcXDerivAt(Double_t t);
                                                                                
                                                                                
   ClassDef(EXKalState,1)      // sample state vector class
};
                                                                                
//=======================================================
// inline functions
//=======================================================
                                                                                
#endif
