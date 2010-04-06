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
#include <iostream>
#include "EXKalState.h"
#include "EXKalSite.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//
// --------------------------------------------------------------------
// Ctors and Dtor
//
EXKalState::EXKalState(Int_t p) 
           : TVKalState(p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, Int_t type, Int_t p) 
           : TVKalState(sv,type,p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       Int_t type, Int_t p) 
           : TVKalState(sv,c,type,p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, const TVKalSite &site,
                       Int_t type, Int_t p) 
           : TVKalState(sv,site,type,p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       const TVKalSite &site, Int_t type, Int_t p) 
           : TVKalState(sv,c,site,type,p)
{
}

//
// --------------------------------------------------------------------
// Implementation of base-class pure virtuals
//

EXKalState * EXKalState::MoveTo(TVKalSite  &to,
                                TKalMatrix &F,
                                TKalMatrix *QPtr) const
{
   Int_t p = GetDimension();
   for (Int_t i=0; i<p; i++) F(i,i) = 1.;
   if (QPtr) {
      QPtr->Zero();
      return (new EXKalState(*this,to,TVKalSite::kPredicted));
   } else {
      return 0;
   }
}

EXKalState & EXKalState::MoveTo(TVKalSite  &to,
                                TKalMatrix &F,
                                TKalMatrix &Q) const
{
   return *MoveTo(to, F, &Q);
}

void EXKalState::DebugPrint() const
{
   cout << " slope     = " << (*this)(0,0) << endl
        << " intercept = " << (*this)(1,0) << endl;
   GetCovMat().DebugPrint(" cov Mat   = ");
}

//
// --------------------------------------------------------------------
// Derived class methods
//

TKalMatrix EXKalState::CalcXAt(Double_t t)
{
   TKalMatrix m(1,1);
   m(0,0) = (*this)(0,0)*t + (*this)(1,0);
   return m;
}

TKalMatrix EXKalState::CalcXDerivAt(Double_t t)
{
   TKalMatrix H(1,2);
   H(0,0) = t;
   H(0,1) = 1.;
   return H;
}
