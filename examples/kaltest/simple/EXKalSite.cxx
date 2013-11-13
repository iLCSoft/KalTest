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
#include <iostream>
#include "EXKalSite.h"
#include "EXKalState.h"
#include "EXHit.h"
#include "TMath.h"

using namespace std;

//_____________________________________________________________________
//  ----------------------------------
//  Sample class for measurement site
//  ----------------------------------
//
ClassImp(EXKalSite)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------
EXKalSite::EXKalSite(Int_t m, Int_t p)
                   :TVKalSite(m,p), fHitPtr(0)
{
}

EXKalSite::EXKalSite(const EXHit & ht, Int_t m, Int_t p)
                   :TVKalSite(m,p), fHitPtr((EXHit *)&ht)
{
   for (Int_t i=0; i<m; i++) {
      GetMeasVec     ()(i,0) = ht.GetX(i);
      GetMeasNoiseMat()(i,i) = TMath::Power(ht.GetDX(i),2);
   }
}
                                                                                
//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods
//  ----------------------------------

TVKalState & EXKalSite::CreateState(const TKalMatrix &sv, Int_t type)
{
   return *(new EXKalState(sv,type));
}

TVKalState & EXKalSite::CreateState(const TKalMatrix &sv, const TKalMatrix &c,
                                    Int_t type)
{
   return *(new EXKalState(sv,c,type));
}

Int_t EXKalSite::CalcExpectedMeasVec(const TVKalState &a, TKalMatrix &h)
{
   h = (*(EXKalState *)&a).CalcXAt(fHitPtr->GetT());
   return 1;
}

Int_t EXKalSite::CalcMeasVecDerivative(const TVKalState &a, TKalMatrix &H)
{
   H = (*(EXKalState *)&a).CalcXDerivAt(fHitPtr->GetT());
   return 1;
}

Bool_t EXKalSite::IsAccepted()
{
   // return kTRUE if this site is accepted by Filter()
   return kTRUE;
}

void EXKalSite::DebugPrint() const
{
   fHitPtr->DebugPrint();
   cerr << " DeltaChi2 = " << GetDeltaChi2() << endl;
}

