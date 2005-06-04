//*************************************************************************
//* ====================
//*  TVKalSite Class
//* ====================
//*
//* (Description)
//*   This is the base class for measurement vector used by Kalman filter.
//* (Requires)
//* 	TKalMatrix
//*   TVKalState
//* (Provides)
//* 	class TVKalSite
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************
//

#include <iostream>
#include "TVKalSite.h"
#include "TVKalState.h"

//_____________________________________________________________________
//  ------------------------------
//  Base Class for measurement vector used by Kalman filter
//  ------------------------------
//
ClassImp(TVKalSite)

TVKalSystem *TVKalSite::fgKalSysPtr = 0;

TVKalSite::TVKalSite(Int_t m, Int_t p)
                   :TObjArray(2),
                    TAttLockable(),
                    fCurStatePtr(0), 
                    fM(m,1),
                    fV(m,m),
                    fH(m,p),
                    fHt(p,m),
                    fResVec(m,1),
                    fR(m,m),
                    fDeltaChi2(0.)
{
   // Create fStateVector at constractor of concreate class:
   // SetStateVector(new TXXXKalmanStateVector(.....))
}

TVKalSite::~TVKalSite()
{
}

//---------------------------------------------------------------
// Filter
//---------------------------------------------------------------

Bool_t TVKalSite::Filter()
{

   // prea and preC should be preset by TVKalState::Propagate()
//std::cerr << "break1" << std::endl;
   TVKalState &prea = GetState(TVKalSite::kPredicted);
   TKalMatrix h = fM;
   if (!CalcExpectedMeasVec(prea,h)) return kFALSE;
//std::cerr << "break2" << std::endl;
   TKalMatrix pull  = fM - h;
   TKalMatrix preC  = GetState(TVKalSite::kPredicted).GetCovMat();
   // Calculate fH and fHt

   if (!CalcMeasVecDerivative(prea,fH)) return kFALSE;
//std::cerr << "break3" << std::endl;
   fHt = TKalMatrix(TKalMatrix::kTransposed, fH);

   // Calculate covariance matrix of residual

   TKalMatrix preR    = fV + fH * preC * fHt;
   TKalMatrix preRinv = TKalMatrix(TKalMatrix::kInverted, preR);

   // Calculate kalman gain matrix

   TKalMatrix preCinv = TKalMatrix(TKalMatrix::kInverted, preC);
   TKalMatrix G       = TKalMatrix(TKalMatrix::kInverted, fV);

   // Calculate filtered state vector

   TKalMatrix curCinv = preCinv + fHt * G * fH; 
   TKalMatrix curC    = TKalMatrix(TKalMatrix::kInverted, curCinv);
   TKalMatrix K       = curC * fHt * G;

   TKalMatrix Kpull  = K * pull;
   TKalMatrix Kpullt = TKalMatrix(TKalMatrix::kTransposed,Kpull);
   TKalMatrix av     = prea + Kpull;
   TVKalState &a     = CreateState(av,curC,TVKalSite::kFiltered);
   TVKalState *aPtr  = &a;

   Add(aPtr);
   SetOwner();

   // Calculate chi2 increment

   fR      = fV - fH * curC *fHt;
   if (!CalcExpectedMeasVec(a,h)) return kFALSE;
//std::cerr << "break4" << std::endl;
   fResVec = fM - h;
   TKalMatrix curResVect = TKalMatrix(TKalMatrix::kTransposed, fResVec);
   fDeltaChi2 = (curResVect * G * fResVec + Kpullt * preCinv * Kpull)(0,0);

   if (IsAccepted()) return kTRUE;
   else              return kFALSE;
}

//---------------------------------------------------------------
// Smooth
//---------------------------------------------------------------

void TVKalSite::Smooth(TVKalSite &pre)
{

   TVKalState &cura  = GetState(TVKalSite::kFiltered);
   TVKalState &prea  = pre.GetState(TVKalSite::kPredicted);
   TVKalState &sprea = pre.GetState(TVKalSite::kSmoothed);

   TKalMatrix curC    = cura.GetCovMat();
   TKalMatrix curFt   = cura.GetPropMat("T");
   TKalMatrix preC    = prea.GetCovMat();
   TKalMatrix spreC   = sprea.GetCovMat();
   TKalMatrix preCinv = TKalMatrix(TKalMatrix::kInverted, preC);
   TKalMatrix curA    = curC * curFt * preCinv;
   TKalMatrix curAt   = TKalMatrix(TKalMatrix::kTransposed, curA);
   TKalMatrix scurC   = curC + curA * (spreC - preC) * curAt;

   TKalMatrix sv = cura + curA * (sprea - prea);
   Add(&CreateState(sv,scurC,TVKalSite::kSmoothed));
   SetOwner();
}
