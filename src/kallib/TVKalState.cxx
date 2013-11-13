//*************************************************************************
//* ====================
//*  TVKalState Class
//* ====================
//*
//* (Description)
//*   This is the base class for state vector used by Kalman filter.
//* (Requires)
//* 	TKalMatrix
//* (Provides)
//* 	class TVKalState
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version
//*
//*************************************************************************
//
#include "TVKalState.h"
#include "TVKalSite.h"
//_____________________________________________________________________
//  ------------------------------
//  Base Class for measurement vector used by Kalman filter
//  ------------------------------
//
ClassImp(TVKalState)

#ifdef __STATE_TIMER__   
Double_t TVKalState::fTime = 0.;
#endif

TVKalState::TVKalState(Int_t type, Int_t p)
                  : TKalMatrix(p,1),
                    fType(type),
                    fSitePtr(0),
                    fF(p,p),
                    fFt(p,p),
                    fQ(p,p),
                    fC(p,p)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
}

TVKalState::TVKalState(const TKalMatrix &sv, 
                       Int_t type, Int_t p)
                  : TKalMatrix(sv),
                    fType(type),
                    fSitePtr(0),
                    fF(p,p),
                    fFt(p,p),
                    fQ(p,p),
                    fC(p,p)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
}

TVKalState::TVKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       Int_t type, Int_t p)
                  : TKalMatrix(sv),
                    fType(type),
                    fSitePtr(0),
                    fF(p,p),
                    fFt(p,p),
                    fQ(p,p),
                    fC(c)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
}

TVKalState::TVKalState(const TKalMatrix &sv, const TVKalSite &site, 
                       Int_t type, Int_t p)
                  : TKalMatrix(sv),
                    fType(type),
                    fSitePtr((TVKalSite *)&site),
                    fF(p,p),
                    fFt(p,p),
                    fQ(p,p),
                    fC(p,p)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
}

TVKalState::TVKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       const TVKalSite &site, Int_t type, Int_t p)
                  : TKalMatrix(sv),
                    fType(type),
                    fSitePtr((TVKalSite *)&site),
                    fF(p,p),
                    fFt(p,p),
                    fQ(p,p),
                    fC(c)
{
   fF.UnitMatrix();
   fFt.UnitMatrix();
}

//-------------------------------------------------------------------
// Propagate
//------------------------------------------------------------------

void TVKalState::Propagate(TVKalSite &to)
{
#ifdef __STATE_TIMER__   
   clock_t t1, t2;
   t1 = clock();
#endif

   // Calculate 
   //    prea:	predicted state vector      : a^k-1_k = f_k-1(a_k-1)
   //    fF:    propagator derivative       : F_k-1   = (@f_k-1/@a_k-1)
   //    fQ:    process noise from k-1 to k : Q_k-1)

   TVKalState &prea    = MoveTo(to,fF,fQ);
   TVKalState *preaPtr = &prea;

   fFt = TKalMatrix(TKalMatrix::kTransposed, fF);

   // Calculate covariance matrix

   TKalMatrix preC = fF * fC * fFt + fQ;

   // Set predicted state vector and covariance matrix to next site

   prea.SetCovMat(preC);
   to.Add(preaPtr);
   to.SetOwner();

#ifdef __STATE_TIMER__   
   t2 = clock();
   fTime += Double_t(t2-t1)/CLOCKS_PER_SEC;
#endif
}
