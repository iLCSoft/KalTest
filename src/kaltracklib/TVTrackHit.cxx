//*************************************************************************
//* ====================
//*  TVTrackHit Class
//* ====================
//*
//* (Description)
//*   Abstract base class to store single hit information.
//* (Requires)
//* (Provides)
//*     class TVTrackHit
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************

#include "TVTrackHit.h"

ClassImp(TVTrackHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

TVTrackHit::TVTrackHit(Int_t m)
        : TKalMatrix(m,2), fDim(m), fXX(), fBfield(30.), fMeasLayerPtr(0)
{
}
      
TVTrackHit::TVTrackHit(const TVMeasLayer &ms,
                             Double_t    *x,
                             Double_t    *dx, 
                       const TVector3    &xx,
                             Double_t     b,
                             Int_t        m)
        : TKalMatrix(m,2), fDim(m), fXX(xx), fBfield(b),
          fMeasLayerPtr((TVMeasLayer *)&ms)
{
   for (Int_t i=0; i<m; i++) {
      (*this)(i,0) = x [i];
      (*this)(i,1) = dx[i];
   }
}

TVTrackHit::TVTrackHit(const TVMeasLayer &ms,
                             Double_t    *x,
                             Double_t    *dx, 
                             Double_t     b,
                             Int_t        m)
        : TKalMatrix(m,2), fDim(m), fXX(0.), fBfield(b),
          fMeasLayerPtr((TVMeasLayer *)&ms)
{
   for (Int_t i=0; i<m; i++) {
      (*this)(i,0) = x [i];
      (*this)(i,1) = dx[i];
   }
}

TVTrackHit::TVTrackHit(const TVTrackHit &hit)
          : TKalMatrix(hit), 
            fDim(hit.fDim),
            fXX(hit.fXX),
            fBfield(hit.fBfield),
            fMeasLayerPtr(hit.fMeasLayerPtr)
{
}

TVTrackHit::~TVTrackHit()
{
}
