//*************************************************************************
//* =====================
//*  TKalTrackSite Class
//* =====================
//*
//* (Description)
//*   Track measurement site class used by Kalman filter.
//* (Requires)
//*     TKalTrackState
//* (Provides)
//*     class TKalTrackSite
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include <iostream>
#include <memory>
//#include "TMath.h"

#include "TKalTrackSite.h"
#include "TKalTrackState.h"
#include "TVTrackHit.h"

#include "TVSurface.h"

using namespace std;

//_____________________________________________________________________
//  ----------------------------------
//  Sample class for measurement site
//  ----------------------------------
//
ClassImp(TKalTrackSite)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------
TKalTrackSite::TKalTrackSite(Int_t m, Int_t p)
              : TVKalSite(m,p), fHitPtr(0), fX0(), fIsHitOwner(kFALSE)
{
}

TKalTrackSite::TKalTrackSite(const TVTrackHit &ht,
                                   Int_t       m, 
                                   Int_t       p)
              : TVKalSite(m,p),
                fHitPtr((TVTrackHit *)&ht), 
                fX0(),
                fIsHitOwner(kFALSE)
{
   for (Int_t i=0; i<m; i++) {
      GetMeasVec     ()(i,0) = ht.GetX(i);
      GetMeasNoiseMat()(i,i) = TMath::Power(ht.GetDX(i),2);
   }
   fX0 = ht.GetMeasLayer().HitToXv(ht); // Set pivot to the hit position
}

TKalTrackSite::~TKalTrackSite()
{
   if (IsHitOwner() && fHitPtr) delete fHitPtr;
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods
//  ----------------------------------

TVKalState & TKalTrackSite::CreateState(const TKalMatrix &sv, Int_t type)
{
   SetOwner();
   return *(new TKalTrackState(sv,*this,type));
}

TVKalState & TKalTrackSite::CreateState(const TKalMatrix &sv,
                                        const TKalMatrix &c,
                                              Int_t       type)
{
   SetOwner();
   return *(new TKalTrackState(sv,c,*this,type));
}

Int_t TKalTrackSite::CalcXexp(const TVKalState &a, 
                                    TVector3   &xx, 
                                    Double_t   &phi) const
{
   std::auto_ptr<TVTrack> hel(&static_cast<const TKalTrackState &>(a).CreateTrack());

   const TVSurface &ms = dynamic_cast<const TVSurface &>(GetHit().GetMeasLayer());
   return ms.CalcXingPointWith(*hel,xx,phi);
}



Int_t TKalTrackSite::CalcExpectedMeasVec(const TVKalState &a, TKalMatrix &h)
{
   Double_t phi = 0.;
   TVector3 xxv;
   if (!CalcXexp(a,xxv,phi)) return 0;	// no hit

   if (a.GetNrows() == 6) h = GetHit().XvToMv(xxv,a(5,0));
   else                   h = GetHit().XvToMv(xxv,0.);

   return 1;
}

Int_t TKalTrackSite::CalcMeasVecDerivative(const TVKalState &a,
                                                 TKalMatrix &H)
{
   // Calculate 
   //    H = (@h/@a) = (@d/@a, @z/@a)^t
   // where
   //        h(a) = (d, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   TVector3      xxv;
   Double_t      phi = 0.;

   if (!CalcXexp(a,xxv,phi)) return 0; // hit on S(x) = 0

   const TVSurface &ms = dynamic_cast<const TVSurface &>(GetHit().GetMeasLayer());
   TKalMatrix    dsdx(ms.CalcDSDx(xxv));   // (@S(x)/@x)

   std::auto_ptr<TVTrack> hel(&static_cast<const TKalTrackState &>(a).CreateTrack());

   TKalMatrix    dxda   = hel->CalcDxDa(phi);  	     // (@x(phi,a)/@a)
   TKalMatrix    dxdphi = hel->CalcDxDphi(phi);	     // (@x(phi,a)/@phi)

   TKalMatrix dphida = dsdx * dxda;
   TKalMatrix dsdphi = dsdx * dxdphi;
   Double_t denom = -dsdphi(0,0);
   dphida *= 1/denom;

   TKalMatrix dxphiada = dxdphi * dphida + dxda; // (@x(phi(a),a)/@a)

   GetHit().GetMeasLayer().CalcDhDa(GetHit(), xxv, dxphiada, H); // H = (@h/@a)

   return 1;
}

Bool_t TKalTrackSite::IsAccepted()
{
   // return kTRUE if this site is accepted by Filter()
#if 1
   Double_t delchi2 = GetDeltaChi2();
   if (delchi2 > 25.) {
      cerr << ">>>> TKalTrackSite::IsAccepted >>>>>>>>>>>>> " << endl
           << " Too big chi2 increment!!! " << endl;
      DebugPrint();
      //return kFALSE;
   }
#endif
   return kTRUE;
}

void TKalTrackSite::DebugPrint() const
{
   cerr << " dchi2 = " << GetDeltaChi2()   << endl;
   cerr << " res_d = " << (*(TKalTrackSite *)this).GetResVec()(0,0) << endl;
   cerr << " res_z = " << (*(TKalTrackSite *)this).GetResVec()(1,0) << endl;
   fHitPtr->DebugPrint();
}

