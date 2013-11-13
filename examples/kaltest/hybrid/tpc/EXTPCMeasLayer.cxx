//*************************************************************************
//* ===================
//*  EXTPCMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXTPCHit.
//* (Requires)
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "EXTPCKalDetector.h"
#include "EXEventGen.h"
#include "TRandom.h"

ClassImp(EXTPCMeasLayer)
                                                                                
EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
                               TMaterial &mout,
                               Double_t   r0,
                               Double_t   lhalf,
                               Double_t   sigmax0,
                               Double_t   sigmax1,
                               Double_t   sigmaz,
                               Bool_t     type,
                         const Char_t     *name)
              : EXVMeasLayer(min, mout, type, name),
                TCylinder(r0, lhalf),
                fSigmaX0(sigmax0),
                fSigmaX1(sigmax1),
                fSigmaZ(sigmaz)
{
}

EXTPCMeasLayer::~EXTPCMeasLayer()
{
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVector3 &xv,
                                        Int_t     side) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = r * phi 
   //     (1,0) = drift distance

   TKalMatrix mv(kMdim,1);
   mv(0,0)  = GetR() * TMath::ATan2(xv.Y(), xv.X());
   mv(1,0)  = (GetLength() * 0.5) - side * xv.Z();
   return mv;
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVTrackHit &vht,
                                  const TVector3   &xv) const
{
   return XvToMv(xv, dynamic_cast<const EXTPCHit &>(vht).GetSide());
}

TVector3 EXTPCMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);

   Double_t phi = ht(0,0) / GetR();
   Double_t z   = ht.GetSide() * (GetLength() * 0.5 - ht(1,0));
   Double_t x   = GetR() * TMath::Cos(phi);
   Double_t y   = GetR() * TMath::Sin(phi);

   return TVector3(x,y,z);
}

void EXTPCMeasLayer::CalcDhDa(const TVTrackHit &vht,
                              const TVector3   &xxv,
                              const TKalMatrix &dxphiada,
                                    TKalMatrix &H)  const
{
   const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);

   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = TMath::Max(5,sdim-1);

   Double_t xv = xxv.X();
   Double_t yv = xxv.Y();
   Double_t xxyy = xv * xv + yv * yv;

   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = - (yv / xxyy) * dxphiada(0,i) 
               + (xv / xxyy) * dxphiada(1,i);
      H(0,i) *= GetR();
      H(1,i) = - ht.GetSide() *  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = ht.GetVdrift();
   }
}

Double_t EXTPCMeasLayer::GetSigmaX(Double_t zdrift) const
{
   return TMath::Sqrt(fSigmaX0 * fSigmaX0 + fSigmaX1 * fSigmaX1 * zdrift);
}

void EXTPCMeasLayer::ProcessHit(const TVector3  &xx,
                                      TObjArray &hits)
{
   Int_t      side = (xx.Z() < 0. ? -1 : 1);
   TKalMatrix h    = XvToMv(xx, side);
   Double_t   rphi = h(0, 0);
   Double_t   d    = h(1, 0);

   Double_t dx = GetSigmaX(d);
   Double_t dz = GetSigmaZ();
   rphi += gRandom->Gaus(0., dx);   // smearing rphi
   d    += gRandom->Gaus(0., dz);   // smearing drift distance

   Double_t v = EXTPCKalDetector::GetVdrift();
   d         += v * EXEventGen::GetT0(); // T0 shift

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = rphi;
   meas [1] = d;
   dmeas[0] = dx;
   dmeas[1] = dz;

   Double_t b = EXTPCKalDetector::GetBfield(xx);
   hits.Add(new EXTPCHit(*this, meas, dmeas, side, v, xx, b));
}

