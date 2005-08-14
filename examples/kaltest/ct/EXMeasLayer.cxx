//*************************************************************************
//* ===================
//*  EXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXHit.
//* (Requires)
//* (Provides)
//*     class EXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXMeasLayer.h"
#include "EXHit.h"

Bool_t   EXMeasLayer::kActive = kTRUE;
Bool_t   EXMeasLayer::kDummy = kFALSE;
//Double_t EXMeasLayer::fgSigmaX = 4.e-4;
//Double_t EXMeasLayer::fgSigmaZ = 4.e-4;
Double_t EXMeasLayer::fgSigmaX = 1.e-2;
Double_t EXMeasLayer::fgSigmaZ = 1.e-2;

ClassImp(EXMeasLayer)
                                                                                
EXMeasLayer::EXMeasLayer(TMaterial &min,      // material inside the layer
                         TMaterial &mout,     // material outside the layer
                         Double_t   r0,       // radius of the cylindrial layer
                         Double_t   lhalf,    // half length of the cylindrical layer
                         Bool_t     isactive) // flag to tell the layer is active
           : TVMeasLayer(min, mout, isactive),
             TCylinder(r0, lhalf)
{
}

EXMeasLayer::~EXMeasLayer()
{
}

TKalMatrix EXMeasLayer::XvToMv(const TVector3 &xv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = r * phi 
   //     (1,0) = z

   TKalMatrix mv(kMdim,1);
   mv(0,0)  = GetR() * TMath::ATan2(xv.Y(), xv.X());
   mv(1,0)  = xv.Z();
   return mv;
}

TKalMatrix EXMeasLayer::XvToMv(const TVTrackHit &,
                               const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXHit &ht = dynamic_cast<const EXHit &>(vht);

   Double_t phi = ht(0,0) / GetR();
   Double_t z   = ht(1,0);
   Double_t x   = GetR() * TMath::Cos(phi);
   Double_t y   = GetR() * TMath::Sin(phi);

   return TVector3(x,y,z);
}

void EXMeasLayer::CalcDhDa(const TVTrackHit &,          // Hit: not used in this example
                           const TVector3   &xxv,       // hit position vector
                           const TKalMatrix &dxphiada,  // @x(\phi(a),a)/@a
                                 TKalMatrix &H)  const  // projector matrix = @h/@a
{
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
      H(1,i) =  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}
