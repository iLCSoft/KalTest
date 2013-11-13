//*************************************************************************
//* ===================
//*  EXBPConeConeMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXBPConeHit.
//* (Requires)
//* (Provides)
//*     class EXBPConeMeasLayer
//* (Update Recored)
//*   2012/01/19  K.Fujii       Original version.
//*
//*************************************************************************
//

#include "EXBPConeMeasLayer.h"
#include "EXBPConeHit.h"
#include "EXBPKalDetector.h"
#include "TRandom.h"
#include "TMath.h"

#include "TVirtualPad.h"
#include "TCONE.h"
#include "TNode.h"
#include "TString.h"

#if 1
#include <iostream>
#endif

ClassImp(EXBPConeMeasLayer)
                                                                                
EXBPConeMeasLayer::EXBPConeMeasLayer(TMaterial &min,
                                     TMaterial &mout,
                                     Double_t   z1,
                                     Double_t   r1,
                                     Double_t   z2,
                                     Double_t   r2,
                                     Double_t   sigmax,
                                     Double_t   sigmaz,
                                     Bool_t     type,
                               const Char_t    *name)
                 : EXVMeasLayer(min, mout, type, name),
                   TCutCone(r1*(z2-z1)/(r2-r1), 
                            r2*(z2-z1)/(r2-r1), 
                               (r2-r1)/(z2-z1),
                            0.,0.,(r2*z1-r1*z2)/(r2-r1)),
                   fZ1(z1),
                   fR1(r1),
                   fZ2(z2),
                   fR2(r2),
                   fSigmaX(sigmax),
                   fSigmaZ(sigmaz)
{
}

EXBPConeMeasLayer::~EXBPConeMeasLayer()
{
}

TKalMatrix EXBPConeMeasLayer::XvToMv(const TVector3 &xxv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = r * phi 
   //     (1,0) = z

   TKalMatrix mv(kMdim,1);
   TVector3 xv = xxv - GetXc();
   Double_t r  = xv.Z()*GetTanA();
    
   mv(0,0)  = r * TMath::ATan2(xv.Y(), xv.X());
   mv(1,0)  = xv.Z();
   return mv;
}

TKalMatrix EXBPConeMeasLayer::XvToMv(const TVTrackHit &,
                                     const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXBPConeMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXBPConeHit &ht = dynamic_cast<const EXBPConeHit &>(vht);

   Double_t r   = ht(1,0) * GetTanA();
   Double_t phi = ht(0,0) / r;
   Double_t x   = GetXc().X() + r * TMath::Cos(phi);
   Double_t y   = GetXc().Y() + r * TMath::Sin(phi);
   Double_t z   = GetXc().Z() + ht(1,0);

   return TVector3(x,y,z);
}

void EXBPConeMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                 const TVector3   &xxv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = TMath::Max(5,sdim-1);

   TVector3 xv = xxv - GetXc();
   Double_t x  = xv.X();
   Double_t y  = xv.Y();
   Double_t z  = xv.Z();
   Double_t xxyy = x * x + y * y;
   Double_t phi  = TMath::ATan2(y, x);
   Double_t tana = GetTanA();
   Double_t r    = z * GetTanA();

   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = - r * (y / xxyy) * dxphiada(0,i) 
               + r * (x / xxyy) * dxphiada(1,i)
               +     tana * phi * dxphiada(2,i);
      H(1,i) =  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}

Bool_t EXBPConeMeasLayer::IsOnSurface(const TVector3 &xx) const
{
    TVector3 xxc = xx - GetXc();
    Double_t r   = xxc.Perp();
    Double_t z   = xxc.Z();
    Double_t s   = (r - GetTanA()*z) * (r + GetTanA()*z);
    const Double_t kTol = 1.e-8;

#if 0
    std::cerr << GetMLName() << ":" << std::endl;
    std::cerr << "s=" << s << " xx=(" << xx.X() << "," << xx.Y() << "," << xx.Z() << ")" << std::endl;
    std::cerr << "bool=" << (TMath::Abs(s) < kTol && ((xx.Z()-fZ1)*(xx.Z()-fZ2) <= 0.)) << std::endl;
    std::cerr << "fZ1=" << fZ1 << " fZ2=" << fZ2 << std::endl;
#endif

    return (TMath::Abs(s) < kTol && ((xx.Z()-fZ1)*(xx.Z()-fZ2) <= 0.));
} 


void EXBPConeMeasLayer::ProcessHit(const TVector3  &xx,
                                         TObjArray &hits)
{
   TKalMatrix h    = XvToMv(xx);
   Double_t   rphi = h(0, 0);
   Double_t   z    = h(1, 0);

   Double_t dx = GetSigmaX();
   Double_t dz = GetSigmaZ();
   rphi += gRandom->Gaus(0., dx);   // smearing rphi
   z    += gRandom->Gaus(0., dz);   // smearing z

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = rphi;
   meas [1] = z;
   dmeas[0] = dx;
   dmeas[1] = dz;

   Double_t b = EXBPKalDetector::GetBfield(xx);
   hits.Add(new EXBPConeHit(*this, meas, dmeas, xx, b));
}

// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXBPConeMeasLayer::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   if (!IsActive()) return;
   if (!GetNodePtr()) {
      const Char_t *name  = GetMLName().Data();
      const Char_t *nname = (GetMLName() + "Node").Data();
      Double_t hlen = TMath::Abs(fZ2-fZ1)/2;
      Double_t r1   = fZ2 > fZ1 ? fR1 : fR2;
      Double_t r2   = fZ2 > fZ1 ? fR2 : fR1;
      TCONE *conep = new TCONE(name, name, "void", hlen, r1, r1, r2, r2);
      conep->SetBit(kCanDelete);
      TNode *nodep = new TNode(nname,nname,name,0.,0.,(fZ1+fZ2)/2);
      nodep->SetLineColor(color);
      nodep->SetLineWidth(0.01);
      SetNodePtr(nodep);
   }
}
