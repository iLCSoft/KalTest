//*************************************************************************
//* ===================
//*  EXVTXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXVTXHit.
//* (Requires)
//* (Provides)
//*     class EXVTXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*   2011/06/17  D.Kamai           Modified to handle ladder structure.
//*************************************************************************
//
#include <iostream>
#include "EXVTXMeasLayer.h"
#include "EXVTXHit.h"
#include "EXVTXKalDetector.h"
#include "TVTrack.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

ClassImp(EXVTXMeasLayer)
                                                                                
EXVTXMeasLayer::EXVTXMeasLayer(TMaterial &min,
                               TMaterial &mout,
			       const TVector3  &center,
			       const TVector3  &normal,
			       Double_t SortingPolicy,
			       Double_t   xiwidth,
			       Double_t   zetawidth,
			       Double_t   xioffset,
			       Double_t   sigmaxi,
                               Double_t   sigmazeta,
                               Bool_t     type,
                         const Char_t    *name)
              : EXVMeasLayer(min, mout, type, name),
                TPlane(center, normal),
		fSortingPolicy(SortingPolicy),
		fXiwidth(xiwidth),
		fZetawidth(zetawidth),
		fXioffset(xioffset),
		fSigmaXi(sigmaxi),
                fSigmaZeta(sigmazeta)
{
}

EXVTXMeasLayer::~EXVTXMeasLayer()
{
}

TKalMatrix EXVTXMeasLayer::XvToMv(const TVector3 &xv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = xi 
   //     (1,0) = zeta

   TKalMatrix mv(kMdim,1);

   mv(0,0) = (xv.X() - GetXc().X())*GetNormal().Y()/GetNormal().Perp() - (xv.Y() - GetXc().Y())*GetNormal().X()/GetNormal().Perp();
   mv(1,0) = xv.Z();

   return mv;
}

TKalMatrix EXVTXMeasLayer::XvToMv(const TVTrackHit &,
                               const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXVTXMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXVTXHit &ht = dynamic_cast<const EXVTXHit &>(vht);

   Double_t z = ht(1,0);
   Double_t x = ht(0,0)*GetNormal().Y()/GetNormal().Perp() + GetXc().X();
   Double_t y =-ht(0,0)*GetNormal().X()/GetNormal().Perp() + GetXc().Y();
   
   return TVector3(x,y,z);
}

void EXVTXMeasLayer::CalcDhDa(const TVTrackHit &vht,
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

   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
     H(0,i) =  (GetNormal().Y() / GetNormal().Perp()) * dxphiada(0,i)
              -(GetNormal().X() / GetNormal().Perp()) * dxphiada(1,i);
     H(1,i) =  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }

}

void EXVTXMeasLayer::ProcessHit(const TVector3  &xx,
                                      TObjArray &hits)
{
   TKalMatrix h    = XvToMv(xx);

   Double_t xi = h(0,0);
   Double_t zeta = h(1,0);
   Double_t dxi = GetSigmaXi();
   Double_t dzeta = GetSigmaZeta();
   xi   += gRandom->Gaus(0., dxi);   // smearing rphi
   zeta += gRandom->Gaus(0., dzeta);   // smearing z

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = xi;
   meas [1] = zeta;
   dmeas[0] = dxi;
   dmeas[1] = dzeta;

   Double_t b = EXVTXKalDetector::GetBfield();
   hits.Add(new EXVTXHit(*this, meas, dmeas, xx, b));
}


Bool_t EXVTXMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  Double_t xi   = (xx.X()-GetXc().X())*GetNormal().Y()/GetNormal().Perp() - (xx.Y()-GetXc().Y())*GetNormal().X()/GetNormal().Perp();
  Double_t zeta = xx.Z();

  if( (xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() < 1e-4){
     if( xi <= GetXiwidth()/2 - GetXioffset() && xi >= - GetXiwidth()/2 -  GetXioffset() && TMath::Abs(zeta) <= GetZetawidth()/2){
      return kTRUE;
    }
    else{
      return kFALSE;
    }
  }
}


// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXVTXMeasLayer::Draw(Int_t color, const Char_t *opt)
{
  Double_t phi0 =  1.570796327;
  Double_t phi;
   if (!gPad) return;
   //if (!IsActive()) return;
   if (!GetNodePtr()){
     const Char_t *name  = GetMLName().Data();
     const Char_t *nname = (GetMLName() + "Node").Data();
     TBRIK *brikp        = new TBRIK(name,name, "Si", GetXiwidth()/2, 0, GetZetawidth()/2);
     brikp->SetBit(kCanDelete);
     phi = GetNormal().Phi();
     Double_t dphi = (Double_t)(phi - phi0);
     Double_t rmat[9]   = {  TMath::Cos(dphi), TMath::Sin(dphi), 0,
			    -TMath::Sin(dphi), TMath::Cos(dphi), 0,
			                   0 ,               0 , 1};
     TRotMatrix *rmatp   = new TRotMatrix( "rmat", "rmat", rmat);
     TNode *nodep        = new TNode(nname,nname,brikp, GetXc().X() - GetXioffset()*TMath::Sin(phi), GetXc().Y() + GetXioffset()*TMath::Cos(phi), 0, rmatp,"");
     nodep->SetLineColor(color);
     nodep->SetLineWidth(0.01);	 
     SetNodePtr(nodep);    

   }

}
