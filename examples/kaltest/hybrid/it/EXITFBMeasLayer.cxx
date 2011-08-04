//*************************************************************************
//* ===================
//*  EXITFBMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXITFBHit.
//* (Requires)
//* (Provides)
//*     class EXITFBMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/07/25  Kim, Youngim      Forward & Backward versoin.
//*
//*   2011/06/30  D.Kamai       Modified to handle turbine-blade-like FTD.
//*************************************************************************
//

#include "EXITFBMeasLayer.h"
#include "EXITFBHit.h"
#include "EXITKalDetector.h"
#include "TVTrack.h"
#include "TRandom.h"
#include <iostream>
#include <iomanip>

#include "TVirtualPad.h"
#include "TTRAP.h"
#include "TNode.h"
#include "TString.h"

using namespace std;
ClassImp(EXITFBMeasLayer)
                                                                                
EXITFBMeasLayer::EXITFBMeasLayer(TMaterial &min,
                                 TMaterial &mout,
                                 const TVector3  &xc,
				 const TVector3  &normal,
				 Double_t   SortingPolicy,
                                 Double_t   rin,
                                 Double_t   rout,
				 Double_t   dxMax,
				 Double_t   dxMin,
                                 Double_t   sigmax,
                                 Double_t   sigmay,
                                 Bool_t     type,
				 Int_t      mode,
				 const Char_t    *name)
               : EXVMeasLayer(min, mout, type, name),
                 TPlane(xc,normal),
		 fSortingPolicy(SortingPolicy),
                 fRin(rin),
                 fRout(rout),
		 fdxMax(dxMax),
		 fdxMin(dxMin),
                 fSigmaX(sigmax),
                 fSigmaY(sigmay),
		 fMode(mode)
{
}

EXITFBMeasLayer::~EXITFBMeasLayer()
{
}

TKalMatrix EXITFBMeasLayer::XvToMv(const TVector3 &xv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = x
   //     (1,0) = y
  
   TKalMatrix mv(kMdim,1);
   mv(0,0)  = Cosalpha()*(Cosphi()*xv.X() - Sinphi()*xv.Y()) + Sinalpha()*(xv.Z() - GetXc().Z());
   mv(1,0)  = Sinphi()*xv.X() + Cosphi()*xv.Y();
   return mv;
}

TKalMatrix EXITFBMeasLayer::XvToMv(const TVTrackHit &,
                                   const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXITFBMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXITFBHit &ht = dynamic_cast<const EXITFBHit &>(vht);
   
   Double_t x = ht(0,0)*Cosphi()*Cosalpha() + ht(1,0)*Sinphi();
   Double_t y = -ht(0,0)*Sinphi()*Cosalpha() + ht(1,0)*Cosphi();
   
   Double_t z = GetXc().Z() + ht(0,0)*Sinalpha();
      
   return TVector3(x, y, z);
}

void EXITFBMeasLayer::CalcDhDa(const TVTrackHit &vht,
                               const TVector3   &xxv,
                               const TKalMatrix &dxphiada,
                                     TKalMatrix &H) const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = 5;

   // Set H = (@h/@a) = (@phi/@a, @r/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
     H(0,i) = Cosalpha()*(Cosphi()*dxphiada(0,i) - Sinphi()*dxphiada(1,i)) + Sinalpha()*dxphiada(2,i);
     H(1,i) = Sinphi()*dxphiada(0,i) + Cosphi()*dxphiada(1,i);
        }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}

Bool_t EXITFBMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  TKalMatrix mv = XvToMv(xx);

  if( TMath::Abs((xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() + (xx.Z()-GetXc().Z())*GetNormal().Z()) < 1e-4){
    if( mv(1,0) <= GetRout() && mv(1,0) >= GetRin() && mv(1,0) >= TMath::Abs(2*GetRout()*mv(0,0)/GetdxMax())){
      if(GetMode()==0 || GetMode()*mv(0,0) >= 0) return kTRUE;
    }
    else{
      return kFALSE;
    }
  }
  return kFALSE;
}

void EXITFBMeasLayer::ProcessHit(const TVector3  &xx,
                                       TObjArray &hits)
{
   TKalMatrix h   = XvToMv(xx);
   Double_t   x = h(0, 0); 
   Double_t   y = h(1, 0); 

   Double_t dx = GetSigmaX();
   Double_t dy = GetSigmaY();
   x += gRandom->Gaus(0., dx);   // smearing x
   y += gRandom->Gaus(0., dy);   // smearing y

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = x;
   meas [1] = y;
   dmeas[0] = dx;
   dmeas[1] = dy;
 
   Double_t b = EXITKalDetector::GetBfield();
   hits.Add(new EXITFBHit(*this, meas, dmeas, xx, b));
}

// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXITFBMeasLayer::Draw(Int_t color, const Char_t *opt)
{
  Double_t upperbase = 0;
  Double_t lowerbase = 0;
  Double_t theta = 0;
  Double_t dy = (GetRout()-GetRin())/2;
  Double_t thick = 0.02; // [cm]
  Double_t z = GetXc().Z();
  {z >0 ? z += thick/2 : z -= thick/2; }
  if (!gPad) return;
  if (!IsActive()) return;
  if (!GetNodePtr()) {
    const Char_t *name  = GetMLName().Data();
    const Char_t *nname = (GetMLName() + "Node").Data();

    if(GetMode()==0){
      upperbase = GetdxMax()/2;
      lowerbase = GetdxMin()/2;
    }else{
      upperbase = GetdxMax()/4;
      lowerbase = GetdxMin()/4;
      theta = GetMode()*TMath::ATan((upperbase-lowerbase)/(2*dy))*180/TMath::Pi();
    }
    TTRAP *trap = new TTRAP(name,name,"Si",dy,theta, 0,
			    thick/2, lowerbase, lowerbase, 0,
			    thick/2, upperbase, upperbase, 0 );      
    trap->SetBit(kCanDelete);
    Double_t rmat[9]   = { Cosphi()*Cosalpha(), -Sinphi()*Cosalpha(), Sinalpha(),
                	   -Cosphi()*Sinalpha(), Sinphi()*Sinalpha(),  Cosalpha(),
			   Sinphi(),           Cosphi(),             0};
    TRotMatrix *rmatp   = new TRotMatrix( "rmat", "rmat", rmat);
    TNode *nodep = new TNode(nname,nname,trap,GetXc().X()+GetMode()*(GetdxMax()+GetdxMin())/8,GetXc().Y(),z+Sinalpha()*GetMode()*(GetdxMax()+GetdxMin())/8,rmatp,"");
    nodep->SetLineColor(color);
    nodep->SetLineWidth(0.01);
    SetNodePtr(nodep);
  }
}
