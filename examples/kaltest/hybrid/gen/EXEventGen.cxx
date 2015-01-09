#include "EXEventGen.h"
#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TPlane.h"
#include "TRandom.h"

//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __FI0__    0.
#define __DZ__     0.
#define __X0__     0.
#define __Y0__     0.
#define __Z0__     0.

ClassImp(EXEventGen)

Double_t EXEventGen::fgT0 = 14.; // [nsec]

THelicalTrack EXEventGen::GenerateHelix(Double_t pt,
                                        Double_t cosmin,
                                        Double_t cosmax)
{
   // ---------------------------
   //  Generate a helical track
   // ---------------------------

   Double_t dr  = __DR__;
   Double_t fi0 = __FI0__ + 2*TMath::Pi()*(gRandom->Uniform()-0.5);
   Double_t cpa = 1. / pt;
   Double_t dz  = __DZ__;
   Double_t cs  = gRandom->Uniform(cosmin, cosmax);
   Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
   Double_t x0  = __X0__;
   Double_t y0  = __Y0__;
   Double_t z0  = __Z0__;

   Double_t b   = dynamic_cast<const EXVKalDetector &>
                 (dynamic_cast<EXVMeasLayer *>
                 (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

   return THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);
}

void EXEventGen::Swim(THelicalTrack &heltrk)
{
   // ---------------------------
   //  Swim track and Make hits
   // ---------------------------

   Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
                            ->GetSortingPolicy()
                         / heltrk.GetRho();

   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;
   Double_t dfisum    = 0.;

   for (Int_t lyr = 0; lyr >= 0; lyr += dlyr) { // loop over layers
      // change direction if it starts looping back
      if (lyr == nlayers - 1) dlyr = -1;

      EXVMeasLayer &ml = *dynamic_cast<EXVMeasLayer *>(fCradlePtr->At(lyr));
      TVSurface    &ms = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyr));
      TVector3 xx;
      Double_t dfis = dfi;
      if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
       || TMath::Abs(dfi) > TMath::Pi()
       || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
         dfi = dfis;
         continue;
      }
      // should use the material behind the surface since dfi is measured 
      // from the last point to the current surface
      Bool_t   dir    = dlyr < 0 ? kTRUE : kFALSE;

      if (fCradlePtr->IsMSOn()) {
         TKalMatrix Qms(5,5);
         ml.CalcQms(dir, heltrk, dfi, Qms);
         Double_t sgphi  = TMath::Sqrt(Qms(1,1));
         Double_t sgtnl  = TMath::Sqrt(Qms(4,4));
         Double_t delphi = gRandom->Gaus(0.,sgphi);
         Double_t deltnl = gRandom->Gaus(0.,sgtnl);
#if 0
         dfi *= 0.5;
         TVector3 x0ms = heltrk.CalcXAt(dfi);
         heltrk.MoveTo(x0ms,dfi);     // M.S. at mid point

         heltrk.ScatterBy(delphi,deltnl);
         dfis = dfi;
#else
         heltrk.ScatterBy(delphi,deltnl); // multiple scattering
         dfis = 0.;
#endif
         // recalculate crossing point
         if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
          || TMath::Abs(dfi) > TMath::Pi()
          || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
            dfi = dfis;
            continue;
         }
         dfis += dfi;
      }
      dfisum += dfi;

      heltrk.MoveTo(xx,dfi);	// move pivot to current hit

      if (fCradlePtr->IsDEDXOn()) {
         TKalMatrix av(5,1);
         heltrk.PutInto(av);
         av(2,0) += ml.GetEnergyLoss(dir, heltrk, dfis); // energy loss
         heltrk.SetTo(av, heltrk.GetPivot());
      }
      if (ml.IsActive() && dynamic_cast<const EXVKalDetector &>(ml.GetParent(kFALSE)).IsPowerOn()) {
         ml.ProcessHit(xx, *fHitBufPtr); // create hit point
      }
      if (lyr == nlayers - 1) break;
   }
}
