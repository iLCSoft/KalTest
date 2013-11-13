#include "EXEventGen.h"
#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TTrackFrame.h"
#include "TPlane.h"
#include "TRandom.h"
#include <iostream>

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

Double_t EXEventGen::fgT0 = 0.; // [nsec]

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

   TVector3 bfield = TBField::GetGlobalBfield(TVector3());
   TTrackFrame globalFrame;
   TTrackFrame localFrame(globalFrame, TVector3(), bfield);
   
   THelicalTrack tiltedHelix = THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,bfield.Mag());
   tiltedHelix.SetFrame(localFrame);

   return tiltedHelix;
}

void EXEventGen::Swim(THelicalTrack &heltrk, Double_t mass)
{
   // ---------------------------
   //  Swim track and Make hits
   // ---------------------------

   const Double_t eps = 1.e-5;
   Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
                            ->GetSortingPolicy()
                         / heltrk.GetRho();

   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;
   Double_t dfisum    = 0.;

   for (Int_t lyr = 0; lyr >= 0 && lyr < nlayers; lyr += dlyr) // loop over layers
   {
      EXVMeasLayer &ml = *dynamic_cast<EXVMeasLayer *>(fCradlePtr->At(lyr));
      TVSurface    &ms = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyr));

      TVector3 crossingPoint;
	  Double_t dfis = dfi;

      if (!ms.CalcXingPointWith(heltrk,crossingPoint,dfi,1,eps))
      {
         dfi = 5*dfis;
         Int_t lyrlast = lyr - dlyr;
         if (lyrlast < 0 || lyrlast >= nlayers - 1) {
            dfi = dfis;
            continue;
         }
         TVSurface    &msl = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyrlast));
         
         if (!msl.CalcXingPointWith(heltrk,crossingPoint,dfi,1,eps))
         {
            dfi = dfis;
         } else {
            if (TMath::Abs(dfi) < TMath::Pi()/4 && TMath::Abs(dfi) > TMath::Abs(dfis)) {
               // start looping back here
               dlyr = -dlyr;
            } else {
               dfi  = dfis;
            }
         }
         continue;
      } else if (TMath::Abs(dfi) > TMath::Pi() ||
                 TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
         dfi = dfis;
         continue;
      }
#if 0
      else {
         std::cerr << " -- hit found: lyr=" << lyr
                   << " r=" << crossingPoint.Perp() << std::endl;
         std::cerr << " xx=(" << crossingPoint.X() << ","
                              << crossingPoint.Y() << ","
                              << crossingPoint.Z() << ")"
                              << std::endl;
      }
#endif
      // should use the material behind the surface since dfi is measured
      // from the last point to the current surface
      Bool_t   dir    = dlyr < 0 ? kTRUE : kFALSE;

      const TMaterial &mat = ml.GetMaterial(dir);
      if (fCradlePtr->IsMSOn()) {
         TKalMatrix Qms(5,5);
         fCradlePtr->CalcQms(mat, heltrk, dfi, Qms, mass);
         Double_t sgphi  = TMath::Sqrt(Qms(1,1));
         Double_t sgtnl  = TMath::Sqrt(Qms(4,4));
         Double_t delphi = gRandom->Gaus(0.,sgphi);
         Double_t deltnl = gRandom->Gaus(0.,sgtnl);

         heltrk.ScatterBy(delphi,deltnl); // multiple scattering
         dfis = 0.;

         // recalculate crossing point
         if (!ms.CalcXingPointWith(heltrk,crossingPoint,dfi,1,eps) || 
			 TMath::Abs(dfi) > TMath::Pi() || 
			 TMath::Abs(dfi + dfisum) > TMath::TwoPi())
		 {
            dfi = dfis;
            continue;
         }

         dfis += dfi;
      }

      dfisum += dfi;

	  //move pivot to current global hit
	  heltrk.MoveTo(crossingPoint, dfi);

      if (fCradlePtr->IsDEDXOn()) 
	  {
         TKalMatrix av(5,1);
         heltrk.PutInto(av);
         av(2,0) += fCradlePtr->GetEnergyLoss(mat, heltrk, dfis, mass); // energy loss
         heltrk.SetTo(av, heltrk.GetPivot());
      }

      if (ml.IsActive() && dynamic_cast<const EXVKalDetector &>(ml.GetParent(kFALSE)).IsPowerOn()) 
	  {
         ml.ProcessHit(crossingPoint, *fHitBufPtr); // create hit point
#if 0
          std::cerr << "lyr: " << lyr << " r=" << crossingPoint.Perp()
          << " xvhit=(" << crossingPoint.X() << ", "
          << crossingPoint.Y() << ", " << crossingPoint.Z() << ")" << std::endl;
#endif
      }
      if (lyr == nlayers - 1) break;
   }
}
