#include "EXEventGen.h"
#include "EXMeasLayer.h"
#include "EXKalDetector.h"
#include "TTrackFrame.h"
#include "TKalDetCradle.h"
#include "TBField.h"
#include "EXHit.h"
#include "TRandom.h"
#include "TBField.h"
#include <iostream>

// -----------------------------------
//  Track Parameters
// -----------------------------------

#define __DR__     0.
#define __FI0__    0.
#define __DZ__     0.
#define __COSMN__ -0.7
#define __COSMX__ +0.7
#define __X0__     0.
#define __Y0__     0.
#define __Z0__     0.

ClassImp(EXEventGen)

Double_t EXEventGen::fgT0 = 0.; // [nsec]

THelicalTrack EXEventGen::GenerateHelix(Double_t pt)
{
   // ---------------------------
   //  Create a helix track
   // ---------------------------

   Double_t dr  = __DR__;
   Double_t fi0 = __FI0__ + 2*TMath::Pi()*(gRandom->Uniform()-0.5);
   fi0 = 3 + gRandom->Uniform(-0.2, 0.2);
   Double_t cpa = 1. / pt;
   Double_t dz  = __DZ__;
   Double_t cs  = gRandom->Uniform(__COSMN__, __COSMX__);
   Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
   tnl = 0.;
   Double_t x0  = __X0__;
   Double_t y0  = __Y0__;
   Double_t z0  = __Z0__;

   //   
   //make a tilted helix
   //

   TVector3 bfield = TBField::GetGlobalBfield(TVector3());

   TTrackFrame globalFrame;
   TTrackFrame localFrame(globalFrame, TVector3(), bfield);

   THelicalTrack tiltedHelix = THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,bfield.Mag());
   tiltedHelix.SetFrame(localFrame);

   return tiltedHelix;
}

void EXEventGen::Swim(THelicalTrack &heltrk)
{
   // ---------------------------
   //  Swim track and make hits
   // ---------------------------

   Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
                         ->GetSortingPolicy()
                         / heltrk.GetRho();

   const Double_t eps = 1.e-5;
   
   Bool_t   is1stloop = kTRUE;
   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;

   for (Int_t lyr = 0; lyr < nlayers && lyr >= 0; lyr += dlyr) 
   {
	  std::cout << "layer: " << lyr << std::endl;
	  EXMeasLayer &ms = *static_cast<EXMeasLayer *>(fCradlePtr->At(lyr));

	  TVector3 crossingPoint;
      Double_t dfis = dfi;

      if(!ms.CalcXingPointWith(heltrk, crossingPoint, dfi, 1, eps)) 
	  {
		 //no crossing point
		 
         if (is1stloop && fHitBufPtr->GetEntries() > 3) {
            is1stloop = kFALSE;
            dlyr      = -1;
            dfi       = dfis;
            continue;
         } 
		 else
		 { 
			 std::cout << "no crossing point, stop generating" << std::endl;
			 break;
		 }
      }

      // should use the material behind the surface since dfi is measured 
      // from the last point to the current surface
      Bool_t   dir    = dlyr < 0 ? kTRUE : kFALSE;

      const TMaterial &mat = ms.GetMaterial(dir);
      if (fCradlePtr->IsMSOn()) 
	  {
         TKalMatrix Qms(5,5);
         ms.CalcQms(dir, heltrk, dfi, Qms);
         Double_t sgphi  = TMath::Sqrt(Qms(1,1));
         Double_t sgtnl  = TMath::Sqrt(Qms(4,4));
         Double_t delphi = gRandom->Gaus(0.,sgphi);
         Double_t deltnl = gRandom->Gaus(0.,sgtnl);

         heltrk.ScatterBy(delphi,deltnl);  // M.S. at start point
         dfis = 0.;

         if (!ms.CalcXingPointWith(heltrk, crossingPoint, dfi, 1, eps)) break;// recalc exact hit
         dfis += dfi;
      }

	  //Move pivot to current hit. In fact we just use the following 
	  //two matrix to ask for transforming in MoveTo.
      heltrk.MoveTo(crossingPoint, dfi);


	  //FIXME:: rotation first, or MS, dE/dx

	  if(fCradlePtr->IsDEDXOn()) 
	  {
		  TKalMatrix av(5,1);
                  heltrk.PutInto(av); 
		  av(2,0) += ms.GetEnergyLoss(dir, heltrk, dfis); // energy loss
		  heltrk.SetTo(av, heltrk.GetPivot());
	  }

      if (ms.IsActive()) 
	  {
         ms.ProcessHit(crossingPoint, *fHitBufPtr);
      }
   }
}

void EXEventGen::DumpHits()
{
	TObjArray& hitsArray = *fHitBufPtr; 

	for(int i=0; i<hitsArray.GetEntries(); i++)
	{
		EXHit &hit = *dynamic_cast<EXHit *>(hitsArray.At(i));

		TVector3 hitPos = hit.GetMeasLayer().HitToXv(hit);
		TVector3 bfield = TBField::GetGlobalBfield(hitPos);

		std::cout << "Hit " << i << ":" 
			 << "x = "  << hitPos.X()
			 << " y ="  << hitPos.Y()
			 << " z ="  << hitPos.Z()
			 << ", r = " << sqrt(hitPos.X()*hitPos.X()+hitPos.Y()*hitPos.Y())
			 << ", B = (" << bfield.X() << "," << bfield.Y() << "," << bfield.Z() << ")" 
			 << std::endl;
	}
}
