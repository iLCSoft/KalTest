#include "EXRungeKuttaEventGen.h"
#include "EXMeasLayer.h"
#include "EXKalDetector.h"
#include "TKalDetCradle.h"
#include "TBField.h"
#include "EXHit.h"
#include "TRandom.h"
#include "TBField.h"
#include <iostream>

ClassImp(EXRungeKuttaEventGen)

Double_t EXRungeKuttaEventGen::fgT0 = 0.; // [nsec]

TRungeKuttaTrack EXRungeKuttaEventGen::GenerateRKTrack(Double_t chg, TVector3 x, TVector3 p)
{
   // Create a Runge-Kutta track
   TRungeKuttaTrack rkTrack(chg, x, p);

   return rkTrack;
}

void EXRungeKuttaEventGen::Swim(TRungeKuttaTrack &rktrk, Double_t mass)
{
   // ---------------------------
   //  Swim track and make hits
   // ---------------------------
   //Bool_t   is1stloop = kTRUE;
   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;

   const Double_t STEP = 0.1; //mm

   TVector3 currentPosition; 
   TVector3 lastPosition = currentPosition;

   TVector3 currentMomentum; 
   TVector3 lastMomentum = currentMomentum;

   EXMeasLayer* lastml = 0;
//#define _DEBUG_
   for (Int_t lyr = 0; lyr < nlayers && lyr >= 0; lyr += dlyr)
   {
      EXMeasLayer& ml = *dynamic_cast<EXMeasLayer *>(fCradlePtr->At(lyr));
	  EXMeasLayer* curml = &ml; 

#ifdef _DEBUG_
	  std::cout << "layer: " << lyr << std::endl;
#endif

	  while(1)
	  { 
#ifdef _DEBUG_
	  std::cout << "in while loop " << std::endl;
#endif
		  rktrk.StepRungeKutta(STEP);
     	  currentPosition = rktrk.GetCurPosition();
    	  currentMomentum = rktrk.GetCurMomentum();
#ifdef _DEBUG_
		  std::cout << "current r: " << currentPosition.Perp() << std::endl;
#endif
    
		  if(lastml!=0 && lastml->CalcS(lastPosition)*lastml->CalcS(currentPosition)<0.)
		  {
			  //loop back
			  dlyr = -1;


			  rktrk.SetCurPosition(lastPosition); 
			  rktrk.SetCurMomentum(lastMomentum);
#ifdef _DEBUG_
    		  std::cout << "loop back" << std::endl;
			  std::cout << "last r: " << lastPosition.Perp() << std::endl;
			  std::cout << "current r: " << currentPosition.Perp() << std::endl;
#endif
			  //break;
			  curml = lastml;
		  } 


    	  if(curml->CalcS(currentPosition)*curml->CalcS(lastPosition)<0.)
    	  {
#ifdef _DEBUG_
    		  std::cout << "Found a root region" << std::endl;
			  std::cout << "last r: " << lastPosition.Perp() << std::endl;
			  std::cout << "current r: " << currentPosition.Perp() << std::endl;
#endif
			  // 
			  //get the exact root
			  //
			  Double_t step = STEP/2.;

			  TVector3 xa = lastPosition;
			  TVector3 pa = lastMomentum;

			  TVector3 xc, pc;

			  const Int_t ITER = 100;
			  for(Int_t n=0; n<ITER; n++)
			  {
				  rktrk.SetCurPosition(xa);
				  rktrk.SetCurMomentum(pa);

				  rktrk.StepRungeKutta(step);
				  xc = rktrk.GetCurPosition();
				  pc = rktrk.GetCurMomentum();

				  if(curml->CalcS(xa)*curml->CalcS(xc)>0.)
				  {
					  xa = xc;
					  pa = pc;
				  }

				  step /= 2.;
			  }

#ifdef _DEBUG_
			  std::cout << "Crossing point:("
				   << xc.X() << ","
				   << xc.Y() << ","
				   << xc.Z() << "), " 
				   << "rho: " << xc.Perp() 
				   << std::endl;
#endif
			  if (ml.IsActive())
			  {	         
				  ml.ProcessHit(xc, *fHitBufPtr); // create hit point
			  }

			  //restore
			  rktrk.SetCurPosition(currentPosition);
			  rktrk.SetCurMomentum(currentMomentum);

    	      lastPosition = currentPosition;
    	      lastMomentum = currentMomentum;
			  
			  break;
    	  }
    	  else
    	  {
    	      lastPosition = currentPosition;
    	      lastMomentum = currentMomentum;
    	  }//if
	  }//while

	  lastml = &ml; 
   }
}

void EXRungeKuttaEventGen::DumpHits()
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
