//-------------------
// VERY TEMPORARY
//-------------------
//#define __LRBYWIRE__
//-------------------
#define __REALISTIC__

#include "TKalDetCradle.h"
#include "EXKalDetector.h"

ClassImp(EXKalDetector)

// ------------------------------------------------------------------
// A very crude model of JLC-CDC
// ------------------------------------------------------------------
//

Double_t kgX0Inv  = 1. / 1.8182e4;        // CO2:C4H10 (90:10)
Double_t kgX0Si   = 1. / 9.36;           // Silicon
Double_t kgX0CFRP = 1. / (42.7/0.1317);  // CFRP

EXKalDetector::EXKalDetector(Int_t    nlayers,
                             Int_t    nwires,
                             Double_t lhalf,
                             Double_t celhw,
                             Double_t tna,
                             Double_t rmin,
                             Double_t rstep,
                             Double_t rt0det)
	: TVKalDetector(nlayers+1), fRT0Det(rt0det)
{
   // prepare materials
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;    // mass number
   Z       = 7.3;                               // atomic number
   density = 1.205e-3;                          // [g/cmm^3]
   radlen  = 3.42e4;                            // [cm]
   TMaterial &air = *new TMaterial("Air", "", A, Z, density, radlen, 0.);

   A       = 28.0855;
   Z       = 14.;
   density = 2.33;
   radlen  = 9.36;
   TMaterial &si = *new TMaterial("Si", "", A, Z, density, radlen, 0.);
   static const Double_t thick    = 0.05616;    // layer thickness

   A       = 12.0107;                           // mass number
   Z       =  6.;                               // atomic number
   density = 0.1317;                            // [g/cmm^3]
   radlen  = 42.7/density;                      // [cm]
   TMaterial &cfrp = *new TMaterial("CFRP", "", A, Z, density, radlen, 0.);

   static Bool_t kActive = kTRUE;
   static Bool_t kDummy  = kFALSE;

   // create inner timing detector
   if (rt0det > 0.) {
         TVector3 we(0., rt0det, -lhalf);
         TVector3 wd(0., 0., 1.);
         Double_t celw = 2*celhw*rt0det/rmin;
         Add(new EXMeasLayer(rt0det      , lhalf, we, wd.Unit(), celw, air, si, kActive));
         Add(new EXMeasLayer(rt0det+thick, lhalf, we, wd.Unit(), celw, si, air, kDummy));
         
         Add(new EXMeasLayer(rt0det+4.5, lhalf, we, wd.Unit(), celw, air, cfrp, kDummy));
         Add(new EXMeasLayer(rt0det+5.0, lhalf, we, wd.Unit(), celw, cfrp, air, kDummy));
   }
   // create measurement layers of central tracker

   Double_t r    = rmin;
   Double_t dfiw = (celhw/2)/rmin;
   for (Int_t slayer=0; slayer<nlayers/nwires; slayer++) {
#ifndef __REALISTIC__
      Int_t sslayer = TMath::Min(slayer/3,2);
#endif
      Double_t tana;
#ifndef __REALISTIC__
      Double_t dfiwr = dfiw;
#else
      Int_t    ncels = (Int_t)(TMath::Pi()*r/celhw);
      Double_t dfiwr = 0.5*TMath::Pi()/ncels;
#endif
#ifndef __LRBYWIRE__
      dfiw   *= -dfiwr/dfiw;
#endif
      switch (slayer%3) {
         case 0:
                 tana = 0.;
                 break;
         case 1:
                 tana = -tna;
                 break;
         case 2:
         default:
                 tana = +tna;
      }
      for (Int_t wire=0; wire<nwires; wire++) {
#ifdef __LRBYWIRE__
         dfiw   *= -dfiwr/dfiw;
#endif
#ifndef __REALISTIC__
         Double_t fact = 1./(sslayer+1);
#else
         Double_t fact = 1.;
#endif
         Double_t celw = 4.*r*TMath::Abs(dfiw)*fact;
         TVector3 we(0., r, -lhalf);
         Double_t sndelphi = tana*lhalf/r;
         Double_t csdelphi = TMath::Sqrt((1-sndelphi)*(1+sndelphi));
         Double_t delphi   = TMath::ASin(sndelphi);
         we.RotateZ(fact*dfiw - delphi);
         TVector3 wd(-tana, 0., 1.);
         wd.RotateZ(fact*dfiw);

         Double_t r0 = r*csdelphi; // r0 = r(z=0) while r = r(z=-lhalf)
         Add(new EXMeasLayer(r0, lhalf, we, wd.Unit(), celw, air, air, kActive));
         r += rstep;
      }
   }
   SetOwner();
}
