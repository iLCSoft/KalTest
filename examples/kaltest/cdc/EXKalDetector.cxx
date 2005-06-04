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
   // create inner timing detector
   if (rt0det > 0.) {
         TVector3 we(0., rt0det, -lhalf);
         TVector3 wd(0., 0., 1.);
         Double_t celw = 2*celhw*rt0det/rmin;
         Add(new EXMeasLayer(rt0det, lhalf, we, wd.Unit(), celw, kgX0Si, kgX0CFRP));
         Add(new EXMeasLayer(rt0det+4.5, lhalf, we, wd.Unit(), celw, kgX0CFRP, kgX0Inv));
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
         Add(new EXMeasLayer(r0, lhalf, we, wd.Unit(), celw, kgX0Inv, kgX0Inv));
         r += rstep;
      }
   }
   SetOwner();
}

Double_t EXKalDetector::CalcRadLength(const TVMeasLayer &from,
                                      const TVMeasLayer &to) const
{
   Double_t rfrom = dynamic_cast<const THype &>(from).GetR0();
   Double_t rto   = dynamic_cast<const THype &>(to).GetR0();
   static const Double_t kTol = 1.e-9;
   if ((rfrom - fRT0Det - kTol) * (rto - fRT0Det - kTol) <= 0.) return 1./kgX0CFRP;
   return 1./kgX0Inv; 
}

Double_t EXKalDetector::CalcSigmaMS0(const TVMeasLayer &vfrom,
                                     const TVMeasLayer &vto,
                                           Double_t     pathlen) const
{
   static const Double_t kMS1  = 0.0136;
   static const Double_t kMS2  = 0.038;

   const EXMeasLayer &from = *dynamic_cast<const EXMeasLayer *> (&vfrom);
   const EXMeasLayer &to   = *dynamic_cast<const EXMeasLayer *> (&vto);

   Double_t xl    = pathlen / CalcRadLength(from,to);
   // ------------------------------------------------------------------
            xl    = TMath::Max(xl,1.e-4); // very crude treatment that
					  // should be improved!
   // ------------------------------------------------------------------

   // ------------------------------------------------------------------
   if (fRT0Det > 0.) {			// VERY TEMPORARY
      Double_t rfrom = from.GetWireEnd().Perp();
      Double_t rto   = to.GetWireEnd().Perp();
      Double_t eps   = 1.e-2;
      if ((fRT0Det+eps - rfrom)*(fRT0Det+eps - rto) <= 0.) {
#if 1
         xl = 0.006;  // CFRP or Field cage
#else
         xl = 0.03;  // CFRP or Field cage
#endif
      }
   }
   // ------------------------------------------------------------------

   Double_t sgms0 = kMS1 * TMath::Sqrt(xl) * (1 + kMS2*TMath::Log(xl));

   if (TKalDetCradle::GetInstance().IsMSOn()) return sgms0;
   else          return 0.;
}
