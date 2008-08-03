
#include "EXTPCKalDetector.h"
#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "TRandom.h"
#include "TMath.h"

#include "TTUBE.h"
#include "TNode.h"
#include "TVirtualPad.h"

Double_t EXTPCKalDetector::fgVdrift = 5.e-3;
ClassImp(EXTPCKalDetector)

EXTPCKalDetector::EXTPCKalDetector(Int_t m)
                : EXVKalDetector(m),
		  fNodePtr(0)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("TPCAir", "", A, Z, density, radlen, 0.);

   A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
   Z       = 16.4;
   density = 0.749e-3;
   radlen  =  1.196e4*2;
   TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);

   A       = 12.0107;
   Z       =  6.;
   density = 0.1317;
   radlen  = 42.7/density;
   TMaterial &cfrp = *new TMaterial("TPCCFRP", "", A, Z, density, radlen, 0.);

   static const Double_t inthick   = 2.1075 * 2.;   // thick of inner shell
   static const Double_t outthick  = 4.1175 * 2.;   // thick of outer shell
   static const Int_t    nlayers   = 200;           // number of layer
   static const Double_t lhalf     = 255.;          // half length
   static const Double_t rmin      = 44.215;        // minimum radius
   static const Double_t rstep     = 0.76775;       // step length of radius
   static const Double_t rtub      = 39.5;          // inner r of support tube
   static const Double_t outerr    = 206.;          // outer radius of TPC
   static const Double_t sigmax0   = 55.e-4;
   static const Double_t sigmax1   = 166.e-4 / 3 / TMath::Sqrt(28);
   static const Double_t sigmaz    = 600.e-4;

   Bool_t active = EXTPCMeasLayer::kActive;
   Bool_t dummy  = EXTPCMeasLayer::kDummy;
   Add(new EXTPCMeasLayer(air, cfrp, rtub, lhalf, sigmax0, sigmax1, sigmaz, dummy));
   Add(new EXTPCMeasLayer(cfrp, gas, rtub+inthick, lhalf, sigmax0, sigmax1, sigmaz, dummy));

   // create measurement layers of central tracker
   Double_t r = rmin;
   for (Int_t layer = 0; layer < nlayers; layer++) {
      Add(new EXTPCMeasLayer(gas, gas, r, lhalf, sigmax0, sigmax1, sigmaz, active));
      r += rstep;
   }
   Add(new EXTPCMeasLayer(gas, cfrp, outerr-outthick, lhalf, sigmax0, sigmax1, sigmaz, dummy));
   Add(new EXTPCMeasLayer(cfrp, air, outerr, lhalf, sigmax0, sigmax1, sigmaz, dummy));

   SetOwner();
}

EXTPCKalDetector::~EXTPCKalDetector()
{
}

// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXTPCKalDetector::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   TNode *nodep = GetNodePtr();
   nodep->cd();

   if (!fNodePtr) {
      EXTPCMeasLayer *inp  = static_cast<EXTPCMeasLayer *>(First());
      EXTPCMeasLayer *outp = static_cast<EXTPCMeasLayer *>(Last());
      Double_t rin  = inp->GetR();
      Double_t rout = outp->GetR();
      Double_t hlen = outp->GetZmax();
      const Char_t *name  = "TPC";
      const Char_t *nname = "TPCNode";
      TTUBE *tubep = new TTUBE(name,name,"void",rin,rout,hlen);
      tubep->SetBit(kCanDelete);
      fNodePtr = new TNode(nname,nname,name);
      fNodePtr->SetLineColor(color);
      fNodePtr->SetLineWidth(0.01);
   }
   EXVKalDetector::Draw(color,opt);
}
