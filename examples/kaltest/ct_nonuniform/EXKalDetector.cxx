#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "EXHit.h"

#include "TRandom.h" // from ROOT
#include "TVirtualPad.h"
#include "TTUBE.h"
#include "TTUBS.h"
#include "TNode.h"

#include <iostream>

Double_t EXKalDetector::fgBfield = 10.;

ClassImp(EXKalDetector)

EXKalDetector::EXKalDetector(Int_t m)
             : EXVKalDetector(m),fNodePtr(0)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;    // mass number
   Z       = 7.3;                               // atomic number
   density = 1.205e-3;                          // [g/cmm^3]
   radlen  = 3.42e4;                            // [cm]
   TMaterial &air = *new TMaterial("Air", "", A, Z, density, radlen, 0.);

   A       = 12.0107;                           // mass number
   Z       =  6.;                               // atomic number
   density = 0.1317;                            // [g/cmm^3]
   radlen  = 42.7/density;                      // [cm]
   //TMaterial &cfrp = *new TMaterial("CFRP", "", A, Z, density, radlen, 0.);

   static const Int_t    nlayers   = 251;        // # sampling layers 
   static const Double_t lhalf     = 2000.;      // half length
   static const Double_t rmin      = 300.;         // r_{min} = radius of 0th layer
   static const Double_t rstep     = 6;         // step in r
   //static const Double_t rcylin    = 43.;       // inner radius of CFRP cylinder 
   //static const Double_t rcylout   = 44.;       // outer radius of CFRP cylinder

   // Create dummy layers of the inner cylinder of the central tracker
   //Add(new EXMeasLayer(air, cfrp, rcylin, lhalf, EXMeasLayer::kDummy));
   //Add(new EXMeasLayer(cfrp, air, rcylout, lhalf, EXMeasLayer::kDummy));

   // Create measurement layers of the central tracker
   Double_t r   = rmin;
   for (Int_t layer = 0; layer < nlayers; layer++) {

      if(layer%1==0)
      {   
          Add(new EXMeasLayer(air, air, r, lhalf, EXMeasLayer::kActive));
      }   
      else
      {   
          Add(new EXMeasLayer(air, air, r, lhalf, EXMeasLayer::kDummy));
      }

      r += rstep;
   }

   SetOwner();
}

EXKalDetector::~EXKalDetector()
{
}

// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXKalDetector::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   TNode *nodep = GetNodePtr();
   nodep->cd();

   if (!fNodePtr) {
      EXMeasLayer *inp  = static_cast<EXMeasLayer *>(First());
      EXMeasLayer *outp = static_cast<EXMeasLayer *>(Last());
      Double_t rin  = inp->GetR();
      Double_t rout = outp->GetR();
      Double_t hlen = outp->GetZmax();
      const Char_t *name  = "DET";
      const Char_t *nname = "DETNode";
      TTUBE *tubep = new TTUBE(name,name,"void",rin,rout,hlen);
      tubep->SetBit(kCanDelete);
      fNodePtr = new TNode(nname,nname,name);
      fNodePtr->SetLineColor(color);
      fNodePtr->SetLineWidth(0.01);
   }

   EXVKalDetector::Draw(color,opt);
}
