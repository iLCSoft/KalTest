
#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "EXHit.h"
#include "TRandom.h"

Double_t EXKalDetector::fgBfield = 30.;

ClassImp(EXKalDetector)

EXKalDetector::EXKalDetector(Int_t m)
             : TVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("Air", "", A, Z, density, radlen, 0.);

   A       = 12.0107;
   Z       =  6.;
   density = 0.1317;
   radlen  = 42.7/density;
   TMaterial &cfrp = *new TMaterial("CFRP", "", A, Z, density, radlen, 0.);

   static const Int_t    nlayers   = 50;
   static const Double_t lhalf     = 200.;
   static const Double_t rmin      = 45.;
   static const Double_t rstep     = 3.;
   static const Double_t rcylin    = 43.;
   static const Double_t rcylout   = 44.;

   // create dummy layers of inner cylinder of central tracker
   Add(new EXMeasLayer(air, cfrp, rcylin, lhalf, EXMeasLayer::kDummy));
   Add(new EXMeasLayer(cfrp, air, rcylout, lhalf, EXMeasLayer::kDummy));

   // create measurement layers of central tracker
   Double_t r   = rmin;
   for (Int_t layer = 0; layer < nlayers; layer++) {
      Add(new EXMeasLayer(air, air, r, lhalf, EXMeasLayer::kActive));
      r += rstep;
   }
   SetOwner();
}

EXKalDetector::~EXKalDetector()
{
}

void EXKalDetector::ProcessHit(const TVector3    &xx,
                               const EXMeasLayer &ms,
                                     TObjArray   &hits)
{
   TKalMatrix h    = ms.XvToMv(xx);
   Double_t   rphi = h(0, 0);
   Double_t   z    = h(1, 0);

   Double_t dx = ms.GetSigmaX();
   Double_t dz = ms.GetSigmaZ();
   rphi += gRandom->Gaus(0., dx);   // smearing rphi
   z    += gRandom->Gaus(0., dz);   // smearing z
                                                                          
   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = rphi;
   meas [1] = z;
   dmeas[0] = dx;
   dmeas[1] = dz;

   hits.Add(new EXHit(ms, meas, dmeas, xx, GetBfield(xx)));
}
