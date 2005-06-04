
#include "EXITKalDetector.h"
#include "EXITMeasLayer.h"
#include "EXITHit.h"
#include "TRandom.h"

ClassImp(EXITKalDetector)

EXITKalDetector::EXITKalDetector(Int_t m)
               : EXVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("ITAir", "", A, Z, density, radlen, 0.);

   A       = 28.0855;
   Z       = 14.;
   density = 2.33;
   radlen  = 9.36;
   TMaterial &si = *new TMaterial("ITSi", "", A, Z, density, radlen, 0.);

   static const Int_t    nlayers  = 4;
   static const Double_t lhalfmin = 18.5;      // layer min half length [cm]
   static const Double_t rmin     = 9.;        // layer radius min
   static const Double_t rstep    = 7.;        // layer radius step
   static const Double_t lstep    = 14.5;      // layer length step
   static const Double_t thick    = 0.05616;   // layer thick

   static const Double_t sigmax   = 1.e-3;
   static const Double_t sigmaz   = 1.e-3;
                                                                                
   Bool_t active = EXITMeasLayer::kActive;
   Bool_t dummy  = EXITMeasLayer::kDummy;

   // create measurement layers of inner tracker
   Double_t r   = rmin;
   Double_t len = lhalfmin;
   for (Int_t layer = 0; layer < nlayers; layer++) {
      Add(new EXITMeasLayer(air, si, r, len, sigmax, sigmaz, active));
      Add(new EXITMeasLayer(si, air, r + thick, len, sigmax, sigmaz, dummy));
      r   += rstep;
      len += lstep;
   }
   SetOwner();
}

EXITKalDetector::~EXITKalDetector()
{
}

void EXITKalDetector::ProcessHit(const TVector3    &xx,
                                 const TVMeasLayer &ms,
                                       TObjArray   &hits)
{
   const EXITMeasLayer &ims = dynamic_cast<const EXITMeasLayer &>(ms);
   TKalMatrix h    = ims.XvToMv(xx);
   Double_t   rphi = h(0, 0);
   Double_t   z    = h(1, 0);

   Double_t dx = ims.GetSigmaX();
   Double_t dz = ims.GetSigmaZ();
   rphi += gRandom->Gaus(0., dx);   // smearing rphi
   z    += gRandom->Gaus(0., dz);   // smearing z

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = rphi;
   meas [1] = z;
   dmeas[0] = dx;
   dmeas[1] = dz;

   hits.Add(new EXITHit(ims, meas, dmeas, xx, GetBfield(xx)));
}
