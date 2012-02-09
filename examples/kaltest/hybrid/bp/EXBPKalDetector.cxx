
#include "EXBPKalDetector.h"
#include "EXBPMeasLayer.h"
#include "EXBPConeMeasLayer.h"
#include "EXBPHit.h"
#include "TRandom.h"
#include <sstream>

ClassImp(EXBPKalDetector)

EXBPKalDetector::EXBPKalDetector(Int_t m)
               : EXVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-8; // [g/cm^3]
   radlen  = 3.42e9;   // [cm]
   TMaterial &vacuum = *new TMaterial("BPVacuum", "", A, Z, density, radlen, 0.);
    
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3; // [g/cm^3]
   radlen  = 3.42e4;   // [cm]
   TMaterial &air = *new TMaterial("BPAir", "", A, Z, density, radlen, 0.);

   A       = 9.012182;
   Z       = 4.;
   density = 1.848; // [g/cm^3] 
   radlen  = 42.70/density; // [cm]
   TMaterial &be = *new TMaterial("BPBe", "", A, Z, density, radlen, 0.);

#if 1
   static const Double_t lhalf     = 17.5;
#else
   static const Double_t lhalf     = 0.01; // very short Be beam pipe for testing
#endif
   static const Double_t rmin      =  1.5;
   static const Double_t thickness = 0.05;

   static const Double_t z1        = lhalf;
   static const Double_t r1        = rmin;
   static const Double_t z2        = 300.;
   static const Double_t r2        =  20.;
    
   static const Double_t sigmax    = 4.e-4;
   static const Double_t sigmaz    = 4.e-4;

   Bool_t active = EXBPMeasLayer::kActive;
   Bool_t dummy  = EXBPMeasLayer::kDummy;

   // create measurement layers of central tracker
   std::stringstream ss;
   ss << "BPC" << std::ends;
   Double_t r   = rmin;
   Double_t len = lhalf;
   Add(new EXBPMeasLayer(vacuum, be , r            , len, sigmax, sigmaz, active, ss.str().data()));
   Add(new EXBPMeasLayer(be    , air, r + thickness, len, sigmax, sigmaz, dummy));
   ss.str("");
   ss.clear();
   ss << "BPMZ" << std::ends;
   Add(new EXBPConeMeasLayer(vacuum,  be, -z1, r1,       -z2, r2      , sigmax, sigmaz, active, ss.str().data()));
   Add(new EXBPConeMeasLayer(be    , air, -z1, r1+thickness, -z2, r2+thickness, sigmax, sigmaz, dummy));
   ss.str("");
   ss.clear();
   ss << "BPPZ" << std::ends;
   Add(new EXBPConeMeasLayer(vacuum,  be, +z1, r1      , +z2, r2      , sigmax, sigmaz, active, ss.str().data()));
   Add(new EXBPConeMeasLayer(be    , air, +z1, r1+thickness, +z2, r2+thickness, sigmax, sigmaz, dummy));
   SetOwner();
}

EXBPKalDetector::~EXBPKalDetector()
{
}
