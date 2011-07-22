
#include "EXITKalDetector.h"
#include "EXITMeasLayer.h"
#include "EXITFBMeasLayer.h"
#include "TRandom.h"
#include <sstream>

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

  ////////// Barrel part of IT /////////////////////

  static const Int_t    nlayers  = 4;
  static const Double_t lhalfmin = 18.5;      // layer min half length [cm]
  static const Double_t rmin     = 9.;        // layer radius min
  static const Double_t rstep    = 7.;        // layer radius step
  static const Double_t lstep    = 14.5;      // layer length step
  static const Double_t thick    = 0.05616;   // layer thick
  static const Double_t sigmax   = 1.e-3;
  static const Double_t sigmaz   = 1.e-3;

  ///////// Forward & Backward Part of IT //////////

  static const Int_t    nlayersfb    = 7;
  static const Int_t    nlyrsfbdmmy  = 10;
  static const Double_t zminfb       = 14.5;      // layer located z-direction [cm]
  static const Double_t zstep        = 14.5;      // layer step
  static const Double_t rinfb        = 2.5;       // layer inner radius
  static const Double_t routfb       = 7.;        // layer outer radius

  Bool_t active = EXVMeasLayer::kActive;
  Bool_t dummy  = EXVMeasLayer::kDummy;

  // create measurement layers of inner tracker
  Double_t r    = rmin;
  Double_t len  = lhalfmin;
  Double_t rifb = rinfb;
  Double_t zfb  = zminfb;
  Double_t rofb = routfb;

  for (Int_t layer = 0; layer < (nlayersfb + nlyrsfbdmmy); layer++) {
    TVector3 xc1(0,0,zfb);         // for forward part
    TVector3 xc2(0,0,zfb+thick);   // for forward part
    TVector3 xc3(-xc1);          // for backward part
    TVector3 xc4(-xc2);          // for backward part

    if (layer < nlayersfb) {
      //  Forward part
      std::stringstream ss;
      std::stringstream ssdm;
      ss << "ITF" << layer << std::ends;
      ssdm << "ITFdm" << layer << std::ends;
      Add(new EXITFBMeasLayer(air, si, xc1, rifb, rofb, sigmax, sigmaz, active,ss.str().data()));
      Add(new EXITFBMeasLayer(si, air, xc2, rifb, rofb, sigmax, sigmaz, dummy, ssdm.str().data()));
      //  Backward part
      ss.clear();
      ssdm.clear();
      ss << "ITB" << layer << std::ends;
      ssdm << "ITBdm" << layer << std::ends;
      Add(new EXITFBMeasLayer(air, si, xc3, rifb, rofb, sigmax, sigmaz, active,ss.str().data()));
      Add(new EXITFBMeasLayer(si, air, xc4, rifb, rofb, sigmax, sigmaz, dummy, ssdm.str().data()));
    } else {      
      Add(new EXITFBMeasLayer(air, air, xc1, rifb, rofb, sigmax, sigmaz, dummy));
      Add(new EXITFBMeasLayer(air, air, xc3, rifb, rofb, sigmax, sigmaz, dummy));
    }
    if (layer < nlayers) {
      std::stringstream ss;
      ss << "IT" << layer << std::ends;
      Add(new EXITMeasLayer(air, si, r, len, sigmax, sigmaz, active, ss.str().data()));
      Add(new EXITMeasLayer(si, air, r + thick, len, sigmax, sigmaz, dummy));
    }
    len += lstep;
    zfb += zstep;
    r   += rstep;

    if      (layer == 1) rifb += 2.5;
    else if (layer == 2) rifb += 1.;
    else if (layer >= 6) rifb += 0.;
    else                 rifb += 2.;

    // layer should be bigger than its inner layer
    // because Rout is used as a sorting policy.
    if      (layer == 3) rofb += 8.;
    else if (layer == 4) rofb += 2.;
    else if (layer >= 5) rofb += 1.e-6;
    else                 rofb += 7.;
  }
  SetOwner();
}

EXITKalDetector::~EXITKalDetector()
{
}
