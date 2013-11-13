#ifndef EXVTXDETECTOR_H
#define EXVTXDETECTOR_H

#include "EXVKalDetector.h"
#include "TMath.h"

class EXVTXKalDetector : public EXVKalDetector {
public:
   EXVTXKalDetector(Int_t m = 100);
  ~EXVTXKalDetector();

private:
   void InitGeometry();
   static const Int_t _nLayer = 6;
   struct GeoData_t {
     Int_t nladder;
     Double_t rmin;  // distance of inner surface of sensitive region from IP
     Double_t dphi;  // azimuthal angle step of each ladder
     Double_t phi0;  // aximuthal angle offset
     std::vector<Double_t> cosphi;  // cos[phi_ladder], cos_phi of each ladder
     std::vector<Double_t> sinphi;  // sin[phi_ladder], sin_phi of each ladder
     Double_t lyrthick;  // sensitive region thickness
     Double_t xiwidth; // ladder width in xi
     Double_t sximin;  // minimum xi of sensitive region.
     Double_t sximax;  // maximum xi of sensitive region
     Double_t hlength; // ladder's half length in z
   };
   std::vector<GeoData_t> _geodata;
  
   ClassDef(EXVTXKalDetector,1)   // Sample hit class
};
#endif
