#ifndef EXITDETECTOR_H
#define EXITDETECTOR_H

#include "EXVKalDetector.h"
#include "TMath.h"

#include <iostream>
class EXITKalDetector : public EXVKalDetector {
public:
  EXITKalDetector(Int_t m = 100);
  ~EXITKalDetector();

private:
    void InitFTDGeometry();

private:
  static const Int_t _nDiskHalf = 7;
  struct FTDGeo_t {
    Int_t nPetal;
    Double_t dphi;
    Double_t phi0;
    Double_t alpha;
    Double_t cosAlpha;
    Double_t sinAlpha;
    std::vector<Double_t> cosphi;  // cos[phi_ladder], cos_phi of each ladder
    std::vector<Double_t> sinphi;  // sin[phi_ladder], sin_phi of each ladder
    Double_t lyrthick;
    Double_t rin;
    Double_t rout;
    Double_t spolicy;
    Double_t zCoord;
    Double_t dxMax;
    Double_t dxMin;
  };
  std::vector<FTDGeo_t> _FTDgeo;
  
  ClassDef(EXITKalDetector,1)   // Sample hit class
};
#endif
