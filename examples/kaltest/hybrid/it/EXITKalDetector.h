#ifndef __EXITDETECTOR__
#define __EXITDETECTOR__

#include "EXVKalDetector.h"
#include "TMath.h"

#include <iostream>
class EXITKalDetector : public EXVKalDetector {
  
  inline void InitFTDGeometry();
  
 public:
  EXITKalDetector(Int_t m = 100);
  ~EXITKalDetector();
  
  ClassDef(EXITKalDetector,1)   // Sample hit class
    
    static const Int_t _nDiskHalf = 7;
  
  struct FTDGeo_t {
    Int_t nPetal;
    Float_t dphi;
    Float_t phi0;
    Float_t alpha;
    Float_t cosAlpha;
    Float_t sinAlpha;
    std::vector<Double_t> cosphi;  // cos[phi_ladder], cos_phi of each ladder
    std::vector<Double_t> sinphi;  // sin[phi_ladder], sin_phi of each ladder
    Float_t lyrthick;
    Float_t rin;
    Float_t rout;
    Float_t spolicy;
    Float_t zCoord;
    Float_t dxMax;
    Float_t dxMin;
  };
  std::vector<FTDGeo_t> _FTDgeo;
  
};

void EXITKalDetector::InitFTDGeometry(){  // geometry of FTD
  // Save frequently used parameters.
  _FTDgeo.resize(_nDiskHalf);
  
  const Int_t npetal = 16;
  
  const Float_t zCoordinate[_nDiskHalf] = { 22.0, 37.1, 64.4, 104.6, 144.7, 184.8, 225.0}; // [cm]
  const Float_t       alpha[_nDiskHalf] = { 6, 6, 4, 4, 4, 4, 4}; //angle of petal to x-y plane  [degree]
  const Float_t         rin[_nDiskHalf] = { 3.9, 4.9, 7.0, 10.0, 13.0, 16.0, 19.0};
  const Float_t        rout[_nDiskHalf] = { 17.41, 17.34, 26.41, 30.25, 30.25, 30.24, 30.23};
  const Float_t     spolicy[_nDiskHalf] = { 17.41, 18.00, 26.41, 30.25, 31.00, 32.00, 33.00};
  const Float_t       dxMax[_nDiskHalf] = { 7.667, 7.667, 12.249, 12.249, 12.249, 12.249, 12.249};
  const Float_t       dxMin[_nDiskHalf] = { 2.294, 2.718, 4.527, 4.192, 5.388, 6.585, 7.782};
    
  const Float_t LayerThickness = 0.02; // [cm]
  
  for(Int_t disk=0; disk<_nDiskHalf; disk++){
    _FTDgeo[disk].nPetal   = npetal;
    _FTDgeo[disk].zCoord   = zCoordinate[disk];
    _FTDgeo[disk].alpha    = (Float_t)(alpha[disk]*TMath::Pi()/180);
    _FTDgeo[disk].cosAlpha = TMath::Cos(_FTDgeo[disk].alpha);
    _FTDgeo[disk].sinAlpha = TMath::Sin(_FTDgeo[disk].alpha);
    _FTDgeo[disk].rin      = rin[disk];
    _FTDgeo[disk].rout     = rout[disk];
    _FTDgeo[disk].spolicy  = spolicy[disk];
    _FTDgeo[disk].dxMax    = dxMax[disk];
    _FTDgeo[disk].dxMin    = dxMin[disk];
    _FTDgeo[disk].dphi     = (Float_t)(2*TMath::Pi()/_FTDgeo[disk].nPetal);
    _FTDgeo[disk].lyrthick = LayerThickness;    
    _FTDgeo[disk].cosphi.resize( _FTDgeo[disk].nPetal );
    _FTDgeo[disk].sinphi.resize( _FTDgeo[disk].nPetal );
    for(Int_t petal=0; petal<_FTDgeo[disk].nPetal; petal++){
      Double_t phi = (Double_t)(_FTDgeo[disk].dphi*petal);
      _FTDgeo[disk].cosphi[petal] = cos(phi);
      _FTDgeo[disk].sinphi[petal] = sin(phi);
    }
  }
}

#endif
