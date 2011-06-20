#ifndef __EXVTXDETECTOR__
#define __EXVTXDETECTOR__

#include "EXVKalDetector.h"
#include "TMath.h"

class EXVTXKalDetector : public EXVKalDetector {
  
 inline void InitGeometry();
  
 public:
  EXVTXKalDetector(Int_t m = 100);
  ~EXVTXKalDetector();
  
  ClassDef(EXVTXKalDetector,1)   // Sample hit class
    
    
    static const Int_t _nLayer = 6;
  struct GeoData_t {
    Int_t nladder;
    Float_t rmin;  // distance of inner surface of sensitive region from IP
    Float_t dphi;  // azimuthal angle step of each ladder
    Float_t phi0;  // aximuthal angle offset
    std::vector<Double_t> cosphi;  // cos[phi_ladder], cos_phi of each ladder
    std::vector<Double_t> sinphi;  // sin[phi_ladder], sin_phi of each ladder
    Float_t lyrthick;  // sensitive region thickness
    Float_t xiwidth; // ladder width in xi
    Float_t sximin;  // minimum xi of sensitive region.
    Float_t sximax;  // maximum xi of sensitive region
    Float_t hlength; // ladder's half length in z
  };
  std::vector<GeoData_t> _geodata;
  
  
};

void EXVTXKalDetector::InitGeometry() // Geometry of 3 doublet vtx
{
  // Save frequently used parameters.
  
  Int_t _maxLadder = 17;
  _geodata.resize(_nLayer);
  
  Float_t nladder[_nLayer]   = { 10, 10, 11, 11, 17, 17}; // Number of ladders in this layer
  Float_t rmin[_nLayer]      = { 1.6, 1.8, 3.7, 3.9, 5.8, 6.0};  // Distance of sensitive area from
  Float_t xiwidth[_nLayer]   = { 1.1, 1.1, 2.2, 2.2, 2.2, 2.2};
  Float_t zetawidth[_nLayer] = { 12.5, 12.5, 25, 25, 25, 25};
  Float_t xioffset[_nLayer]  = { -0.1575320104, -0.1575320104, -0.1531923469, -0.1531923469, -0.2293817688, -0.2293817688};


  for(Int_t ly=0;ly<_nLayer;ly++){
    if( _maxLadder < _geodata[ly].nladder ) { _maxLadder = _geodata[ly].nladder; }

    _geodata[ly].nladder = nladder[ly];
    _geodata[ly].rmin    = rmin[ly];
    _geodata[ly].phi0    = 1.570796327; // phi offset.
    _geodata[ly].dphi    = 2*TMath::Pi()/_geodata[ly].nladder;
    _geodata[ly].xiwidth = xiwidth[ly];
    _geodata[ly].sximin  = -xioffset[ly] - xiwidth[ly]/2.0;
    _geodata[ly].sximax  = -xioffset[ly] + xiwidth[ly]/2.0;
    _geodata[ly].hlength = zetawidth[ly]/2;
    _geodata[ly].lyrthick= 5e-3;
    _geodata[ly].cosphi.resize( _geodata[ly].nladder );
    _geodata[ly].sinphi.resize( _geodata[ly].nladder );
    for(Int_t ld = 0;ld<_geodata[ly].nladder;ld++) {
      Double_t phi = (Double_t)(_geodata[ly].phi0 + _geodata[ly].dphi*ld);
      _geodata[ly].cosphi[ld] = cos(phi);
      _geodata[ly].sinphi[ld] =-sin(phi); // direction of z axis is fliped
    }
  }
}

#endif
