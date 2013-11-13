
#include "EXVTXKalDetector.h"
#include "EXVTXMeasLayer.h"
#include "EXVTXHit.h"
#include "TRandom.h"
#include "TMath.h"
#include <sstream>

ClassImp(EXVTXKalDetector)

EXVTXKalDetector::EXVTXKalDetector(Int_t m)
                : EXVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("VTXAir", "", A, Z, density, radlen, 0.);

   A       = 28.0855;
   Z       = 14.;
   density = 2.33;   
   radlen  = 9.36;
   TMaterial &si = *new TMaterial("VTXSi", "", A, Z, density, radlen, 0.);
   
   Bool_t active = EXVTXMeasLayer::kActive;
   Bool_t dummy  = EXVTXMeasLayer::kDummy;
   Bool_t activeness[2] = { active, dummy};
   
//====================================   
//     geometry of 3 doublet vtx 
//====================================

   InitGeometry();

   static const Double_t sigmaxi   = 1.44e-4;
   static const Double_t sigmazeta = 1.44e-4;
   
   for (Int_t layer = 0; layer < _nLayer; layer++) {
      Int_t    nladder   = _geodata[layer].nladder;
      Double_t rmin      = _geodata[layer].rmin;
      Double_t sximin    = _geodata[layer].sximin;
      Double_t sximax    = _geodata[layer].sximax;
      Double_t sxiwidth  = sximax-sximin;
      Double_t xioffset  = (sximin+sximax)/2;
      Double_t hlength   = _geodata[layer].hlength;
      Double_t lyrthick  = _geodata[layer].lyrthick;
      Double_t plyhwidth = rmin*TMath::Tan( TMath::Pi() / nladder );
      static const Double_t eps = 1e-6;
     
      for (Int_t ladder = 0; ladder < nladder; ladder++) {
         Double_t cosphi  =_geodata[layer].cosphi[ladder];
         Double_t sinphi  =_geodata[layer].sinphi[ladder];
       
         TVector3 xc( rmin*cosphi, rmin*sinphi, 0);
         TVector3 normal( cosphi, sinphi, 0);
         std::stringstream ss;
         ss << "VTX" << "_ly" << layer << "_ld" << ladder << std::ends;

         Bool_t inner = activeness[layer%2];
         Bool_t outer = activeness[(layer+1)%2];

         if (layer%2 == 0) { // Side drop of ladder0 is defined after the last ladder,
            if (ladder == 0) {   // bacause Side drop of ladder0 is outer than the last ladder.

               Add(new EXVTXMeasLayer(air, si, xc, normal, rmin + 2*ladder*eps, plyhwidth - sximin, 2*hlength, (plyhwidth + sximin)/2, sigmaxi, sigmazeta, inner,ss.str().data()));
               xc.SetXYZ( (rmin+lyrthick)*cosphi, (rmin+lyrthick)*sinphi, 0);
               Add(new EXVTXMeasLayer(si, air, xc, normal, rmin + (2*ladder+1)*eps, plyhwidth - sximin, 2*hlength, (plyhwidth + sximin)/2, sigmaxi, sigmazeta, outer));
	   
               xc.SetXYZ( rmin*cosphi, rmin*sinphi, 0);
               Add(new EXVTXMeasLayer(air, si, xc, normal, rmin + 2*nladder*eps, sxiwidth + sximin - plyhwidth, 2*hlength, (sximax + plyhwidth)/2, sigmaxi, sigmazeta, inner,ss.str().data()));
               xc.SetXYZ( (rmin+lyrthick)*cosphi, (rmin+lyrthick)*sinphi, 0);
               Add(new EXVTXMeasLayer(si, air, xc, normal, rmin + (2*nladder+1)*eps, sxiwidth + sximin -plyhwidth, 2*hlength, (sximax + plyhwidth)/2, sigmaxi, sigmazeta, outer));
            } else {
	   
               Add(new EXVTXMeasLayer(air, si, xc, normal, rmin + 2*ladder*eps, sxiwidth, 2*hlength, xioffset, sigmaxi, sigmazeta, inner,ss.str().data()));
               xc.SetXYZ( (rmin+lyrthick)*cosphi, (rmin+lyrthick)*sinphi, 0);
               Add(new EXVTXMeasLayer(si, air, xc, normal, rmin + (2*ladder+1)*eps, sxiwidth, 2*hlength, xioffset, sigmaxi, sigmazeta, outer));
            }
         } else {
            Add(new EXVTXMeasLayer(air, si, xc, normal, rmin + 2*ladder*eps, sxiwidth, 2*hlength, xioffset, sigmaxi, sigmazeta, inner,ss.str().data()));
            xc.SetXYZ( (rmin+lyrthick)*cosphi, (rmin+lyrthick)*sinphi, 0);
            Add(new EXVTXMeasLayer(si, air, xc, normal, rmin + (2*ladder+1)*eps, sxiwidth, 2*hlength, xioffset, sigmaxi, sigmazeta, outer));
         }
      }
   }
   SetOwner();			
}

EXVTXKalDetector::~EXVTXKalDetector()
{
}

// ---------------
//  Geometry Data
// ---------------
void EXVTXKalDetector::InitGeometry() // Geometry of 3 doublet vtx
{
    // Save frequently used parameters.
    
    Int_t _maxLadder = 17;
    _geodata.resize(_nLayer);
    
    Double_t nladder[_nLayer]   = { 10, 10, 11, 11, 17, 17}; // Number of ladders in this layer
    Double_t rmin[_nLayer]      = { 1.6, 1.8, 3.7, 3.9, 5.8, 6.0};  // Distance of sensitive area from
    Double_t xiwidth[_nLayer]   = { 1.1, 1.1, 2.2, 2.2, 2.2, 2.2};
    Double_t zetawidth[_nLayer] = { 12.5, 12.5, 25, 25, 25, 25};
    Double_t xioffset[_nLayer]  = { -0.1575320104, -0.1575320104, -0.1531923469, -0.1531923469, -0.2293817688, -0.2293817688};
    
    
    for (Int_t ly=0;ly<_nLayer;ly++) {
        if ( _maxLadder < _geodata[ly].nladder ) { _maxLadder = _geodata[ly].nladder; }
        
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
        for (Int_t ld = 0;ld<_geodata[ly].nladder;ld++) {
            Double_t phi = (Double_t)(_geodata[ly].phi0 + _geodata[ly].dphi*ld);
            _geodata[ly].cosphi[ld] = cos(phi);
            _geodata[ly].sinphi[ld] =-sin(phi); // direction of z axis is fliped
        }
    }
}
