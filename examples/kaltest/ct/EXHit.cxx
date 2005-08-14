#include "EXHit.h"
#include "EXMeasLayer.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXHit::EXHit(Int_t m)
     : TVTrackHit(m)
{
}
      
EXHit::EXHit(const EXMeasLayer &ms,    // measurement layer
                   Double_t    *x,     // coordinate array
                   Double_t    *dx,    // coordinate error array
                   Double_t     b,     // magnetic field
                   Int_t        m)     // dimension of meas. vector
     : TVTrackHit(ms, x, dx, b, m)
{
}

EXHit::~EXHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   const EXMeasLayer &ms
                 = dynamic_cast<const EXMeasLayer &>(GetMeasLayer());
   TKalMatrix h  = ms.XvToMv(xv);
   Double_t r    = ms.GetR();
   Double_t phih = h(0,0) / r;
   Double_t phim = (*this)(0,0) / r;
   Double_t dphi = phih - phim;

   static Double_t kPi    = TMath::Pi();
   static Double_t kTwoPi = 2 * kPi;

   while (dphi < -kPi) dphi += kTwoPi;
   while (dphi >  kPi) dphi -= kTwoPi;

   h(0,0)  = r * (phim + dphi);
   h(1,0) += 0.;

   return h;
}

void EXHit::DebugPrint(Option_t *) const
{
   cerr << "------------------- Site Info -------------------------" << endl;

   for (Int_t i=0; i<GetDimension(); i++) {
      Double_t x  = (*this)(i,0);
      Double_t dx = (*this)(i,1);
      cerr << " x[" << i << "] = " << setw(8) << setprecision(5) << x
           << "    "
           << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
           << setprecision(7)
           << resetiosflags(ios::showpoint)
           << endl;
   }
   cerr << "-------------------------------------------------------" << endl;
}
