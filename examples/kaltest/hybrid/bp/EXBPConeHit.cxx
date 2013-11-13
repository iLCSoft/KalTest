#include "EXBPConeHit.h"
#include "EXBPConeMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXBPConeHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXBPConeHit::EXBPConeHit(Int_t m)
           : TVTrackHit(m)
{
}
      
EXBPConeHit::EXBPConeHit(const EXBPConeMeasLayer &ms,
                               Double_t          *x,
                               Double_t          *dx, 
                         const TVector3          &xx,
                               Double_t           b,
                               Int_t              m)
           : TVTrackHit(ms, x, dx, b, m),
             fXX(xx)
{
}

EXBPConeHit::~EXBPConeHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXBPConeHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   const EXBPConeMeasLayer &ms
                 = dynamic_cast<const EXBPConeMeasLayer &>(GetMeasLayer());
   TKalMatrix h  = ms.XvToMv(xv);
   TVector3 xxv  = xv - ms.GetXc();
   Double_t tana = ms.GetTanA();
   Double_t r    = xxv.Z() * tana;
   Double_t rm   = ((*this)(1,0)-ms.GetXc().Z()) * tana;

   Double_t phih = h(0,0) / r;
   Double_t phim = (*this)(0,0) / rm;
   Double_t dphi = phih - phim;

   static const Double_t kPi    = TMath::Pi();
   static const Double_t kTwoPi = 2 * kPi;

   while (dphi < -kPi) dphi += kTwoPi;
   while (dphi >  kPi) dphi -= kTwoPi;

   h(0,0)  = r * (phim + dphi);
   h(1,0) += 0.;

   return h;
}

void EXBPConeHit::DebugPrint(Option_t *) const
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
