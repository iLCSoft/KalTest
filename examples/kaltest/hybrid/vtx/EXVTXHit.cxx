#include "EXVTXHit.h"
#include "EXVTXMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXVTXHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXVTXHit::EXVTXHit(Int_t m)
        : TVTrackHit(m)
{
}
      
EXVTXHit::EXVTXHit(const EXVTXMeasLayer &ms,
                         Double_t       *x,
                         Double_t       *dx, 
                   const TVector3       &xx,
                         Double_t        b,
                         Int_t           m)
        : TVTrackHit(ms, x, dx, b, m),
          fXX(xx)
{
}

EXVTXHit::~EXVTXHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXVTXHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   return dynamic_cast<const EXVTXMeasLayer &>(GetMeasLayer()).XvToMv(xv);
}

void EXVTXHit::DebugPrint(Option_t *) const
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
