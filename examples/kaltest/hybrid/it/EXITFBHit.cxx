#include "EXITFBHit.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXITFBHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXITFBHit::EXITFBHit(Int_t m)
         : TVTrackHit(m)
{
}
      
EXITFBHit::EXITFBHit(const EXITFBMeasLayer &ms,
                           Double_t        *x,
                           Double_t        *dx, 
                     const TVector3        &xx,
                           Double_t         b,
                           Int_t            m)
         : TVTrackHit(ms, x, dx, b, m),
           fXX(xx)
{
}

EXITFBHit::~EXITFBHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXITFBHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   return dynamic_cast<const EXITFBMeasLayer &>(GetMeasLayer()).XvToMv(xv);
}

void EXITFBHit::DebugPrint(Option_t *) const
{
   cerr << "------------------- Site Info -------------------------" << endl;
   cerr << " x  = " << setw(8) << setprecision(5) << (*this)(0,0)
        << " dx = " << setw(6) << setprecision(2) << (*this)(0,1) << endl
        << " y  = " << setw(8) << setprecision(5) << (*this)(1,0)
        << " dy = " << setw(6) << setprecision(2) << (*this)(1,1) << endl
        << " z  = " << setprecision(7) << resetiosflags(ios::showpoint)
                    << dynamic_cast<const EXITFBMeasLayer &>
                      (GetMeasLayer()).GetXc().Z()                << endl;
   cerr << "-------------------------------------------------------" << endl;
}
