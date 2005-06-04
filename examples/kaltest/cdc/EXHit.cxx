#include "EXHit.h"
#include "EXMeasLayer.h"

ClassImp(EXHit)

//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

EXHit::EXHit(Int_t m)
        : TVTrackHit(m), fLR(-1), fCellNo(0), fVdrift(0.7e-3)
{
}
      
EXHit::EXHit(const EXMeasLayer &ms,
                   Double_t    *x,
                   Double_t    *dx, 
                   Int_t        lr,
                   Int_t        cn,
                   Double_t     v,
             const TVector3    &xx,
                   Double_t     b,
                   Int_t        m)
        : TVTrackHit(ms,x,dx,xx,b,m), fLR(lr), fCellNo(cn), fVdrift(v)
{
}

EXHit::~EXHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix EXHit::XvToMv(const TVector3 &xv, Double_t t0) const
{
   TKalMatrix h = GetMeasLayer().XvToMv(*this,xv);

#ifndef __TPC__
   h(0,0) += fVdrift * t0 * fLR; // T0 shift
#else
   h(1,0) += fVdrift * t0;       // T0 shift
#endif
   return h;
}

void EXHit::DebugPrint(Option_t *) const
{
   cerr << "------------------- Site Info -------------------------" << endl
        << " cel  = " << setw(4) << fCellNo  << "  "
        << setiosflags(ios::showpos)
        << " LR ="    << setw(3) << fLR      << "  "
        << resetiosflags(ios::showpos)
        << setiosflags(ios::showpoint)
        << " xw = ("  << setw(7) << setprecision(5) << GetWireEnd().X()
        << ","        << setw(7) << setprecision(5) << GetWireEnd().Y()
        << ","        << setw(7) << setprecision(5) << GetWireEnd().Z()
        << ")"        << endl;

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
