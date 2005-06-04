#include <iostream>
#include "TMatrixD.h"

using namespace std;

class EXHit : public TMatrixD {
public:
   EXHit(Int_t m = 1) : TMatrixD(m,2), fDim(m), fT(0.) {}
   EXHit(Double_t t, Double_t* x, Double_t *dx, Int_t m = 1) 
        : TMatrixD(m,2), fDim(m), fT(t)
   {
      for (Int_t i=0; i<m; i++) {
         (*this)(i,0) = x [i];
         (*this)(i,1) = dx[i];
      } 
   }

   ~EXHit() {}

   Double_t GetX (Int_t i) const { return (*this)(i,0); } 
   Double_t GetDX(Int_t i) const { return (*this)(i,1); }
   Double_t GetT ()        const { return fT; } 

   void DebugPrint() const
   {
      for (Int_t i=0; i<fDim; i++) {
         Double_t x  = (*this)(i,0);
         Double_t dx = (*this)(i,1);
         cerr << "t = " << fT << " x = " << x << " dx = " << dx << endl;
      } 
   }

private:
   Int_t    fDim;         // dimension of coordinate space
   Double_t fT;           // independent variable

   ClassDef(EXHit,1)      // Sample hit class
};
