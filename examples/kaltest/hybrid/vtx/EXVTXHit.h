#ifndef __EXVTXHIT__
#define __EXVTXHIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXVTXMeasLayer.h"

class EXVTXHit : public TVTrackHit {
public:
   EXVTXHit(Int_t m = kMdim);

   EXVTXHit(const EXVTXMeasLayer &ms,
                  Double_t       *x,
                  Double_t       *dx,
            const TVector3       &xx,
                  Double_t        b,
                  Int_t           m = kMdim);

   virtual ~EXVTXHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;
        
   inline  const TVector3 GetExactX() const { return fXX;     }

private:
   TVector3 fXX;   // exact hit position

   ClassDef(EXVTXHit,1)      // Sample hit class
};

#endif
