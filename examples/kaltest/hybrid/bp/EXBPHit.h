#ifndef EXBPHIT_H
#define EXBPHIT_H

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXBPMeasLayer.h"

class EXBPHit : public TVTrackHit {
public:
   EXBPHit(Int_t m = kMdim);

   EXBPHit(const EXBPMeasLayer &ms,
                 Double_t      *x,
                 Double_t      *dx,
           const TVector3      &xx,
                 Double_t       b,
                 Int_t          m = kMdim);

   virtual ~EXBPHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;
        
   inline  const TVector3 GetExactX() const { return fXX;     }

private:
   TVector3 fXX;   // exact hit position

   ClassDef(EXBPHit,1)      // Sample hit class
};

#endif
