#ifndef EXBPCONEHIT_H
#define EXBPCONEHIT_H

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXBPConeMeasLayer.h"

class EXBPConeHit : public TVTrackHit {
public:
   EXBPConeHit(Int_t m = kMdim);

   EXBPConeHit(const EXBPConeMeasLayer &ms,
                     Double_t          *x,
                     Double_t          *dx,
               const TVector3          &xx,
                     Double_t           b,
                     Int_t              m = kMdim);

   virtual ~EXBPConeHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;
        
   inline  const TVector3 GetExactX() const { return fXX;     }

private:
   TVector3 fXX;   // exact hit position

   ClassDef(EXBPConeHit,1)      // Sample hit class
};

#endif
