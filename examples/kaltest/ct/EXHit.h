#ifndef __EXHIT__
#define __EXHIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXMeasLayer.h"

class EXHit : public TVTrackHit {
public:
   EXHit(Int_t m = kMdim);

   EXHit(const EXMeasLayer &ms,
               Double_t    *x,
               Double_t    *dx,
         const TVector3    &xx,
               Double_t     b,
               Int_t        m = kMdim);

   virtual ~EXHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;

private:

   ClassDef(EXHit,1)      // Sample hit class
};

#endif
