#ifndef __EXITHIT__
#define __EXITHIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXITMeasLayer.h"

class EXITHit : public TVTrackHit {
public:
   EXITHit(Int_t m = kMdim);

   EXITHit(const EXITMeasLayer &ms,
                 Double_t      *x,
                 Double_t      *dx,
           const TVector3      &xx,
                 Double_t       b,
                 Int_t          m = kMdim);

   virtual ~EXITHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;

private:

   ClassDef(EXITHit,1)      // Sample hit class
};

#endif
