#ifndef __EXITFBHIT__
#define __EXITHFBIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXITFBMeasLayer.h"

class EXITFBHit : public TVTrackHit {
public:
   EXITFBHit(Int_t m = kMdim);

   EXITFBHit(const EXITFBMeasLayer &ms,
                   Double_t        *x,
                   Double_t        *dx,
             const TVector3        &xx,
                   Double_t         b,
                   Int_t            m = kMdim);

   virtual ~EXITFBHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;

private:

   ClassDef(EXITFBHit,1)      // Sample hit class
};

#endif
