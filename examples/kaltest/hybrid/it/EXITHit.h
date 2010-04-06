#ifndef __EXITHIT__
#define __EXITHIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXITMeasLayer.h"

class EXITHit : public TVTrackHit {
public:
#ifdef TWO_DIM
   EXITHit(Int_t m = kMdim);

   EXITHit(const EXITMeasLayer &ms,
                 Double_t      *x,
                 Double_t      *dx,
           const TVector3      &xx,
                 Double_t       b,
                 Int_t          m = kMdim);
#else
   EXITHit(Int_t m = 1);

   EXITHit(const EXITMeasLayer &ms,
                 Double_t      *x,
                 Double_t      *dx,
           const TVector3      &xx,
                 Double_t       b,
                 Int_t          m = 1);
#endif

   virtual ~EXITHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;

   inline  const TVector3 GetExactX() const { return fXX;     }

private:
   TVector3 fXX;        // exact hit position

   ClassDef(EXITHit,1)      // Sample hit class
};

#endif
