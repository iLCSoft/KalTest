#ifndef __EXTPCHIT__
#define __EXTPCHIT__

#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "EXTPCMeasLayer.h"

class EXTPCHit : public TVTrackHit {
public:
   EXTPCHit(Int_t m = kMdim);

   EXTPCHit(const EXTPCMeasLayer &ms,
                  Double_t       *x,
                  Double_t       *dx,
                  Int_t           side,
                  Double_t        v,
            const TVector3       &xx,
                  Double_t        b,
                  Int_t           m = kMdim);

   virtual ~EXTPCHit();

   virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;

   virtual void       DebugPrint(Option_t *opt = "")           const;

   Int_t    GetSide  () const { return fSide;   }
   Double_t GetVdrift() const { return fVdrift; }

private:
   Int_t    fSide;
   Double_t fVdrift;

   ClassDef(EXTPCHit,1)      // Sample hit class
};

#endif
