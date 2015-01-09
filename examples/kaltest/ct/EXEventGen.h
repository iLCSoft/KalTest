#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "TKalDetCradle.h"
#include "THelicalTrack.h"
#include "TStraightTrack.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle &cradle, TObjArray &hitbuf) : fCradlePtr(&cradle), fHitBufPtr(&hitbuf) { }
   virtual ~EXEventGen() { }

   THelicalTrack  GenerateHelix(Double_t pt = 1.);
   TStraightTrack GenerateStraightTrack(Double_t p = 5.);

   void          Swim(THelicalTrack  &heltrk);
   void          Swim(TStraightTrack &trk);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;   // pointer to detector system
   TObjArray     *fHitBufPtr;   // pointer to hit array

   static Double_t  fgT0;       // t0

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
