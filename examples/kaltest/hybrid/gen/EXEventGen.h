#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "TKalDetCradle.h"
#include "THelicalTrack.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt = 1.);
   void          Swim(THelicalTrack &heltrk);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
