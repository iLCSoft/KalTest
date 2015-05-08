#ifndef __EXRKEVENTGEN__
#define __EXRKEVENTGEN__

#include "TKalDetCradle.h"
#include "TRKTrack.h"

class EXRKEventGen {
public:
   EXRKEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXRKEventGen() {}

   TRKTrack GenerateRKTrack(Double_t chg, TVector3 x, TVector3 p);

   void          Swim(TRKTrack &rktrk,
   Double_t      mass = 0.13957018);

   void          DumpHits();

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;   // pointer to detector system
   TObjArray     *fHitBufPtr;   // pointer to hit array

   static Double_t  fgT0;       // t0

   ClassDef(EXRKEventGen,1)   // Event Generator
};

#endif
