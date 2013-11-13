#ifndef EXEVENTGEN_H
#define EXEVENTGEN_H

#include "TKalDetCradle.h"
#include "THelicalTrack.h"
#include "TBField.h"

class EXEventGen {
public:
   EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt,
                               Double_t cosmin,
                               Double_t cosmax);
   void          Swim(THelicalTrack &heltrk,
                      Double_t       mass = 0.13957018);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
