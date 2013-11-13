#ifndef __EXRUNGEKUTTAEVENTGEN__
#define __EXRUNGEKUTTAEVENTGEN__

#include "TKalDetCradle.h"
#include "TRungeKuttaTrack.h"

class EXRungeKuttaEventGen {
public:
   EXRungeKuttaEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXRungeKuttaEventGen() {}

   TRungeKuttaTrack GenerateRKTrack(Double_t chg, TVector3 x, TVector3 p);

   void          Swim(TRungeKuttaTrack &rktrk,
   Double_t      mass = 0.13957018);

   void          DumpHits();

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;   // pointer to detector system
   TObjArray     *fHitBufPtr;   // pointer to hit array

   static Double_t  fgT0;       // t0

   ClassDef(EXRungeKuttaEventGen,1)   // Event Generator
};

#endif
