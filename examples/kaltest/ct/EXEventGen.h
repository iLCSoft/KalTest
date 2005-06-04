#ifndef __EXEVENTGEN__
#define __EXEVENTGEN__

#include "EXToyGLD.h"
#include "THelicalTrack.h"

class EXEventGen {
public:
   EXEventGen(EXToyGLD &cradle) : fCradlePtr(&cradle) { }
   virtual ~EXEventGen() { }

   THelicalTrack GenerateHelix(Double_t pt = 1.);
   void          Swim(THelicalTrack &heltrk);

private:
   EXToyGLD *fCradlePtr;

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
