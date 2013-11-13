#ifndef __EXVTXDETECTOR__
#define __EXVTXDETECTOR__

#include "EXVKalDetector.h"

class EXVTXKalDetector : public EXVKalDetector {
public:
   EXVTXKalDetector(Int_t m = 100);
   ~EXVTXKalDetector();

   ClassDef(EXVTXKalDetector,1)   // Sample hit class
};

#endif
