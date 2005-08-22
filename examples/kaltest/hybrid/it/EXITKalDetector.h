#ifndef __EXITDETECTOR__
#define __EXITDETECTOR__

#include "EXVKalDetector.h"

class EXITKalDetector : public EXVKalDetector {
public:
   EXITKalDetector(Int_t m = 100);
   ~EXITKalDetector();

   ClassDef(EXITKalDetector,1)   // Sample hit class
};

#endif
