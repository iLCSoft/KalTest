#ifndef EXBPKALDETECTOR_H
#define EXBPKALDETECTOR_H

#include "EXVKalDetector.h"

class EXBPKalDetector : public EXVKalDetector {
public:
   EXBPKalDetector(Int_t m = 100);
   ~EXBPKalDetector();

   ClassDef(EXBPKalDetector,1)   // Sample hit class
};

#endif
