#ifndef __EXITDETECTOR__
#define __EXITDETECTOR__

#include "EXVKalDetector.h"

class EXITKalDetector : public EXVKalDetector {
public:
   EXITKalDetector(Int_t m = 100);
   ~EXITKalDetector();

   void ProcessHit(const TVector3    &xx,
                   const TVMeasLayer &ms,
                         TObjArray   &hits);

   ClassDef(EXITKalDetector,1)   // Sample hit class
};

#endif
