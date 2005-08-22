#ifndef __EXVDETECTOR__
#define __EXVDETECTOR__

#include "TVector3.h"
#include "TVKalDetector.h"

class TVMeasLayer;

class EXVKalDetector : public TVKalDetector {
public:
   EXVKalDetector(Int_t m = 100);
   ~EXVKalDetector();

   static Double_t GetBfield (const TVector3 &xx = TVector3(0.))
                             { return fgBfield; }

private:
   static Double_t fgBfield;   // magnetic field [kG]

   ClassDef(EXVKalDetector,1)   // Sample hit class
};

#endif
