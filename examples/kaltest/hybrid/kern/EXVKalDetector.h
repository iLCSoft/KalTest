#ifndef __EXVDETECTOR__
#define __EXVDETECTOR__

#include "TVector3.h"
#include "TVKalDetector.h"

class TVMeasLayer;

class EXVKalDetector : public TVKalDetector {
public:
   EXVKalDetector(Int_t m = 100);
   ~EXVKalDetector();

   virtual void ProcessHit(const TVector3    &xx,
                           const TVMeasLayer &ms, 
                                 TObjArray   &hits) = 0;

   Double_t     GetBfield (const TVector3 &xx = TVector3(0.)) const
                           { return fgBfield; }

private:
   static Double_t fgBfield;   // magnetic field [kG]

   ClassDef(EXVKalDetector,1)   // Sample hit class
};

#endif
