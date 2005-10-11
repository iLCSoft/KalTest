#ifndef __EXVDETECTOR__
#define __EXVDETECTOR__

#include "TVector3.h"
#include "TVKalDetector.h"

class TVMeasLayer;

class EXVKalDetector : public TVKalDetector {
public:
   EXVKalDetector(Int_t m = 100);
   virtual ~EXVKalDetector();

   inline virtual Bool_t IsPowerOn() const { return fIsPowerOn;   }
   inline virtual void   PowerOn  ()       { fIsPowerOn = kTRUE;  }
   inline virtual void   PowerOff ()       { fIsPowerOn = kFALSE; }

   static Double_t GetBfield (const TVector3 &xx = TVector3(0.))
                             { return fgBfield; }

private:
   Bool_t fIsPowerOn;
   static Double_t fgBfield;   // magnetic field [kG]

   ClassDef(EXVKalDetector,1)   // Sample hit class
};

#endif
