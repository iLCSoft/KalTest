#ifndef __EXDETECTOR__
#define __EXDETECTOR__

#include "TVector3.h"         // from ROOT

#include "EXVKalDetector.h"   
#include "EXMeasLayer.h"

class EXKalDetector : public EXVKalDetector {
public:
   // Ctor and Dtor
   EXKalDetector(Int_t m = 100);
   ~EXKalDetector();

   // Utility methods

   static Double_t GetBfield (const TVector3 &xx = TVector3(0.))
                               { return fgBfield; }

   void  Draw(Int_t color, const Char_t *opt = "");

private:
   static Double_t fgBfield;   // magnetic field [kG]
   TNode *fNodePtr;            // node pointer

   ClassDef(EXKalDetector,1)   // Sample hit class
};

#endif
