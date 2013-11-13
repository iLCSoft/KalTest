#ifndef __EXDETECTOR__
#define __EXDETECTOR__

#include <iostream>

#include "TVector3.h"

#include "TVKalDetector.h"
#include "EXMeasLayer.h"

using namespace std;

class EXKalDetector : public TVKalDetector {
public:
   EXKalDetector(Int_t    nlayers = 50,
                 Int_t    nwires  = 5,
                 Double_t lhalf   = 150.,
                 Double_t celhw   = 2.5,
                 Double_t tna     = 0.05,
                 Double_t rmin    = 45.,
                 Double_t rstep   = 2.2,
                 Double_t rt0det = -999.);

   ~EXKalDetector() {}

   static Double_t GetBfield (const TVector3 &xx = TVector3(0.))
                           { return fgBfield; }

private:
   Double_t fRT0Det;         // radius of T0 detector

   static Double_t fgBfield; //! magnetic field [kG]

   ClassDef(EXKalDetector,1)      // Sample hit class
};

#endif
