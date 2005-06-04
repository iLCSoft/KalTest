#ifndef __EXDETECTOR__
#define __EXDETECTOR__

#include <iostream>

#include "TVector3.h"

#include "TVKalDetector.h"
#include "EXMeasLayer.h"

using namespace std;

class EXKalDetector : public TVKalDetector {
public:
   EXKalDetector(Int_t m = 100) : TVKalDetector(m) {}
   EXKalDetector(Int_t    nlayers,
                 Int_t    nwires,
                 Double_t lhalf,
                 Double_t celhw,
                 Double_t tna,
                 Double_t rmin,
                 Double_t rstep,
                 Double_t rt0det = -999.);

   ~EXKalDetector() {}

   Double_t CalcRadLength(const TVMeasLayer &from,
                          const TVMeasLayer &to) const;
   Double_t CalcSigmaMS0 (const TVMeasLayer &from,
                          const TVMeasLayer &to,
                                Double_t     pathlen) const;

private:
   Double_t fRT0Det;   // radius of T0 detector

   ClassDef(EXKalDetector,1)      // Sample hit class
};

#endif
