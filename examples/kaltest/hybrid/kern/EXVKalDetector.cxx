#include "EXVKalDetector.h"
#include "TVKalDetector.h"

Double_t EXVKalDetector::fgBfield = 30.;

ClassImp(EXVKalDetector)

EXVKalDetector::EXVKalDetector(Int_t m)
             : TVKalDetector(m),
               fIsPowerOn(kTRUE)
{
}

EXVKalDetector::~EXVKalDetector()
{
}


