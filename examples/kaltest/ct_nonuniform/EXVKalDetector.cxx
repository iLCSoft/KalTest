#include "EXVKalDetector.h"
//#include "EXVMeasLayer.h"
#include "TVKalDetector.h"
//#include "TTUBE.h"
#include "TTUBS.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"

#include <iostream>

Double_t EXVKalDetector::fgBfield  = 10.;
TNode   *EXVKalDetector::fgNodePtr = 0;

ClassImp(EXVKalDetector)

EXVKalDetector::EXVKalDetector(Int_t m)
             : TVKalDetector(m),
               fIsPowerOn(kTRUE)
{
}

EXVKalDetector::~EXVKalDetector()
{
}

TNode *EXVKalDetector::GetNodePtr()
{
   if (!fgNodePtr) {
      TRotMatrix* rot = new TRotMatrix("rotm","rotm", 10.,80.,10.,80.,10.,80.);
      new TTUBE("Det","Det","void",1820.,1820.,2020.);
      fgNodePtr = new TNode("World","World","Det",0.,0.,0.,"rotm");
   }
   return fgNodePtr;
}

void EXVKalDetector::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   TNode *nodep = GetNodePtr();
   nodep->cd();
   TIter next(this);
   TObject *objp;
   while ((objp = next())) {
      TAttDrawable *dp = dynamic_cast<TAttDrawable *>(objp);
      if (dp) dp->Draw(color, opt); 
   }
   nodep->Draw("pad same");
}
