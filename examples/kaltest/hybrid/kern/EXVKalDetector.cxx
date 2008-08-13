#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TVKalDetector.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"

Double_t EXVKalDetector::fgBfield  = 30.;
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
      new TRotMatrix("rotm","rotm", 10.,80.,10.,80.,10.,80.);
      new TTUBE("Det","Det","void",210.,210.,260.);
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
   TAttDrawable *objp;
   while ((objp = dynamic_cast<TAttDrawable *>(next()))) {
      if (objp) objp->Draw(color, opt); 
   }
   nodep->Draw("pad same");
}
