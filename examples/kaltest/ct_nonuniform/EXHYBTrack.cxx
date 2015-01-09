//*************************************************************************
//* =================
//*  EXHYBTrack Class
//* =================
//*
//* (Description)
//*   Track class for Kalman filter
//* (Requires)
//*     TVKalSystem
//* (Provides)
//*     class EXHYBTrack
//* (Update Recored)
//*   2003/09/30  Y.Nakashima   Original version.
//*   2005/02/23  A.Yamaguchi   Added a new data member, fMass.
//*   2005/08/25  K.Fujii       Added drawable attribute.
//*
//*************************************************************************
                                                                                
#include "EXHYBTrack.h"
#include "TKalTrackSite.h"     // from KalTrackLib
#include "TVirtualPad.h"       // from ROOT
#include "TPolyMarker3D.h"     // from ROOT
#define __DEBUG__
#ifdef __DEBUG__
#include <iostream>
#include "TBField.h"
#endif
//_________________________________________________________________________
//  ------------------------------
//   EXHYBTrack: Kalman rack class
//  ------------------------------

ClassImp(EXHYBTrack)
                                                                                
//_________________________________________________________________________
//  ----------------------------------
//   Utility Method
//  ----------------------------------
//_________________________________________________________________________
// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXHYBTrack::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad || !GetEntries()) return;
   gPad->cd();

   TPolyMarker3D *pm3dp = new TPolyMarker3D(this->GetEntries());
   pm3dp->SetBit(kCanDelete);
   pm3dp->SetMarkerColor(color);
   pm3dp->SetMarkerStyle(6);

   Int_t nhits = 0;
   TIter next(this);
   TKalTrackSite *sitep = 0;
   while ((sitep = static_cast<TKalTrackSite *>(next()))) { 
      TVector3 pos = sitep->GetGlobalPivot();
      pm3dp->SetPoint(nhits, pos.X(), pos.Y(), pos.Z());
#ifdef __DEBUG__
      TVector3 bv = TBField::GetGlobalBfield(pos);
	  std::cerr << "hit " << nhits << ": r = " << pos.Perp()
                << " B(" << pos.X()
		        << ", " << pos.Y()
		        << ", " << pos.Z()
                << ") = (" << bv.X() << ", "
                << bv.Y() << ", "
                << bv.Z() << ")"
                << std::endl;
      std::cerr << "a=(" << sitep->GetCurState()[0][0] << ", "
                << sitep->GetCurState()[1][0] << ", "
                << sitep->GetCurState()[2][0] << ", "
                << sitep->GetCurState()[3][0] << ", "
                << sitep->GetCurState()[4][0] << ", "
                << ")" << std::endl;
      std::cerr << "(ndf,chi2)=("
                << GetNDF()  << ", "
                << GetChi2() << ")"
                << std::endl << std::endl;
#endif
      nhits++;
   }
   pm3dp->Draw();
   gPad->Update();
}
