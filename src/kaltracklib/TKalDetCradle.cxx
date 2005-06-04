//*************************************************************************
//* =====================
//*  TKalDetCradle Class
//* =====================
//*
//* (Description)
//*   A singleton class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//*     TObjArray
//*     TVKalDetector
//* (Provides)
//*     class TKalDetCradle
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi	Original version.
//*
//*************************************************************************

#include "TKalDetCradle.h"
#include "TVMeasLayer.h"
#include "TVKalDetector.h"
#include "TKalTrackSite.h"
#include "TVSurface.h"
#include "TObjNum.h"
#include <memory>

ClassImp(TKalDetCradle)

TKalDetCradle::TKalDetCradle(Int_t n)
             : TObjArray(n), fIsForward(kTRUE), fIsMSON(kTRUE), fDone(kFALSE)
{
}

TKalDetCradle::~TKalDetCradle()
{
}

void TKalDetCradle::Install(TVKalDetector &det)
{
   TIter next(&det);
   TObject *mlp = 0;  // measment layer pointer
   while ((mlp = next())) {
      Add(mlp);
      dynamic_cast<TAttElement *>(mlp)->SetParentPtr(&det);
      det.SetParentPtr(this);
   }
   fDone = kFALSE;
}

void TKalDetCradle::CalcTable(const TKalTrackSite &from,
                              const TKalTrackSite &to)
{
   if (!fDone) Update();
   fMeasLayerTable.Clear();
   fDPhiTable.Delete();

   Int_t fridx = from.GetHit().GetMeasLayer().GetIndex();
   Int_t toidx = to.GetHit().GetMeasLayer().GetIndex();
   Int_t di    = fridx > toidx ? -1 : 1;
   fIsForward  = di > 0 ? kTRUE : kFALSE;
   std::auto_ptr<TVTrack> hel(&static_cast<TKalTrackState &>
                             (from.GetCurState()).CreateTrack());
   TVector3 xx;
   Double_t phi     = 0.;
   Double_t philast = 0.;

   for (Int_t i=fridx; (di>0 && i<=toidx-di)||(di<0 && i>=toidx-di); i += di) {
      if(dynamic_cast<TVSurface*>(At(i+di))->CalcXingPointWith(*hel, xx, phi)) {
         fMeasLayerTable.Add(At(i));
         fDPhiTable.Add(new TObjNum(phi - philast));
         philast = phi;
      }
   }
}

void TKalDetCradle::Update()
{
   fDone = kTRUE;

   UnSort();
   Sort();

   TIter next(this);
   TVMeasLayer *mlp = 0;
   Int_t i = 0;
   while ((mlp = dynamic_cast<TVMeasLayer *>(next()))) {
      mlp->SetIndex(i++);
   }
}
