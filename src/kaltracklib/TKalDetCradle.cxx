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
//*   2005/02/23  A.Yamaguchi  	Original version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*
//*************************************************************************

#include "TKalDetCradle.h"   // from KalTrackLib
#include "TVMeasLayer.h"     // from KalTrackLib
#include "TVKalDetector.h"   // from KalTrackLib
#include "TKalTrackSite.h"   // from KalTrackLib
#include "TKalTrackState.h"  // from KalTrackLib
#include "TVSurface.h"       // from GeomLib
#include <memory>            // from STL

ClassImp(TKalDetCradle)

//_________________________________________________________________________
//  ----------------------------------
//   Ctors and Dtor
//  ----------------------------------

TKalDetCradle::TKalDetCradle(Int_t n)
             : TObjArray(n), fIsMSON(kTRUE), fIsDEDXON(kTRUE), fDone(kFALSE)
{
}

TKalDetCradle::~TKalDetCradle()
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________
// -----------------
//  Install
// -----------------
//    installs a sub-detector into this cradle.
//
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

//_________________________________________________________________________
// -----------------
//  Transport
// -----------------
//    transports state (sv) from site (from) to site (to), taking into 
//    account multiple scattering and energy loss and updates state (sv),
//    propagator matrix (F), and process noise matrix (Q).
//
void TKalDetCradle::Transport(const TKalTrackSite  &from,  // site from
                              const TKalTrackSite  &to,    // site to
                                    TKalMatrix     &sv,    // state vector
                                    TKalMatrix     &F,     // propagator matrix
                                    TKalMatrix     &Q)     // process noise matrix
{
   // ---------------------------------------------------------------------
   //  Sort measurement layers in this cradle if not
   // ---------------------------------------------------------------------
   if (!fDone) Update();
	
   // ---------------------------------------------------------------------
   //  Locate sites from and to in this cradle
   // ---------------------------------------------------------------------
   Int_t  fridx = from.GetHit().GetMeasLayer().GetIndex(); // index of site from
   Int_t  toidx = to.GetHit().GetMeasLayer().GetIndex();   // index of site to
   Int_t  di    = fridx > toidx ? -1 : 1;                  // layer increment
   Bool_t isout = di > 0 ? kTRUE : kFALSE;  // out-going or in-coming

   std::auto_ptr<TVTrack> help(&static_cast<TKalTrackState &>
                              (from.GetCurState()).CreateTrack()); // tmp track
   TVTrack &hel = *help;
	
   TVector3 xx;                // expected hit position vector
   Double_t fid     = 0.;      // deflection angle from the last hit

   Int_t sdim = sv.GetNrows();                // # track parameters
   for (Int_t p=0; p<sdim; p++) F(p, p) = 1.; // initialize F to unity
   TKalMatrix DF(sdim, sdim);                 // propagator matrix segment

   // ---------------------------------------------------------------------
   //  Loop over layers and transport sv, F, and Q step by step
   // ---------------------------------------------------------------------
   Int_t ifr = fridx;
   Int_t ito = fridx + di;
   while ((ifr-fridx)*(ifr-toidx)<=0 && (ito-fridx)*(ito-toidx)<=0) {
      if(dynamic_cast<TVSurface *>(At(ito))->CalcXingPointWith(hel, xx, fid)) {
         const TVMeasLayer   &ml  = *dynamic_cast<TVMeasLayer *>(At(ifr));
         TKalMatrix Qms(sdim, sdim);
         if (IsMSOn()) ml.CalcQms(isout, hel, fid, Qms); // Qms for this step

         hel.MoveTo(xx, fid, &DF);         // move to next expected hit
         if (sdim == 6) DF(5, 5) = 1.;     // t0 stays the same
         F = DF * F;                       // update F
         TKalMatrix DFt  = TKalMatrix(TMatrixD::kTransposed, DF);

         Q = DF * (Q + Qms) * DFt;         // transport Q to next expected hit 

         if (IsDEDXOn()) {
            hel.PutInto(sv);                              // copy hel to sv
            sv(2,0) += ml.GetEnergyLoss(isout, hel, fid); // correct for dE/dx
            hel.SetTo(sv, hel.GetPivot());                // save sv back to hel
         }
         ifr = ito;
      }
      ito += di;
   }
   // ---------------------------------------------------------------------
   //  Move pivot from last expected hit to actural hit at site to
   // ---------------------------------------------------------------------
   hel.MoveTo(to.GetPivot(), fid, &DF); // move pivot to actual hit (to)
   hel.PutInto(sv);                     // save updated hel to sv
   F = DF * F;                          // update F accordingly
}

//_________________________________________________________________________
// -----------------
//  Update
// -----------------
//    sorts meaurement layers according to layer's sorting policy
//    and puts index to layers from inside to outside.
//
void TKalDetCradle::Update()
{
   fDone = kTRUE;

   UnSort();   // unsort
   Sort();     // sort layers according to sorting policy

   TIter next(this);
   TVMeasLayer *mlp = 0;
   Int_t i = 0;
   while ((mlp = dynamic_cast<TVMeasLayer *>(next()))) {
      mlp->SetIndex(i++);
   }
}
