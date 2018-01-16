#ifndef TRKTRACK_H
#define TRKTRACK_H
//*************************************************************************
//* ====================
//*  TRKTrack Class
//* ====================
//*
//* (Description)
//*   A class to implement a Runge-Kutta track object.
//*
//*
//* (Requires)
//*     TVTrack, TEveTrackPropagator
//* (Provides)
//*     class TRKTrack
//*
//*************************************************************************
//
#include "TParticle.h"
#include "TEveVSDStructs.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"

#include "TVector3.h"

#include "TVMeasLayer.h"

//_____________________________________________________________________
//  -----------------------------------
//  Runge-Kutta Track Class
//  -----------------------------------
class TRKTrack : public TEveTrackPropagator {
public:
   // Ctors and Dtor
   TRKTrack(Double_t chg, TVector3 x, TVector3 p);
   TRKTrack(const TRKTrack&);
   TRKTrack& operator=(const TRKTrack&) = default ;

   virtual ~TRKTrack() {}

   // Utility methods
   void StepRungeKutta(Double_t step);
 
   Double_t GetCharge() const;
 
   void SetCharge(Int_t chg);

   inline void SetCurPosition(TVector3 x) { fPosition = x; }
   inline void SetCurMomentum(TVector3 p) { fMomentum = p; }

   inline TVector3 GetCurPosition() { return fPosition; }
   inline TVector3 GetCurMomentum() { return fMomentum; }

private:
   TVector3 fPosition{};
   TVector3 fMomentum{};

   TEveMagField* fMagField{};

   ClassDef(TRKTrack,1)   
};

#endif
