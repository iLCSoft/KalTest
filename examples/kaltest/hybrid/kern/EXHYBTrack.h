#ifndef __EXHYBTRACK__
#define __EXHYBTRACK__
//*************************************************************************
//* =================
//*  EXHYBTrack Class
//* =================
//*
//* (Description)
//*   Hybrid track class for Kalman filter
//* (Requires)
//*     TKalTrack
//* (Provides)
//*     class EXHYBTrack
//* (Update Recored)
//*   2005/08/26  K.Fujii       Original version.
//*
//*************************************************************************
                                                                                
#include "TKalTrack.h"         // from KalTrackLib
#include "TAttDrawable.h"      // from Utils

//_________________________________________________________________________
//  ------------------------------
//   EXHYBTrack: Kalman Track class
//  ------------------------------
                                                                                
class EXHYBTrack : public TKalTrack, public TAttDrawable {
public:
   EXHYBTrack(Int_t n = 1) : TKalTrack(n) {}
   ~EXHYBTrack() {}

   using TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt);

   ClassDef(EXHYBTrack,1)  // Hybrid track class for Kalman Filter
};

#endif
