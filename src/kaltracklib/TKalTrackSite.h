#ifndef TKALTRACKSITE_H
#define TKALTRACKSITE_H
//*************************************************************************
//* =====================
//*  TKalTrackSite Class
//* =====================
//*
//* (Description)
//*   Track measurement site class used by Kalman filter.
//* (Requires)
//*     TKalTrackState
//* (Provides)
//*     class TKalTrackSite
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2004/09/17  K.Fujii           Added ownership flag.
//*   2010/04/06  K.Fujii           Added a setter for the pivot and a
//*                                 condition object
//*
//*************************************************************************

#include "TVector3.h"    // from ROOT
#include "TVKalSite.h"   // from KalLib
#include "TVTrackHit.h"  // from KalTrackLib

#include "TTrackFrame.h" // from KalTrackLib
#include "TBField.h"     // from Bfield

class TVKalState;
class TKalTrackState;
class TKalFilterCond;

//_________________________________________________________________________
//  ---------------------------------
//  Class for Kalman measurement site
//  ---------------------------------
//
class TKalTrackSite : public TVKalSite {
public:
   TKalTrackSite(const TKalTrackSite&) = default ;
   TKalTrackSite& operator=(const TKalTrackSite&) = default ;
   TKalTrackSite(Int_t m = kMdim, Int_t p = kSdim);
   TKalTrackSite( const TVTrackHit &ht, 
                        Int_t       p = kSdim);
   ~TKalTrackSite();

   Int_t        CalcExpectedMeasVec  (const TVKalState &a, TKalMatrix &h);
   Int_t        CalcMeasVecDerivative(const TVKalState &a, TKalMatrix &H);
   Bool_t       IsAccepted();

   void         DebugPrint() const;

   inline       TTrackFrame GetFrame() const { return fFrame; }

   inline  void SetFrame(TTrackFrame frame) { fFrame = frame; }

   inline const TVTrackHit & GetHit    () const { return *fHitPtr;             }
   inline       Bool_t       IsInB     () const { return GetBfield() != 0.;    }
   inline       Bool_t       IsHitOwner() const { return fIsHitOwner;          }

   inline       void   SetPivot     (const TVector3 &x0)  { fX0 = x0;          }
   inline       void   SetHitOwner  (Bool_t b=kTRUE)      { fIsHitOwner = b;   }
   inline       void   SetFilterCond(TKalFilterCond *cp)  { fCondPtr    = cp;  }

   inline       TVector3     GetPivot      () const { return fX0;     }
                TVector3     GetLocalPivot () const;

   inline       TVector3     GetGlobalPivot() const { return fX0;     }
   inline       void         SetBfield (Double_t b) { fBfield = b;    }

   Double_t     GetBfield () const 
   { 
	   if(!TBField::IsUsingUniformBfield()) {
		   return fBfield;
	   } 
	   else {
		   return fHitPtr->GetBfield();
	   }
   }

private:
   TVKalState & CreateState(const TKalMatrix &sv, Int_t type = 0);
   TVKalState & CreateState(const TKalMatrix &sv,
                            const TKalMatrix &C,
                                  Int_t       type = 0);
   Int_t        CalcXexp   (const TVKalState &a, 
                                  TVector3   &xx,
                                  Double_t   &phi) const;

private:
   const TVTrackHit     *fHitPtr{};     // pointer to corresponding hit
         TVector3        fX0{};         // global pivot
         Bool_t          fIsHitOwner{}; // true if site owns hit
         TKalFilterCond *fCondPtr{};    // pointer to filter condition object

         TTrackFrame  fFrame{};         // site specific local frame

         Double_t        fBfield{};

   ClassDef(TKalTrackSite,1)  // sample measurement site class
};
                                                                                
#endif
