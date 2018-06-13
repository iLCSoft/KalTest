#ifndef TKALDETCRADLE_H
#define TKALDETCRADLE_H
//*************************************************************************
//* =====================
//*  TKalDetCradle Class
//* =====================
//*
//* (Description)
//*   A sigleton to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObjArray
//* 	TVKalDetector
//* (Provides)
//* 	class TKalDetCradle
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi	 Original Version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*   2010/04/06  K.Fujii        Modified Transport() to allow a 1-dim hit,
//*                              for which pivot is at the xpected hit.
//*
//*************************************************************************

#include "TObjArray.h"     // from ROOT
#include "TAttElement.h"   // from Utils
#include "TKalMatrix.h"    // from KalTrackLib
#include "TKalTrack.h"     // from KalTrackLib
#include <memory>          // from STL

class TKalTrackSite;
class TVKalDetector;
class TVMeasLayer;
class TVTrack;
class TMaterial;

//_____________________________________________________________________
//  ------------------------------
//   Detector system class
//  ------------------------------

class TKalDetCradle : public TObjArray, public TAttElement {
public:
   TKalDetCradle(Int_t n = 1);
   virtual ~TKalDetCradle();

   // Utility methods
   virtual void Install(TVKalDetector &det);

   inline virtual void   SwitchOnMS   ()       { fIsMSON = kTRUE;    }
   inline virtual void   SwitchOffMS  ()       { fIsMSON = kFALSE;   }
   inline virtual void   SwitchOnDEDX ()       { fIsDEDXON = kTRUE;  }
   inline virtual void   SwitchOffDEDX()       { fIsDEDXON = kFALSE; }
   inline virtual void   Close        ()       { fIsClosed = kTRUE;
                                                 Update();           }
   inline virtual void   Reopen       ()       { fIsClosed = kFALSE; }
   inline virtual Bool_t IsMSOn       () const { return fIsMSON;     }
   inline virtual Bool_t IsDEDXOn     () const { return fIsDEDXON;   }
   inline virtual Bool_t IsClosed     () const { return fIsClosed;   }

   void Transport(const TKalTrackSite  &from, // site from
                        TKalTrackSite  &to,   // site to
                        TKalMatrix     &sv,   // state vector
                        TKalMatrix     &F,    // propagator matrix
                        TKalMatrix     &Q);   // process noise matrix

   int  Transport(const TKalTrackSite  &from, // site from
                  const TVMeasLayer    &to,   // layer to reach
                        TVector3       &x0,   // intersection point
                        TKalMatrix     &sv,   // state vector
                        TKalMatrix     &F,    // propagator matrix
                        TKalMatrix     &Q);   // process noise matrix

   int  Transport(const TKalTrackSite  &from, // site from
                  const TVMeasLayer    &to,   // layer to reach
                        TVector3       &x0,   // intersection point
                        TKalMatrix     &sv,   // state vector
                        TKalMatrix     &F,    // propagator matrix
                        TKalMatrix     &Q,    // process noise matrix
			std::unique_ptr<TVTrack> &help);// pointer to updated track object

   int  Transport2(const TKalTrackSite  &from,   // site from
                  const TVMeasLayer     &to,     // layer to reach
                        TVector3        &x0,     // intersection point
                        TKalMatrix      &sv,     // state vector
                        TKalMatrix      &F,      // propagator matrix
                        TKalMatrix      &Q,      // process noise matrix
			std::unique_ptr<TVTrack>    &help);  // pointer to updated track object

   static void SetUseRungeKuttaTrack(Bool_t b) { fUseRKTrack = b; }

private:
   void Update();

private:
   Bool_t    fIsMSON{};         //! switch for multiple scattering
   Bool_t    fIsDEDXON{};       //! switch for energy loss
   Bool_t    fDone{};           //! flag to tell if sorting done
   Bool_t    fIsClosed{};       //! flag to tell if cradle closed

   static Bool_t   fUseRKTrack;

   ClassDef(TKalDetCradle,1)  // Base class for detector system
};

#endif
