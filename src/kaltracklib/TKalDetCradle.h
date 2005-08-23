#ifndef __TKALDETCRADLE__
#define __TKALDETCRADLE__
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
//*   2005/02/23  A.Yamaguchi	   Original Version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*
//*************************************************************************

#include "TObjArray.h"     // from ROOT
#include "TAttElement.h"   // from Utils
#include "TKalMatrix.h"    // from KalTrackLib

class TKalTrackSite;
class TVKalDetector;

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
   inline virtual Bool_t IsMSOn       () const { return fIsMSON;     }
   inline virtual Bool_t IsDEDXOn     () const { return fIsDEDXON;   }

   void Transport(const TKalTrackSite  &from, // site from
                  const TKalTrackSite  &to,   // sit to
                        TKalMatrix     &sv,   // state vector
                        TKalMatrix     &F,    // propagator matrix
                        TKalMatrix     &Q);   // process noise matrix
private:
   void Update();

private:
   Bool_t    fIsMSON;         //! switch for multiple scattering
   Bool_t    fIsDEDXON;       //! switch for energy loss
   Bool_t    fDone;           //! flag to tell if sorting done

   ClassDef(TKalDetCradle,1)  // Base class for detector system
};

#endif
