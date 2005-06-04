#ifndef __EXTOYGLD__
#define __EXTOYGLD__
//*************************************************************************
//* =====================
//*  EXToyGLD Class
//* =====================
//*
//* (Description)
//*   A sigleton to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObjArray
//* 	TVKalDetector
//* (Provides)
//* 	class EXToyGLD
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi	Original Version.
//*
//*************************************************************************

#include "TKalDetCradle.h"

//_____________________________________________________________________
//  ------------------------------
//  Detector system class
//  ------------------------------
//

class EXToyGLD : public TKalDetCradle {
public:
   EXToyGLD(Int_t n = 1);
   virtual ~EXToyGLD();

   const TObjArray &GetHits() const { return fHits;   }
         TObjArray &GetHits()       { return fHits;   }
         void       Reset  ()       { fHits.Delete(); }

private:
   TObjArray fHits;

   ClassDef(EXToyGLD,1)  // Base class for detector system
};

//=======================================================
// inline functions, if any
//=======================================================

#endif
