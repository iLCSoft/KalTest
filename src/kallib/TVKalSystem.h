#ifndef __TVKALSYSTEM__
#define __TVKALSYSTEM__
//*************************************************************************
//* ===================
//*  TVKalSystem Class
//* ===================
//*
//* (Description)
//*   Base class for Kalman filtering class
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalSystem
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************

#include "TObjArray.h"
#include "TVKalSite.h"

//_____________________________________________________________________
//  ------------------------------
//  Kalman Filtering class
//  ------------------------------
//
class TKalMatrix;

class TVKalSystem : public TObjArray {
public:

   // Ctors and Dtor

   TVKalSystem(Int_t n = 1);
   virtual ~TVKalSystem();

   // Utility methods

   virtual Bool_t AddAndFilter(TVKalSite &next);
   virtual void   SmoothBackTo(Int_t k);
   virtual void   SmoothAll();

   inline  void Add(TObject *obj);

   // Getters

   inline virtual TVKalSite  & GetCurSite()     { return *fCurSitePtr; }
   inline virtual TVKalState & GetState(TVKalSite::EStType t) 
                                   { return fCurSitePtr->GetState(t); }
   inline virtual Double_t     GetChi2() { return fChi2; }
          virtual Int_t        GetNDF (Bool_t self = kTRUE);
   
   // Setters

private:
   TVKalSite   *fCurSitePtr;  // pointer to current site
   Double_t     fChi2;        // current total chi2
   
   ClassDef(TVKalSystem,1)  // Base class for Kalman Filter
};

//=======================================================
// inline functions
//=======================================================

void TVKalSystem::Add(TObject *obj)
{
   TObjArray::Add(obj); 
   fCurSitePtr = (TVKalSite *)obj;
}

#endif
