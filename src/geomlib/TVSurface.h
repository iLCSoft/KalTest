#ifndef __TVSURFACE__
#define __TVSURFACE__
//*************************************************************************
//* ====================
//*  TVSurface Class
//* ====================
//*
//* (Description)
//*   This is the base class for various solids.
//* (Requires)
//*     TObject;
//* (Provides)
//*     class TVSurface
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*   2005/02/23  K.Fujii       Added new methods, Compare() and
//*                             GetSortingPolicy().
//*
//*************************************************************************
//
#include "TObject.h"
#include "TMatrixD.h"
#include "TVector3.h"

class TVTrack;
//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any surface
//  -----------------------------------

class TVSurface : public TObject {
public:

   virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                            TVector3 &xx,
                                            Double_t &phi,
                                            Double_t  eps = 1.e-8) const;
   virtual Double_t CalcS            (const TVector3 &xx) const = 0;
   virtual TMatrixD CalcDSDx         (const TVector3 &xx) const = 0;
   virtual Bool_t   IsOnSurface      (const TVector3 &xx) const = 0;
   virtual Bool_t   IsOutside        (const TVector3 &xx) const = 0;

   virtual Double_t GetSortingPolicy ()                   const = 0;

   virtual Int_t    Compare   (const TObject *obj) const;
   virtual Bool_t   IsSortable()                   const { return kTRUE; }

private:
 
   ClassDef(TVSurface,1)      // Base class for any surface
};

#endif
