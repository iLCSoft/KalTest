#ifndef __EXMEASLAYER__
#define __EXMEASLAYER__
//*************************************************************************
//* ===================
//*  EXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXHit.
//* (Requires)
//* (Provides)
//*     class EXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include <iostream>
#include "TVector3.h"

#include "TKalMatrix.h"
#include "THype.h"
#include "TVMeasLayer.h"
#include "KalTrackDim.h"

using namespace std;

class TVTrackHit;

class EXMeasLayer : public TVMeasLayer, public THype {
public:

   // Ctors and Dtor

   EXMeasLayer(Double_t r0, Double_t lhalf, TVector3 we, TVector3 wd,
               Double_t cellw, Double_t fwdx0inv, Double_t bwdx0inv)
         : TVMeasLayer(fwdx0inv, bwdx0inv),
           THype(r0,lhalf,wd.Perp()/wd.Z()), fWireEnd(we), fWireDir(wd),
           fCellWidth(cellw), 
           fWirePhi(TMath::ATan2(we.Y(),we.X())), 
           fDphi(cellw/we.Perp()) 
   {
   }
   virtual ~EXMeasLayer() {}

   // Getters

   inline       TVector3     GetWireEnd  () const { return fWireEnd;     }
   inline       TVector3     GetWireDir  () const { return fWireDir;     }
   inline       Double_t     GetCellWidth() const { return fCellWidth;   }

   inline       TVector3     GetWireEnd(Int_t cellno) const 
   {
      TVector3 we = fWireEnd;
      we.RotateZ(cellno * fDphi);
      return we;
   }
   inline       TVector3     GetWireDir(Int_t cellno) const 
   {
      TVector3 wd = fWireDir;
      wd.RotateZ(cellno * fDphi);
      return wd;
   }

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv)     const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht)     const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)      const;

   // Utiliy Methods

   virtual Int_t      CalcCellNo(const TVector3   &xv)     const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv,
                                       Int_t      &lr,
                                       Int_t      &cellno) const;
   virtual TVector3   MvToXv    (const TKalMatrix &mv,
                                       Int_t       lr,
                                       Int_t       cellno) const;

private:
   TVector3    fWireEnd;	// wire end
   TVector3    fWireDir;	// wire direction
   Double_t    fCellWidth;      // cell width at the wire end
   Double_t    fWirePhi;	// wire phi   at the wire end
   Double_t    fDphi;		// cell width at the wire end in radian

   ClassDef(EXMeasLayer,1) 	// Sample measurement layer class
};

#endif
