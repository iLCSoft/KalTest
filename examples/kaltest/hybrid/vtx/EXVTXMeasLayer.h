#ifndef __EXVTXMEASLAYER__
#define __EXVTXMEASLAYER__
//*************************************************************************
//* ===================
//*  EXVTXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXVTXHit.
//* (Requires)
//* (Provides)
//*     class EXVTXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*   2011/06/17  D.Kamai           Modified to handle ladder structure.
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>
class TVTrackHit;

class EXVTXMeasLayer : public EXVMeasLayer, public TPlane {
public:
   // Ctors and Dtor

   EXVTXMeasLayer(TMaterial &min,
                  TMaterial &mout,
                  const TVector3  &center,
                  const TVector3  &normal,
                  Double_t   SortingPolicy,
                  Double_t   xiwidth,
                  Double_t   zetawidth,
                  Double_t   xioffset,
                  Double_t   sigmaxi,
                  Double_t   sigmazeta,
                  Bool_t     type = EXVMeasLayer::kActive,
            const Char_t    *name = "FPCCDVTXML");
   virtual ~EXVTXMeasLayer();

   // Parrent's pure virtuals that must be implemented

   virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                 const TVector3   &xv) const;
   virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   virtual void       CalcDhDa  (const TVTrackHit &ht,
                                 const TVector3   &xv,
                                 const TKalMatrix &dxphiada,
                                       TKalMatrix &H)  const;
   
   
   inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
   
   virtual void       ProcessHit(const TVector3   &xx,
                                       TObjArray  &hits);

   Double_t GetSortingPolicy() const { return fSortingPolicy; }
   Double_t GetXiwidth() const { return fXiwidth; }
   Double_t GetZetawidth() const { return fZetawidth; }
   Double_t GetXioffset() const { return fXioffset; }
   Double_t GetSigmaXi() const { return fSigmaXi; }
   Double_t GetSigmaZeta() const { return fSigmaZeta; }

   using TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt);
   
 private:
   Double_t fSortingPolicy;
   Double_t fXiwidth;
   Double_t fZetawidth;
   Double_t fXioffset;
   Double_t fSigmaXi;  // sigma_x
   Double_t fSigmaZeta;  // sigma_z
   
   ClassDef(EXVTXMeasLayer,1)   // Sample measurement layer class   
   
};

#endif
