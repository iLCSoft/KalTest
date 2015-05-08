#ifndef TRUNGEKUTTATRACK_H
#define TRUNGEKUTTATRACK_H
//*************************************************************************
//* ====================
//*  TRungeKuttaTrack Class
//* ====================
//*
//* (Description)
//*   A class to implement a Runge-Kutta track object.
//*
//*
//* (Requires)
//* (Provides)
//*     class TRungeKuttaTrack
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TVMeasLayer.h"

//_____________________________________________________________________
//  -----------------------------------
//  Runge-Kutta Track Class
//  -----------------------------------
class TRungeKuttaTrack {
public:
   // Ctors and Dtor
   TRungeKuttaTrack(Double_t q, TVector3 x, TVector3 p);
   TRungeKuttaTrack(const TRungeKuttaTrack&);
   virtual ~TRungeKuttaTrack() {}

   //Step function
   void StepRungeKutta(Double_t step);
   //void StepRungeKutta(Double_t step,
   //                    Double_t* vect, 
   //			    	 Double_t* vout);
   void StepRungeKutta(Double_t  step, 
                       TVector3& vx, 
                       TVector3& vp);

 
   //Setter
   inline void SetCurPosition(TVector3 x) { fPosition = x; }
   inline void SetCurMomentum(TVector3 p) { fMomentum = p; }
   inline void SetCharge     (Double_t q) { fCharge   = q; }

   //Getter
   inline  TVector3 GetCurPosition() const { return fPosition; }
   inline  TVector3 GetCurMomentum() const { return fMomentum; }
   inline  Double_t GetCharge()      const { return fCharge;   }

private:
   TVector3 fPosition;
   TVector3 fMomentum;
   Double_t fCharge;

   ClassDef(TRungeKuttaTrack,1)   
};

#endif
