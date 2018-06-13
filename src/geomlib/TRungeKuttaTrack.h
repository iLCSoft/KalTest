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
#include "THelicalTrack.h"
#include "TVSurface.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Runge-Kutta Track Class
//  -----------------------------------
class TRungeKuttaTrack : public THelicalTrack {

public:
   // Ctors and Dtor
   
   TRungeKuttaTrack(Double_t     dr    = 0.,
                    Double_t     phi0  = 0.,
                    Double_t     kappa = 1.e-5,
                    Double_t     dz    = 0.,
                    Double_t     tanl  = 0.,
                    Double_t     x0    = 0.,
                    Double_t     y0    = 0.,
                    Double_t     z0    = 0.,
                    Double_t     b     = 3.,
					TTrackFrame *fp    = 0);

   TRungeKuttaTrack(const TMatrixD    &a,
                    const TVector3    &x0,
                          Double_t     b  = 3.,
						  TTrackFrame *fp = 0);
                          

   TRungeKuttaTrack(const TVector3 &x1,
                    const TVector3 &x2,
                    const TVector3 &x3,
                          Double_t b   = 3.,
                          Bool_t   dir = kIterForward);

   
   TRungeKuttaTrack(const TMatrixD    &a,
                    const TVector3    &x0,     //x0 is a local pivot
                          Double_t     b,
                    const TTrackFrame &frame);

   TRungeKuttaTrack(const Double_t  chg,
		            const TVector3 &x,
					const TVector3 &p);

   virtual ~TRungeKuttaTrack() {}

   void SetToTrack(THelicalTrack& heltrack) const;
   void SetFromTrack(THelicalTrack& heltrack);

   // Utility methods
   virtual void MoveTo(const TVector3 &globalPivot, 
                             Double_t &step,     
		                     TMatrixD  *FPtr = 0,    
		                     TMatrixD  *F12Ptr = 0,    
		                     TMatrixD  *CPtr = 0);

   virtual void MoveTo(const TVector3    &globalPivot, 
                             Double_t    &step,     
		                     TKalMatrix  &FPtr);

   virtual TVector3 CalcXAt   (Double_t h) const;
   virtual TMatrixD CalcDxDa  (Double_t h) const;
   virtual TMatrixD CalcDxDphi(Double_t h) const;
   
   inline virtual Double_t GetPath(Double_t p) const { return p; }

   // Setter
   inline void SetCurPosition(TVector3 x) { fPosition = x; }
   inline void SetCurMomentum(TVector3 p) { fMomentum = p; }
   inline void SetCharge     (Double_t q) { fCharge   = q; }
   inline void SetAlpha      (Double_t a) { fAlpha    = a; }

   // Getter
   inline TVector3 GetCurPosition() const { return fPosition; }
   inline TVector3 GetCurMomentum() const { return fMomentum; }

   inline Double_t GetCharge()      const { return fCharge;   }

   // Take a step 
   void StepRungeKutta(Double_t  step, 
                       TVector3& vx, 
                       TVector3& vp) const;

   // Take a step, update current position and momentum of this track
   void StepRungeKutta(Double_t step);

   // Update momentum and position from current state vector
   void UpdatePX();

   // Calculate state vector from momentum and position
   void PXToSV(TVector3& p, TVector3& x, TKalMatrix& sv);

private:
   void CalcAXK(vector<TVector3>& A,
			    vector<TVector3>& X,
			    vector<TVector3>& K,
				Double_t          h) const;

   void CalcDKDh(vector<TKalMatrix>& DKDh, vector<TVector3>& K, Double_t h) const;
   void CalcDKDA(vector<TKalMatrix>& DKDA, vector<TVector3>& K, Double_t h) const;

   void CalcDKDa(vector<TKalMatrix>& DKDa, 
		         vector<TVector3>&   K, 
				 Double_t            h) const;
   
   void CalcDKDx(vector<TKalMatrix>& DKDx, 
		         vector<TVector3>&   K, 
				 Double_t            h) const;

   TKalMatrix CalcDpxDa() const;
   TKalMatrix CalcDpDa () const;

   TKalMatrix CalcDpxDpx(Double_t  h) const;
   TKalMatrix CalcDxDp  (Double_t  h) const;
   TKalMatrix CalcDaDpx (TVector3& p) const;

   void CalcXPAt(Double_t step, TVector3& vx, TVector3& vp) const;

   void SVToP (TVector3& p);

   TVector3 GetLocalBfield(TVector3 x0) const;

private:
   Double_t fCharge;     

   TVector3 fPosition{};
   TVector3 fMomentum{};

   Double_t fSkappa{};          // sign of kappa
   Double_t fLambda{};

   ClassDef(TRungeKuttaTrack,1)   
};

#endif // TRUNGEKUTTATRACK_H
