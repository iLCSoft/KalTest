#ifndef TTRACKFRAME_H
#define TTRACKFRAME_H
//*************************************************************************
//* ===================
//*  TTrackFrame Class
//* ===================
//*
//* (Description)
//*    TTrackFrame is a class to transform 3D vector and state vector,
//*    and propagation matrix.
//* (Requires)
//* 	TRotation
//*     TVector3
//* (Provides)
//* 	class TTrackFrame
//*
//*************************************************************************

#include "TKalMatrix.h"
#include "TRotation.h"
#include "TVector3.h"


class TTrackFrame { 
public:

   enum TRType {
	    kLocalToLocal = 0,
		kLocalToGlobal,
		kGlobalToLocal
   };

   TTrackFrame(){}

   //Copy constructor
   TTrackFrame(const TTrackFrame &orig);

   //Build frame by last frame, shift vector and B field. The B field determines local
   //rotation matrix.
   TTrackFrame(const TTrackFrame &lastFrame, const TVector3 &v, const TVector3 &b);

   //TTrackFrame& operator=( const TTrackFrame& frame );

   virtual ~TTrackFrame() {}

   //Transform 3D position vector
   TVector3 Transform(const TVector3 &v, TRType tt = kLocalToLocal) const;

   TVector3 TransformBfield(const TVector3 &b, TRType tt = kLocalToLocal);

   //Transform TKalMatrix(such as state vector, propagation matrix, and covariance matrix),
   //this is generally local to local transformation.
   void       Transform(TKalMatrix* sv, TKalMatrix* Fr = 0);

   inline TRotation GetRotation() const { return fRotMat; }
   inline TVector3  GetShift   () const { return fShift;  }

   inline void  SetRotation (const TRotation& r) { fRotMat = r; }
   inline void  SetShift    (const TVector3&  v) { fShift  = v; }

private:

   TKalMatrix CalcRotationMatrix(const TKalMatrix& b);

   void       CalcDtrDt         (TKalMatrix& dtrdt, const TKalMatrix& R);
   void       CalcDappDtr       (TKalMatrix& dappdtr, TKalMatrix& tr, Double_t drhosign,
                                 Double_t cpasign, Int_t  mode = 0);
   void       CalcDtDap         (TKalMatrix& dtdap, TKalMatrix& a, Double_t cpasign, Int_t mode = 0);

private:
   //From the global frame to this frame
   TRotation  fRotMat{};
   TVector3   fShift{};

   //From last local frame to this frame
   TRotation  fDeltaRotMat{};
   TVector3   fDeltaShift{};

   ClassDef(TTrackFrame,1)
};
#endif
