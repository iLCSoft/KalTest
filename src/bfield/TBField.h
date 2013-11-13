#ifndef TBFIELD_H
#define TBFIELD_H
//*************************************************************************
//* =====================
//*  TBField Class
//* =====================
//*
//* (Description)
//*   A sigleton to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObject
//* (Provides)
//* 	class TBField
//* (Update Recored)
//*   2013/01/31  Bo Li	 Original Version.
//*
//*************************************************************************

#include "TObject.h"     // from ROOT
#include "TVector3.h"
//#include "TEveTrackPropagator.h"     // from ROOT


//_____________________________________________________________________
//  ------------------------------
//   Detector system class
//  ------------------------------

class TBField : public TObject {
public:
   TBField();
   virtual ~TBField();

   // Utility methods
   static TVector3 GetGlobalBfield(const TVector3& globalPosition);

   //static void SetBfieldPtr(TEveMagField* b){ fMagField = b; }
   //static void SetBfield   (TVector3      b){ fBfield   = b; }
   static void SetUseUniformBfield(Bool_t b) { fUseUniformBfield = b; }
   static void SetBfieldCoeff(Double_t k)    { fFieldCoeff       = k; }

private:
   //static    TVector3      fBfield;     //a constant b field

   //static    TEveMagField* fMagField;   //a field map
   static Bool_t   fUseUniformBfield;
   static Double_t fFieldCoeff;

   ClassDef(TBField,1)  // Base class for detector system
};

#endif
