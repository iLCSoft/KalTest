//*************************************************************************
//* =====================
//*  TBField Class
//* =====================
//*
//* (Description)
//*   A singleton class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//*     TObject
//* (Provides)
//*     class TBField
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi  	Original version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*   2010/04/06  K.Fujii        Modified Transport() to allow a 1-dim hit,
//*                              for which pivot is at the expected hit.
//*
//*************************************************************************

#include <memory>      // from STL
#include <cmath>
#include <iostream>    // from STL

#include "TBField.h"   // from KalTrackLib

ClassImp(TBField)

TEveMagField* TBField::fField = 0;

//Bool_t   TBField::fUseUniformBfield = kFALSE;
Bool_t   TBField::fUseUniformBfield = kTRUE;
Double_t TBField::fFieldCoeff       = 1.;

//_________________________________________________________________________
//  ----------------------------------
//   Ctors and Dtor
//  ----------------------------------

TBField::TBField()
{
}

TBField::~TBField()
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________

TVector3 TBField::GetGlobalBfield(const TVector3& globalPosition)
{
	if(!fUseUniformBfield && fField!=0) {
		TEveVectorD vv = fField->GetFieldD(globalPosition.X(), 
				                          globalPosition.Y(), 
										  globalPosition.Z());

		return TVector3(vv[0], vv[1], vv[2]);
	}

	// ILD uniform magnetic field
	Double_t B0   = 3.5; 

	if(fUseUniformBfield)
	{
		return TVector3(0,0,B0);
	}

	// an artificial non-uniform magnetic field
    const Double_t zmax    = 3000.;
    const Double_t rmax    = 3000.;
    const Double_t coeffxy = fFieldCoeff / (zmax * rmax);
    const Double_t bmax    = 3.5; // B at the origin

    if(fUseUniformBfield)
    {
        return TVector3(0,0,bmax);
    }

    Double_t xg = globalPosition.X();
    Double_t yg = globalPosition.Y();
    Double_t zg = globalPosition.Z();

    Double_t bx = bmax * coeffxy * zg * xg; 
    Double_t by = bmax * coeffxy * zg * yg; 
    Double_t bz = bmax * (1. - coeffxy * zg * zg);

	TVector3 bfield(bx, by, bz);

	return bfield;
}
