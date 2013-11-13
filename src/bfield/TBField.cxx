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

#include "TBField.h"   // from KalTrackLib
#include <memory>            // from STL
#include <cmath>
#include <iostream>          // from STL

ClassImp(TBField)

//TEveMagField* TBField::fMagField = 0;
//TVector3      TBField::fBfield   = TVector3(0,0,3);

Bool_t   TBField::fUseUniformBfield = kFALSE;
Double_t TBField::fFieldCoeff       = 3.;

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
#if 0	   
	 if(fUseUniformBfield)
	 {
		 return 10*TVector3(0,0,3);
	 }

	 //const Double_t coeffxy = 0.0008;
	 //const Double_t coeffxy = 0.00001;
     const Double_t coeffxy = 0.01;
     //const Double_t coeffxy = 0.001;

	 Double_t bx = coeffxy * globalPosition.X();
	 Double_t by = coeffxy * globalPosition.Y();

	 Double_t bt = sqrt(bx*bx+by*by);
	 Double_t bz = sqrt(fabs((3 - bt)*(3+bt)));
	 //bz = 0;

	 TVector3 bfield(bx,by,bz);
	 bfield *= 10;
#endif

#if 0
	 Double_t B0   = 3.5;

	 //return TVector3(0,0,B0);

	 if(fUseUniformBfield)
	 {
		 return TVector3(0,0,B0);
	 }

	 Double_t R    = 1500.; //mm
	 Double_t Br   = globalPosition.Perp2()/R/R;
	 Double_t Bz   = B0 * sqrt(1 - Br*Br);
	 Br *= B0;

     Double_t phi = globalPosition.Phi();
	 Double_t Bx  = Br * cos(phi);
	 Double_t By  = Br * sin(phi);
	 Bx = By = 0;

	 TVector3 bfield(Bx, By, Bz);
#else
    const Double_t zmax    = 3000.;
    const Double_t rmax    = 3000.;
    const Double_t coeffxy = fFieldCoeff / (zmax * rmax);
    const Double_t bmax    = 3.; // B at the origin

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
#endif

#if 0
	 std::cout << "b field at " << "(" 
		  << globalPosition.X() << ","
		  << globalPosition.Y() << ","
		  << globalPosition.Z() << "):   ("
		  << bfield.X() << "," << bfield.Y() << "," << bfield.Z() << ")" << std::endl;
#endif

#if 0
	 if(globalPosition.Perp()<650)
		 return TVector3(0,0,3.5);
	 else
		 return TVector3(0,0,1);
#endif

	 return bfield;
}
