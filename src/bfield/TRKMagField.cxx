//*************************************************************************
//* ====================
//*  TRKMagField Class
//* ====================
//*
//* (Description)
//* (Requires)
//* (Provides)
//*
//*************************************************************************
//
#include "TRKMagField.h"
#include "TVector3.h"
#include "TBField.h"

#include <iostream>
#include <iomanip>

ClassImp(TRKMagField)

//_____________________________________________________________________
//  -----------------------------------
//  Runge-Kutta magnetic field Class
//  -----------------------------------
//

TEveVectorD TRKMagField::GetFieldD(double x, double y, double z) const
{ 
	//const Double_t CM2MM = 10.;
	const Double_t CM2MM = 1.;

	TVector3 bfield = TBField::GetGlobalBfield(TVector3(x*CM2MM, y*CM2MM, z*CM2MM));

    return TEveVectorD(bfield.x(), bfield.y(), bfield.z());
}

TVector3 TRKMagField::GetField(TVector3 x) 
{ 
	//const Double_t CM2MM = 10.;
	const Double_t CM2MM = 1.;

	TVector3 bfield = TBField::GetGlobalBfield( x*CM2MM );

    return bfield;
}
