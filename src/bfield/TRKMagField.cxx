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
#if 0
virtual TRKMagField::~TRKMagField()
{
	std::cout << "Delete TRKMagField" << std::endl;

	TEveMagField::~TEveMagField();
}
#endif

TEveVectorD TRKMagField::GetFieldD(double x, double y, double z) const
{ 
	//const Double_t CM2MM = 10.;
	const Double_t CM2MM = 1.;

	TVector3 bfield = TBField::GetGlobalBfield(TVector3(x*CM2MM, y*CM2MM, z*CM2MM));
#if 0
	std::cout << std::setprecision(15) << "(" << x << "," << y << "," << z << "), "
		 << "bfield: (" 
		 << bfield.X() << ", " 
		 << bfield.Y() << ", "
		 << bfield.Z() << ")"
		 << std::endl;
#endif
    return TEveVectorD(bfield.x(), bfield.y(), bfield.z());
}

TVector3 TRKMagField::GetField(TVector3 x) 
{ 
	//const Double_t CM2MM = 10.;
	const Double_t CM2MM = 1.;

	TVector3 bfield = TBField::GetGlobalBfield( x*CM2MM );
#if 0
	std::cout << std::setprecision(15) << "(" << x << "," << y << "," << z << "), "
		 << "bfield: (" 
		 << bfield.X() << ", " 
		 << bfield.Y() << ", "
		 << bfield.Z() << ")"
		 << std::endl;
#endif
    return bfield;
}
