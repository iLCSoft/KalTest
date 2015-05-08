#ifndef TRKMAGFIELD_H
#define TRKMAGFIELD_H
//*************************************************************************
//* ====================
//*  TRKMagField Class
//* ====================
//*
//* (Description)
//*   A class to implement a magnetic field for Runge-Kutta track object.
//*
//*
//* (Requires)
//*     TEveMagField 
//* (Provides)
//*     class TRKMagField
//*
//*************************************************************************
//

#include "TEveTrackPropagator.h"

class TRKMagField : public TEveMagField
{
public:
    TRKMagField() {}
	//virtual ~TRKMagField();
 
    virtual TEveVectorD  GetFieldD(double x, double y, double z) const;
    static  TVector3     GetField ( TVector3 x);

private:

   ClassDef(TRKMagField,1)   
};

#endif

