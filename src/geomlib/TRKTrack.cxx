//*************************************************************************
//* ====================
//*  TRKTrack Class
//* ====================
//*
//* (Description)
//* (Requires)
//* (Provides)
//*
//*************************************************************************
//
#include <iostream>
#include "TRKTrack.h"
#include "TBField.h"
#include "TRKMagField.h"

ClassImp(TRKTrack)

TRKTrack::TRKTrack(Double_t chg, TVector3 x, TVector3 p)
	     :fPosition(x),fMomentum(p)
{
	SetCharge(int(chg));
	SetMagFieldObj(new TRKMagField());
}

TRKTrack::TRKTrack(const TRKTrack&  orig)
	     :fPosition(orig.fPosition), fMomentum(orig.fMomentum)
{
	SetCharge(int(orig.GetCharge()));
	SetMagFieldObj(new TRKMagField());
}

// Utility methods
void TRKTrack::StepRungeKutta(Double_t step)
{
    const Double_t MM2CM = 0.1;
	const Double_t CM2MM = 10;

	//fPosition.Print();
	//fMomentum.Print();

	//momentum
	Double_t p   = fMomentum.Mag();
	Double_t px1 = fMomentum.X()/p;
	Double_t py1 = fMomentum.Y()/p;
	Double_t pz1 = fMomentum.Z()/p;
	
	//input and output vectors
	Double_t vect[7] = {fPosition.X()*MM2CM, fPosition.Y()*MM2CM, fPosition.Z()*MM2CM, 
		                px1, py1, pz1, p};
	Double_t vout[7];

	//step length
    Double_t aStep =  step * MM2CM;

	//RK
	TEveTrackPropagator::StepRungeKutta(aStep, vect, vout);
	
	fMomentum.SetXYZ(vout[3]*vout[6], vout[4]*vout[6], vout[5]*vout[6]);
	fPosition.SetXYZ(vout[0]*CM2MM, vout[1]*CM2MM, vout[2]*CM2MM);

	//std::cout << std::endl;
}

Double_t TRKTrack::GetCharge() const
{
    return TEveTrackPropagator::fH.fCharge;
}

void TRKTrack::SetCharge(Int_t chg)
{
    TEveTrackPropagator::fH.fCharge = chg;
}
