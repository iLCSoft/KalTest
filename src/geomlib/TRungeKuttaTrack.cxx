//*************************************************************************
//* ====================
//*  TRungeKuttaTrack Class
//* ====================
//*
//* (Description)
//* (Requires)
//* (Provides)
//*
//*************************************************************************
//
#include <iostream>
#include <cmath>

#include "TRungeKuttaTrack.h"
#include "TBField.h"
#include "TRKMagField.h"

ClassImp(TRungeKuttaTrack)

//_____________________________________________________________________
//  -----------------------------------
//  Helical Track Class
//  -----------------------------------
//
//_____________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//

static Double_t hk[4]  = {0., 0.5,   0.5,   1. };
static Double_t h2k[4] = {0., 0.125, 0.125, 0.5};

//All the parameters are local
TRungeKuttaTrack::TRungeKuttaTrack(Double_t dr,
                             Double_t phi0,
                             Double_t kappa,
                             Double_t dz,
                             Double_t tanl,
                             Double_t x0,
                             Double_t y0,
                             Double_t z0,
                             Double_t b,
							 TTrackFrame *fp)
             : THelicalTrack(dr,phi0,kappa,dz,tanl, x0,y0,z0, b),fCharge(1.)
{
    if (fp) fFrame = *fp;
    UpdatePX();
}

TRungeKuttaTrack::TRungeKuttaTrack(const TMatrixD &a,
                                   const TVector3 &x0,
                                         Double_t  b,
								      TTrackFrame *fp)
             : THelicalTrack(a, x0, b, fp), fCharge(1)
{
}

TRungeKuttaTrack::TRungeKuttaTrack(const TMatrixD    &a,
                                   const TVector3    &x0,     //x0 is a local pivot
                                         Double_t     b,
                                   const TTrackFrame &frame)
             : THelicalTrack(a, x0, b, frame), fCharge(1.)
{
}

TRungeKuttaTrack::TRungeKuttaTrack(const TVector3 &x1, 
                                   const TVector3 &x2,
                                   const TVector3 &x3,
                                         Double_t  b,
                                         Bool_t    dir)
                 : fCharge(1)
{
   SetMagField(b);
   CalcStartHelix(x1,x2,x3,dir);
}


TRungeKuttaTrack::TRungeKuttaTrack(const Double_t chg,
		                           const TVector3 &x,
								   const TVector3 &p)
                 :fCharge(chg),fPosition(x),fMomentum(p)
{	
	fSkappa = TMath::Sign(1., chg);
}

//_____________________________________________________________________
//  ----------------
//  Utility methods
//  ----------------

void TRungeKuttaTrack::SetToTrack(THelicalTrack& heltrack) const
{
	TKalMatrix sv(kSdim, 1);

	PutInto(sv);
	heltrack.SetTo(sv, GetPivot());
	heltrack.SetFrame(GetFrame());
	heltrack.SetMagField(GetMagField());
}

void TRungeKuttaTrack::SetFromTrack(THelicalTrack& heltrack)
{
	TKalMatrix sv(kSdim, 1);

	heltrack.PutInto(sv);
	SetTo(sv, heltrack.GetPivot());
	SetAlpha(heltrack.GetPtoR());
	SetFrame(heltrack.GetFrame());
	SetMagField(heltrack.GetMagField());

	UpdatePX();
}

void TRungeKuttaTrack::MoveTo(const TVector3    &globalPivot, 
                                    Double_t    &step,     
	                                TKalMatrix  &rkDF)
{
	TVector3 xv0to = fFrame.Transform(globalPivot, TTrackFrame::kGlobalToLocal);

	TKalMatrix F12(6,5);
	MoveTo(xv0to, step, &rkDF, &F12);

	double p0 = fPhi0;
	TVector3 pv = GetCurPosition();

	double point_x = pv.X();
	double point_y = pv.Y();

	double rho = GetRho();
	double cpa = GetKappa();
	double tanl = GetTanLambda();
	double tanp = tan(p0);
	double sinp = sin(p0);
	double cosp = cos(p0);

	double xc = point_x + rho * cosp;
	double yc = point_y + rho * sinp;

	double dxcdkp = -rho / cpa * cosp;
	double dycdkp = -rho / cpa * sinp;
	
	double dxcdpp = rho * (-sinp);
	double dycdpp = rho * cosp;

    ////////////////////////////////////////// 
    double dxpdk = F12(3,2);
	double dypdk = F12(4,2);
	double dkpdk = rkDF(2,2);
	double dpdk  = rkDF(1,2);    //@phi^{prime}/@kappa

	double dzdk  = F12(5,2);

	double dxcdk = dxpdk + dxcdkp*dkpdk + dxcdpp*dpdk;
	double dycdk = dypdk + dycdkp*dkpdk + dycdpp*dpdk;
	double dppdk = cosp*cosp*(dycdk/(xc-point_x) - tanp*dxcdk/(xc-point_x));

	double dxpdphi0 = F12(3,1); 
	double dypdphi0 = F12(4,1);
	double dkpdphi0 = rkDF(2,1);
	double dphipdphi0 = rkDF(1,1);

	double dxcdphi0 = dxpdphi0 + dxcdkp*dkpdphi0 + dxcdpp*dphipdphi0;
	double dycdphi0 = dypdphi0 + dycdkp*dkpdphi0 + dycdpp*dphipdphi0;
	double dppdphi0 = cosp*cosp *(dycdphi0/(xc-point_x) - tanp*dxcdphi0/(xc-point_x));

	double ddrdphi0 = dxcdphi0*cosp - (xc-point_x)*sinp*dppdphi0 +
			                dycdphi0*sinp + (yc-point_y)*cosp*dppdphi0 +
							rho/cpa* dkpdphi0;

	double ddrdk = dxcdk*cosp - (xc-point_x)*sinp*dppdk +
		            dycdk*sinp + (yc-point_y)*cosp*dppdk +
					rho/cpa* dkpdk;


	double ddzdk = dzdk - rho * tanl * (dppdk-dpdk);

	double dxpddr = F12(3,0);
	double dypddr = F12(4,0);
	double dkpddr = rkDF(2,0);
	double dphipddr = rkDF(1,0);

	double dxcddr = dxpddr + dxcdkp*dkpddr + dxcdpp*dphipddr;
	double dycddr = dypddr + dycdkp*dkpddr + dycdpp*dphipddr; 
	double dphi0ddr = cosp*cosp*(dycddr/(xc-point_x) - tanp * dxcddr/(xc-point_x)); 

	double ddrddr = dxcddr * cosp - (xc-point_x) * sinp * dphi0ddr 
			            + dycddr * sinp + (yc-point_y) * cosp * dphi0ddr 
						+ rho/cpa * dkpddr;

	double dzpddr = F12(5,0);
	double ddzddr = dzpddr - rho * tanl * (dphi0ddr - dphipddr);

	double dzpdphi0 = F12(5,1);
    double dphi0pdphi0 = cosp*cosp*(dycdphi0/(xc-point_x) - tanp * dxcdphi0/(xc-point_x));
	double ddzdphi0 = dzpdphi0 - rho * tanl * (dphi0pdphi0 - 1.);

	double dxpddz = F12(3,3);
	double dypddz = F12(4,3);
	double dkpddz = rkDF(2,3);
	double dppddz = rkDF(1,3);

	double dxcddz = dxpddz + dxcdkp*dkpddz + dxcdpp*dppddz;
	double dycddz = dypddz + dycdkp*dkpddz + dycdpp*dppddz;
	double dphi0ddz = cosp*cosp*(dycddz/(xc-point_x) - tanp * dxcddz/(xc-point_x));

	double ddrddz = dxcddz * cosp - (xc-point_x) * sinp * dphi0ddz 
                    + dycddz * sinp + (yc-point_y) * cosp * dphi0ddz 
                    + rho/cpa * dkpddz;

	double dzddz  = F12(5,3);
	double ddzddz = dzddz - rho * tanl * (dphi0ddz-dppddz);

	double dxpdtanl = F12(3,4);
	double dypdtanl = F12(4,4);
	double dkpdtanl = rkDF(2,4);
	double dppdtanl = rkDF(1,4);

	double dxcdtanl = dxpdtanl + dxcdkp*dkpdtanl + dxcdpp*dppdtanl;
	double dycdtanl = dypdtanl + dycdkp*dkpdtanl + dycdpp*dppdtanl;
	double dphi0dtanl = cosp*cosp*(dycdtanl/(xc-point_x) - tanp * dxcdtanl/(xc-point_x));

	double dzdtanl  = F12(5,4);
	double ddzdtanl = dzdtanl - rho * tanl * (dphi0dtanl-dppdtanl);

	double dkddr = dkpddr;

	rkDF(0,0) = ddrddr;
	rkDF(0,1) = ddrdphi0;
	rkDF(0,2) = ddrdk;
	rkDF(0,3) = ddrddz;
	rkDF(0,4) = ddrddz;

	rkDF(1,0) = dphi0ddr;
	rkDF(1,1) = dppdphi0;
	rkDF(1,2) = dppdk;
	rkDF(1,3) = dphi0ddz;
	rkDF(1,4) = dphi0dtanl;

	rkDF(2,0) = dkddr; //1*dkddr
	rkDF(2,1) = dkpdphi0; //1*dkpdphi0
	rkDF(2,2) = dkpdk; //1*dkpdk
	rkDF(2,3) = dkpddz; //1*dkpddz
	rkDF(2,4) = dkpdtanl; //1*dkpdtanl

	rkDF(3,0) = ddzddr;
	rkDF(3,1) = ddzdphi0;
	rkDF(3,2) = ddzdk;
	rkDF(3,3) = ddzddz;
	rkDF(3,4) = ddzdtanl;

    // dtanlda is trival, so rkDF(4,:) is zero

	/////////////////////////////
	// transformation of frame
	////////////////////////////
    TKalMatrix av(5,1);
	PutInto(av);

    Int_t sdim = rkDF.GetNrows();
    TKalMatrix Fr(sdim,sdim);
	
	TVector3 globalBfield = TBField::GetGlobalBfield(globalPivot);
    TTrackFrame frameOfNewPivot(fFrame, xv0to, globalBfield);
	
    SetFrame(frameOfNewPivot);
    SetMagField(globalBfield.Mag());

    TVector3 pivot = frameOfNewPivot.Transform(xv0to, TTrackFrame::kLocalToLocal); 

	frameOfNewPivot.Transform(&av, &Fr);
	SetTo(av, pivot);
    
	rkDF = Fr * rkDF;
}

void TRungeKuttaTrack::MoveTo(const TVector3 & /*globalPivot */ , // new global pivot
                                    Double_t &step,        // step, an input parameter
     		                        TMatrixD  *FPtr,       // propagator matrix
     		                        TMatrixD  *F12Ptr,     // intermediate propagator matrix 
		                            TMatrixD  *CPtr)       // covariance matrix
{
   // --------------------------------------------------
   // Use Runge-Kutta method to propagate
   // --------------------------------------------------

   //vector to store the position and momentum after stepping
   TVector3 vx, vp;
   TKalMatrix av(5,1);

   //the track position and momentum is calculated, 
   //and both should be global. Then the new values
   //are got by stepping the track.

   TKalMatrix dpda  = CalcDpxDa();

   TKalMatrix dpdp  = CalcDpxDpx(step);

   CalcXPAt(step, vx, vp);
   SetCurPosition(vx);
   SetCurMomentum(vp);

   TKalMatrix dadp = CalcDaDpx(vp);

   PXToSV(vp, vx, av);
   
   SetTo(av, vx);

   if (!FPtr && !CPtr) {
      return;
   }

   TMatrixD Fdummy(5,5);
   TMatrixD &F = FPtr ? *FPtr : Fdummy;

   TMatrixD F12dummy(6,5);
   TMatrixD &F12 = F12Ptr ? *F12Ptr : F12dummy;

   TMatrixD temp(5,5);

   temp.ResizeTo(F);
   temp = dadp * dpdp * dpda;
   F = temp;

   if(F.GetNrows() == 6) F(5,5) = 1.;

   F12 = dpdp * dpda;

   if (!CPtr) {
      return;
   }

   // ---------------------------------------------------
   // Calculate C' = C^k-1_k
   // ---------------------------------------------------
   TMatrixD &C  = *CPtr;
   TMatrixD  Ft = TMatrixD(TMatrixD::kTransposed, F);
   TMatrixD  Cp = F * C * Ft;
   C = Cp;
}

TVector3 TRungeKuttaTrack::CalcXAt(Double_t step) const
{
    TVector3 vx = GetCurPosition();	
	TVector3 vp = GetCurMomentum();

	//Runge-Kutta method 
	StepRungeKutta(step, vx, vp);

    //Transform local coodinate to global
    TVector3 globalX = fFrame.Transform(vx, TTrackFrame::kLocalToGlobal);
	
	return globalX;
}

TMatrixD TRungeKuttaTrack::CalcDxDa(Double_t h) const
{
	TKalMatrix dxdpx  = CalcDxDp(h);
	dxdpx.ResizeTo(3,6);
	dxdpx(0,3) = dxdpx(1,4) = dxdpx(2,5) = 1.;

	TKalMatrix dpxda = CalcDpxDa();
	TKalMatrix dxda  = dxdpx * dpxda;

    //Tramsform
    const TRotation& rotMat = fFrame.GetRotation();
    TKalMatrix invRotMat(rotMat.Inverse());
    dxda = invRotMat * dxda;

    return dxda;
}

TMatrixD TRungeKuttaTrack::CalcDxDphi(Double_t h) const
{
    vector<TKalMatrix> DKDh;
    vector<TVector3> K;

	CalcDKDh(DKDh, K, h);

	TVector3 um = GetCurMomentum().Unit();

	TKalMatrix dxdh = TKalMatrix::ToKalMat(GetCurMomentum().Unit()) + 
		              h/3. * TKalMatrix::ToKalMat(K[0] + K[1] + K[2]) + 
					  h*h/6. * (DKDh[0] + DKDh[1] + DKDh[2]);

    //Transform
    const TRotation& rotMat = fFrame.GetRotation();
    TKalMatrix invRotMat(rotMat.Inverse());
    dxdh = invRotMat * dxdh;
	
	return dxdh;
}

TKalMatrix TRungeKuttaTrack::CalcDpDa() const
{        
	// The state vector, a, could be the predicted or the filtered one.
	
    Double_t chg = 1.;    //the absolute value of particle charge
 
    Double_t phi0 = fPhi0;
    Double_t cpa  = fKappa;
    Double_t tanl = fTanL;
 
    Double_t sinphi0 = sin(phi0);
    Double_t cosphi0 = cos(phi0);
    Double_t cpa2    = cpa*cpa;
    Double_t acpa    = fabs(cpa);

	TKalMatrix dpda(3,5);

    dpda(0,1) = -chg/acpa*cosphi0;
    dpda(0,2) = chg*fSkappa/cpa2*sinphi0;
   
    dpda(1,1) = -chg/acpa*sinphi0;
    dpda(1,2) = -chg*fSkappa/cpa2*cosphi0;
   
    dpda(2,2) = -chg*fSkappa/cpa2*tanl;
    dpda(2,4) = chg/acpa; 

	return dpda;
}
    
TKalMatrix TRungeKuttaTrack::CalcDpxDa() const
{        
	// The state vector, a, could be the predicted or the filtered one.
	
    Double_t chg = 1.;    //the absolute value of particle charge
 
    Double_t phi0 = fPhi0;
    Double_t cpa  = fKappa;
    Double_t tanl = fTanL;
	Double_t drho = fDrho;

    Double_t sinphi0 = sin(phi0);
    Double_t cosphi0 = cos(phi0);
    Double_t cpa2    = cpa*cpa;
    Double_t acpa    = fabs(cpa);

	TKalMatrix dpxda(6,5);
	TKalMatrix dpda(3,5);
	dpxda.Zero();
	dpda.Zero();

    dpda(0,1) = dpxda(0,1) = -chg/acpa*cosphi0;
    dpda(0,2) = dpxda(0,2) = chg*fSkappa/cpa2*sinphi0;
   
    dpda(1,1) = dpxda(1,1) = -chg/acpa*sinphi0;
    dpda(1,2) = dpxda(1,2) = -chg*fSkappa/cpa2*cosphi0;
   
    dpda(2,2) = dpxda(2,2) = -chg*fSkappa/cpa2*tanl;
    dpda(2,4) = dpxda(2,4) = chg/acpa; 

	dpxda(3,0) = cosphi0;
	dpxda(3,1) = -drho * sinphi0;

	dpxda(4,0) = sinphi0;
	dpxda(4,1) = drho * cosphi0;

	dpxda(5,3) = 1.;

	return dpxda;
}

TKalMatrix TRungeKuttaTrack::CalcDpxDpx(Double_t h) const
{
	vector<TKalMatrix> DKDa;
	vector<TKalMatrix> DKDx;
	vector<TVector3> K;

	CalcDKDa(DKDa, K, h);
	CalcDKDx(DKDx, K, h);

	TKalMatrix dpdp(3,3);
	TKalMatrix dppdp(3,3);

	TKalMatrix u(3,3);
	u.UnitMatrix();

	dpdp = u + 1./6 * h * (DKDa[0] + 2. * DKDa[1] + 2. * DKDa[2] + DKDa[3]);
	
	dppdp = dpdp;
	
	TKalMatrix dpdx(3,3);
	dpdx.Zero();

    TKalMatrix dxdp(3,3);
	dxdp = CalcDxDp(h);

	TKalMatrix dxdx(3,3);
    dxdx.UnitMatrix();

	TKalMatrix dpxdpx(6,6);
	
	dpdp.ResizeTo(dpxdpx);
	dpxdpx = dpdp;

	for(Int_t i=0; i<3; i++)
	{
		for(Int_t j=0; j<3; j++)
		{
			dpxdpx(i,j+3) = dpdx(i,j);
		}
	}

	for(Int_t i=0; i<3; i++)
	{
		for(Int_t j=0; j<3; j++)
		{
			dpxdpx(i+3,j) = dxdp(i,j);
		}
	}

	for(Int_t i=0; i<3; i++)
	{
		for(Int_t j=0; j<3; j++)
		{
			dpxdpx(i+3,j+3) = dxdx(i,j);
		}
	}

	return dpxdpx;
}

TKalMatrix TRungeKuttaTrack::CalcDaDpx(TVector3& p) const
{ 
   TKalMatrix dadpx(5,6);
   TKalMatrix dadp(5,3);

   Double_t chg  = 1.;

   Double_t px = p.X();
   Double_t py = p.Y();
   Double_t pz = p.Z();

   Double_t pt  = sqrt(px*px + py*py);
   Double_t pt2 = pt*pt;
   Double_t pt3 = pt*pt2;

   Double_t phi0 = atan2(-p.X(), p.Y());
   Double_t cosphi0 = cos(phi0);

   dadpx(0,3) = 1./cosphi0;

   dadp(1,0) = dadpx(1,0) = -py/pt2;
   dadp(1,1) = dadpx(1,1) = px/pt2;

   dadp(2,0) = dadpx(2,0) = -fSkappa*chg*px/pt3;
   dadp(2,1) = dadpx(2,1) = -fSkappa*chg*py/pt3;

   dadpx(3,5) = 1.;

   dadp(4,0) = dadpx(4,0) = -px*pz/pt3;
   dadp(4,1) = dadpx(4,1) = -py*pz/pt3;
   dadp(4,2) = dadpx(4,2) = 1./pt;

   return dadpx;
}

TVector3 TRungeKuttaTrack::GetLocalBfield(TVector3 x0) const
{
	TVector3 globalx0 = fFrame.Transform(x0, TTrackFrame::kLocalToGlobal);
	TVector3 globalbf = TRKMagField::GetField(globalx0);
	return fFrame.TransformBfield(globalbf, TTrackFrame::kGlobalToLocal);
}

void TRungeKuttaTrack::CalcAXK(vector<TVector3>& A, 
		                       vector<TVector3>& X,
							   vector<TVector3>& K,
							   Double_t          h) const
{
	A.clear();
	X.clear();
	K.clear();

	TVector3 a0 = GetCurMomentum().Unit();
	TVector3 x0 = GetCurPosition();

	Double_t lambda = fCharge * 1./GetCurMomentum().Mag() * 2.99792458e-4;

	TVector3 ai;
	TVector3 xi;

	TVector3 ki = lambda * a0.Cross(GetLocalBfield(x0));

	for(Int_t i=0; i<4; i++)
	{ 
		ai = a0 + hk[i] * h * ki;
		xi = x0 + hk[i] * h * a0 + h2k[i] * h * h * ki;
	    
		ki = lambda * ai.Cross(GetLocalBfield(xi));
	
		A.push_back(ai);
		X.push_back(xi);
		K.push_back(ki);
	}
}

void TRungeKuttaTrack::CalcDKDA(vector<TKalMatrix>& DKDA, 
		                        vector<TVector3>&   K, 
								Double_t            h) const
{
	DKDA.clear();

	vector<TVector3> A;
	vector<TVector3> X;

	CalcAXK(A, X, K, h);

	TKalMatrix dkda(3,3);

	for(Int_t i=0; i<4; i++)
	{
		dkda.Zero();

		TVector3 b = GetLocalBfield(X[i]); 

		dkda(0,1) =  b.Z();
	   	dkda(0,2) = -b.Y();
	   	dkda(1,0) = -b.Z();
	   	dkda(1,2) =  b.X();
	   	dkda(2,0) =  b.Y();
	   	dkda(2,1) = -b.X();

		dkda *= fCharge * 1./GetCurMomentum().Mag() * 2.99792458e-4;

		DKDA.push_back(dkda);	
    }
}

void TRungeKuttaTrack::CalcDKDa(vector<TKalMatrix>& DKDa, 
		                        vector<TVector3>&   K, 
								Double_t            h) const
{	
	DKDa.clear();

	vector<TKalMatrix> DKDA;

	CalcDKDA(DKDA, K, h);

	TKalMatrix dada(3,3);
	TKalMatrix dkda(3,3);
	
	TKalMatrix um(3,3);
	um.UnitMatrix();

	dada.UnitMatrix();
	dkda = DKDA[0]*dada;
	DKDa.push_back(dkda);

	for(Int_t i=1; i<4; i++)
	{
		dada = um + hk[i] * h * dkda;
		dkda = DKDA[i] * dada;
		DKDa.push_back(dkda);
	}
}

void TRungeKuttaTrack::CalcDKDx(vector<TKalMatrix>& DKDx, 
		                        vector<TVector3>&   K, 
								Double_t            h) const
{	
	DKDx.clear();

	vector<TKalMatrix> DKDA;

	CalcDKDA(DKDA, K, h);

	TKalMatrix dadx(3,3);
	TKalMatrix dkdx(3,3);
	
	DKDx.push_back(dkdx);

	for(Int_t i=1; i<4; i++)
	{
		dadx = hk[i] * h * dkdx;
		dkdx = DKDA[i] * dadx;
		DKDx.push_back(dkdx);
	}
}

void TRungeKuttaTrack::CalcDKDh(vector<TKalMatrix>& DKDh, vector<TVector3>& K, Double_t h) const 
{
	DKDh.clear();

	vector<TKalMatrix> DKDA;

	CalcDKDA(DKDA, K, h);

	TKalMatrix dadh(3,1);
	TKalMatrix dkdh(3,1);
	
	dkdh.Zero();
	DKDh.push_back(dkdh);

	for(Int_t i=1; i<4; i++)
	{
		dadh = hk[i] * TKalMatrix::ToKalMat(K[i-1]) + hk[i] * h * dkdh;
		dkdh = DKDA[i] * dadh;
		DKDh.push_back(dkdh);
	}
}	



TKalMatrix TRungeKuttaTrack::CalcDxDp(Double_t h) const
{
	Double_t p = GetCurMomentum().Mag();

	vector<TKalMatrix> DKDa;
	vector<TVector3> K;//not used

	CalcDKDa(DKDa, K, h);
	
	TKalMatrix um(3,3);
	um.UnitMatrix();

	TKalMatrix dxdp = h/p * um + h * h /(6 * p) * (DKDa[0] + DKDa[1] + DKDa[2]);

	return dxdp;
}

void TRungeKuttaTrack::UpdatePX()
{
   //get the local track position
   Double_t csf0 = TMath::Cos(fPhi0);
   Double_t snf0 = TMath::Sin(fPhi0);
   Double_t x    = fX0.X() + fDrho * csf0;
   Double_t y    = fX0.Y() + fDrho * snf0;
   Double_t z    = fX0.Z() + fDz;         
   TVector3 xv(x,y,z);
   
   TVector3 pv;
   SVToP(pv);

   SetCurMomentum(pv);
   SetCurPosition(xv);
}

void TRungeKuttaTrack::StepRungeKutta(Double_t  step, 
					                  TVector3& vx, 
	             				      TVector3& vp) const
{
  //vx: vector of position
  //vp: vector of momentum

  const double khalf  = 1./2.;
  const double ksixth = 1./6.;

  //unit: mm, GeV/c and Tesla
  const double kec = 2.99792458e-4;

  double h, hhalf, hsixth;

  double tl     = 0.;
  double s      = step;
  double rest;

  double p     = vp.Mag();
  double kappa = kec * fCharge / p;

#ifdef __DEBUG__
  cout << "******************" << endl;
  cout << "Charge : " << fCharge << endl;
  cout << "p : " << p << endl;
  cout << "kec : " << kec << endl;
  cout << "******************" << endl;
#endif  

  //magnetic field
  TVector3 bf; 

  //fourth order Runge-Kutta 
  TVector3 K1, K2, K3, K4;

  //tangent of track
  TVector3 alpha;

  alpha = vp.Unit();

  if(fabs(s)<=1.e-10)
  { 
#ifdef __DEBUG__
	  cout << "stepping return" << endl;
#endif

	  return;
  }

#ifdef __DEBUG__
  cout << "vp org = " << endl;
  vp.Print();
  alpha.Print();
#endif

  do {
    rest  = step - tl;
    if (fabs(s) > fabs(rest)) s = rest;

	h        = s;
    hhalf    = khalf  * h;
	hsixth   = ksixth * h;

	bf = GetLocalBfield(vx);

	K1 = kappa * alpha.Cross(bf);

    bf = GetLocalBfield(vx + hhalf * (alpha + khalf * hhalf * K1));

	K2 = kappa * (alpha + hhalf * K1).Cross(bf);
	K3 = kappa * (alpha + hhalf * K2).Cross(bf);

    bf = GetLocalBfield(vx + h * (alpha + hhalf * K3));
	K4 = kappa * (alpha + h * K3).Cross(bf);

	vx += (alpha + (K1 + K2 + K3) * hsixth) * h;
	alpha += (K1 + K4 + 2. * (K2 + K3)) * hsixth;

#ifdef __DEBUG__
	TVector3 deltav = (K1 + K4 + 2. * (K2 + K3)) * hsixth;
	cout << "the delta vector is :" << endl;
	deltav.Print();

	cout << "the mag is : " << endl;
	bf.Print();
	cout << "kappa = " << kappa << endl;
	cout << "list of K:" << endl;
	K1.Print();
	K2.Print();
	K3.Print();
	K4.Print();
#endif

	//cout << "ALPHA" << endl;
    //alpha.Print();
	alpha = alpha.Unit();
    //alpha.Print();

    tl += s;
	rest = step - tl;

	//cout << "rest = " << rest << endl;

    if(step < 0.)rest = -rest;
    if(rest <= 1.e-5 * fabs(step))
    {
       vp = alpha * p;
       return;
    }
  }while(1);
}

void TRungeKuttaTrack::StepRungeKutta(Double_t step)
{
#ifdef __DEBUG__
	cout << "TRungeKuttaTrack::StepRungeKutta()" << endl;
#endif

	TVector3 xv = GetCurPosition();
	TVector3 xp = GetCurMomentum();

#ifdef __DEBUG__
	cout << "Step: " << step << endl;
	cout << "x and p are " << endl;
	xv.Print();
	xp.Print();
#endif

	StepRungeKutta(step, xv, xp);

	SetCurPosition(xv);
	SetCurMomentum(xp);

#ifdef __DEBUG__
	cout << "End of TRungeKuttaTrack::StepRungeKutta()" << endl;
#endif
}

void TRungeKuttaTrack::CalcXPAt(Double_t step, TVector3& vx, TVector3& vp) const
{
    vx = GetCurPosition();
	vp = GetCurMomentum();

	//Runge-Kutta method 
	StepRungeKutta(step, vx, vp);
}

void TRungeKuttaTrack::PXToSV(TVector3& p, TVector3& /* x */, TKalMatrix& sv)
{ 
	// This function is for the case of that pivot is on the track, 
	// so the first and the third element of state vector (i.e. drho and dz) are both zero.
	
	const Double_t PI2 = TMath::TwoPi();
    	
	sv(0,0) = 0.;
	sv(1,0) = atan2(-p.X(), p.Y());
	sv(2,0) = 1./p.Perp() * fSkappa;
	sv(3,0) = 0.;
	sv(4,0) = 1./p.Perp() * p.Z();
    
	while(sv(1,0)<0.0) sv(1,0) += PI2;
	while(sv(1,0)>PI2) sv(1,0) -= PI2;
}

void TRungeKuttaTrack::SVToP(TVector3& p)
{ 
	// we assume the absolute value of charge for a particle is |Q| = 1.
	
	fSkappa = TMath::Sign(1., fKappa);
	fCharge = fSkappa;

	Double_t cpa  = fKappa;
	Double_t phi0 = fPhi0;
	Double_t tanl = fTanL; 

	Double_t pt   = 1./fabs(cpa);

	Double_t sinphi0 = sin(phi0);
	Double_t cosphi0 = cos(phi0);
	
	p.SetXYZ(-pt*sinphi0, pt*cosphi0, pt*tanl);
}
