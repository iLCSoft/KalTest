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

TRungeKuttaTrack::TRungeKuttaTrack(Double_t q, TVector3 x, TVector3 p)
	     :fPosition(x),
		  fMomentum(p),
		  fCharge  (q)
{
}

TRungeKuttaTrack::TRungeKuttaTrack(const TRungeKuttaTrack&  orig)
	     :fPosition(orig.fPosition), 
		  fMomentum(orig.fMomentum),
		  fCharge  (orig.GetCharge())
{
}

void TRungeKuttaTrack::StepRungeKutta(Double_t  step, 
					                  TVector3& vx, 
	             				      TVector3& vp)
{
  //vx: vector of position
  //vp: vector of momentum

  const double khalf  = 1./2.;
  const double ksixth = 1./6.;

  //unit: cm,gev/c and kgauss
  const double kec = 2.9979251e-4;

  double h, hhalf, hsixth;

  double tl     = 0.;
  double s      = step;
  double rest;

  double p     = vp.Mag();
  double kappa = kec * fCharge / p;

  //magnetic field
  TVector3 bf; 

  //fourth order Runge-Kutta 
  TVector3 K1, K2, K3, K4;

  //tangent of track
  TVector3 alpha;

  alpha = vp.Unit();

  do {
    rest  = step - tl;
    if (fabs(s) > fabs(rest)) s = rest;

	h        = s;
    hhalf    = khalf  * h;
	hsixth   = ksixth * h;

	bf = TRKMagField::GetField(vx);

	K1 = kappa * alpha.Cross(bf);

    bf = TRKMagField::GetField(vx + hhalf * (alpha + khalf * hhalf * K1));

	K2 = kappa * (alpha + hhalf * K1).Cross(bf);
	K3 = kappa * (alpha + hhalf * K2).Cross(bf);

    bf = TRKMagField::GetField(vx + h * (alpha + hhalf * K3));
	K4 = kappa * (alpha + h * K3).Cross(bf);

	vx += (alpha + (K1 + K2 + K3) * hsixth) * h;
	alpha += (K1 + K4 + 2. * (K2 + K3)) * hsixth;
	alpha = alpha.Unit();

    tl += s;
	rest = step - tl;

    if(step < 0.)rest = -rest;
    if(rest < 1.e-5 * fabs(step))
    {
       vp = alpha * p;
       return;
    }
  }while(1);
}


#if 0
void TRungeKuttaTrack::StepRungeKutta(Double_t step,
                                         Double_t* vect, 
										 Double_t* vout)
{
  ///	*  Input parameters						 *
  ///	*	CHARGE    Particle charge				 *
  ///	*	STEP	  Step size					 *
  ///	*	VECT	  Initial co-ords,direction cosines,momentum	 *
  ///	*  Output parameters						 *
  ///	*	VOUT	  Output co-ords,direction cosines,momentum	 *
  ///	*  User routine called  					 *
  ///	*	CALL GUFLD(X,F) 					 *

  //const Int_t kix  = 0;
  //const Int_t kiy  = 1;
  //const Int_t kiz  = 2;
  //const Int_t kipx = 3;
  //const Int_t kipy = 4;
  //const Int_t kipz = 5;
  //Double_t g1, g2, g3, g4, g5, g6;
  //Double_t f1, f2, f3, f4;
  //Double_t hrho, tet, norm, hp, rho1, sint, cost;

  Double_t h2, h4, f[4];
  Double_t /* xyzt[3], */ a, b, c, ph,ph2;
  Double_t secxs[4],secys[4],seczs[4];
  //Double_t hxp[3];
  Double_t dxt, dyt, dzt;
  //Double_t ang2;

  //Double_t est;
  Double_t at, bt, ct, cba;

  Double_t x;
  Double_t y;
  Double_t z;

  Double_t xt;
  Double_t yt;
  Double_t zt;

  // const Int_t maxit = 1992;
  const Int_t maxit  = 500;
  //const Int_t maxcut = 11;

  //const Double_t hmin   = 1e-4; // !!! MT ADD,  should be member
  //const Double_t kdlt   = 1e-3; // !!! MT CHANGE from 1e-4, should be member
  //const Double_t kdlt32 = kdlt/32.;
  const Double_t kthird = 1./3.;
  const Double_t khalf  = 0.5;
  const Double_t kec    = 2.9979251e-3;

  //const Double_t kpisqua = 9.86960440109;

  // *.
  // *.    ------------------------------------------------------------------
  // *.
  // *             this constant is for units cm,gev/c and kgauss
  // *
  Int_t iter = 0;
  Int_t ncut = 0;
  for(Int_t j = 0; j < 7; j++)
    vout[j] = vect[j];

  Double_t pinv   = kec * fCharge / vect[6];
  Double_t tl     = 0.;
  Double_t h      = step;
  Double_t rest;

  TVector3 bf;

  do {
    rest  = step - tl;
    if (TMath::Abs(h) > TMath::Abs(rest))
       h = rest;

	bf = TRKMagField::GetField(TVector3(vect[0], vect[1], vect[2]));
    f[0] = bf.X();
    f[1] = bf.Y();
    f[2] = bf.Z();

    // * start of integration
    x      = vout[0];
    y      = vout[1];
    z      = vout[2];
    a      = vout[3];
    b      = vout[4];
    c      = vout[5];

    h2     = khalf * h;
    h4     = khalf * h2;
    ph     = pinv * h;
    ph2    = khalf * ph;
    secxs[0] = (b * f[2] - c * f[1]) * ph2;
    secys[0] = (c * f[0] - a * f[2]) * ph2;
    seczs[0] = (a * f[1] - b * f[0]) * ph2;
    //ang2 = (secxs[0]*secxs[0] + secys[0]*secys[0] + seczs[0]*seczs[0]);
    //if (ang2 > kpisqua) break;

    dxt    = h2 * a + h4 * secxs[0];
    dyt    = h2 * b + h4 * secys[0];
    dzt    = h2 * c + h4 * seczs[0];
    xt     = x + dxt;
    yt     = y + dyt;
    zt     = z + dzt;

#if 0
    // * second intermediate point
    est = TMath::Abs(dxt) + TMath::Abs(dyt) + TMath::Abs(dzt);
    if (est > h) {
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }
#endif

    // xyzt[0] = xt;
    // xyzt[1] = yt;
    // xyzt[2] = zt;

	bf = TRKMagField::GetField(TVector3(xt, yt, zt));
    f[0] = bf.X();
    f[1] = bf.Y();
    f[2] = bf.Z();

    at     = a + secxs[0];
    bt     = b + secys[0];
    ct     = c + seczs[0];

    secxs[1] = (bt * f[2] - ct * f[1]) * ph2;
    secys[1] = (ct * f[0] - at * f[2]) * ph2;
    seczs[1] = (at * f[1] - bt * f[0]) * ph2;
    at     = a + secxs[1];
    bt     = b + secys[1];
    ct     = c + seczs[1];
    secxs[2] = (bt * f[2] - ct * f[1]) * ph2;
    secys[2] = (ct * f[0] - at * f[2]) * ph2;
    seczs[2] = (at * f[1] - bt * f[0]) * ph2;
    dxt    = h * (a + secxs[2]);
    dyt    = h * (b + secys[2]);
    dzt    = h * (c + seczs[2]);
    xt     = x + dxt;
    yt     = y + dyt;
    zt     = z + dzt;
    at     = a + 2.*secxs[2];
    bt     = b + 2.*secys[2];
    ct     = c + 2.*seczs[2];

#if 0
    est = TMath::Abs(dxt)+TMath::Abs(dyt)+TMath::Abs(dzt);
    if (est > 2.*TMath::Abs(h)) {
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }
#endif

    // xyzt[0] = xt;
    // xyzt[1] = yt;
    // xyzt[2] = zt;

	bf = TRKMagField::GetField(TVector3(xt, yt, zt));
    f[0] = bf.X();
    f[1] = bf.Y();
    f[2] = bf.Z();

    z      = z + (c + (seczs[0] + seczs[1] + seczs[2]) * kthird) * h;
    y      = y + (b + (secys[0] + secys[1] + secys[2]) * kthird) * h;
    x      = x + (a + (secxs[0] + secxs[1] + secxs[2]) * kthird) * h;

    secxs[3] = (bt*f[2] - ct*f[1])* ph2;
    secys[3] = (ct*f[0] - at*f[2])* ph2;
    seczs[3] = (at*f[1] - bt*f[0])* ph2;
    a      = a+(secxs[0]+secxs[3]+2. * (secxs[1]+secxs[2])) * kthird;
    b      = b+(secys[0]+secys[3]+2. * (secys[1]+secys[2])) * kthird;
    c      = c+(seczs[0]+seczs[3]+2. * (seczs[1]+seczs[2])) * kthird;

#if 0
    est    = TMath::Abs(secxs[0]+secxs[3] - (secxs[1]+secxs[2]))
      + TMath::Abs(secys[0]+secys[3] - (secys[1]+secys[2]))
      + TMath::Abs(seczs[0]+seczs[3] - (seczs[1]+seczs[2]));

    if (est > kdlt && TMath::Abs(h) > hmin) {
		std::cout << "haha" << std::endl;
      if (ncut++ > maxcut) break;
      h *= khalf;
      continue;
    }
#endif

    ncut = 0;
    // * if too many iterations, go to helix
    if (iter++ > maxit) break;

    tl += h;
    //if (est < kdlt32) h *= 2.;
    cba    = 1./ TMath::Sqrt(a*a + b*b + c*c);
    vout[0] = x;
    vout[1] = y;
    vout[2] = z;
    vout[3] = cba*a;
    vout[4] = cba*b;
    vout[5] = cba*c;
    rest = step - tl;
    if (step < 0.) rest = -rest;
    if (rest < 1.e-5*TMath::Abs(step))
    {
       //Float_t dot = (vout[3]*vect[3] + vout[4]*vect[4] + vout[5]*vect[5]);
       //fH.fPhi += TMath::ACos(dot);
       return;
    }

  } while(1);

  // angle too big, use helix

#if 0
  std::cout << "angle too big, use helix" << std::endl;
  f1  = f[0];
  f2  = f[1];
  f3  = f[2];
  f4  = TMath::Sqrt(f1*f1+f2*f2+f3*f3);
  rho = -f4*pinv;
  tet = rho * step;

  hnorm = 1./f4;
  f1 = f1*hnorm;
  f2 = f2*hnorm;
  f3 = f3*hnorm;

  hxp[0] = f2*vect[kipz] - f3*vect[kipy];
  hxp[1] = f3*vect[kipx] - f1*vect[kipz];
  hxp[2] = f1*vect[kipy] - f2*vect[kipx];

  hp = f1*vect[kipx] + f2*vect[kipy] + f3*vect[kipz];

  rho1 = 1./rho;
  sint = TMath::Sin(tet);
  cost = 2.*TMath::Sin(khalf*tet)*TMath::Sin(khalf*tet);

  g1 = sint*rho1;
  g2 = cost*rho1;
  g3 = (tet-sint) * hp*rho1;
  g4 = -cost;
  g5 = sint;
  g6 = cost * hp;

  vout[kix] = vect[kix] + g1*vect[kipx] + g2*hxp[0] + g3*f1;
  vout[kiy] = vect[kiy] + g1*vect[kipy] + g2*hxp[1] + g3*f2;
  vout[kiz] = vect[kiz] + g1*vect[kipz] + g2*hxp[2] + g3*f3;

  vout[kipx] = vect[kipx] + g4*vect[kipx] + g5*hxp[0] + g6*f1;
  vout[kipy] = vect[kipy] + g4*vect[kipy] + g5*hxp[1] + g6*f2;
  vout[kipz] = vect[kipz] + g4*vect[kipz] + g5*hxp[2] + g6*f3;

  //fH.fPhi += tet;
#endif
}
#endif

#if 0
void TRungeKuttaTrack::StepRungeKutta(Double_t step)
{
    const Double_t MM2CM = 0.1;
	const Double_t CM2MM = 10;

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
	StepRungeKutta(aStep, vect, vout);
	
	fMomentum.SetXYZ(vout[3]*vout[6], vout[4]*vout[6], vout[5]*vout[6]);
	fPosition.SetXYZ(vout[0]*CM2MM, vout[1]*CM2MM, vout[2]*CM2MM);

	//std::cout << std::endl;
}
#else
void TRungeKuttaTrack::StepRungeKutta(Double_t step)
{
    const Double_t MM2CM = 1.;//0.1
	const Double_t CM2MM = 1.;//10

    TVector3 vx = MM2CM*fPosition;	
	TVector3 vp = fMomentum;

	//step length
    Double_t aStep =  step * MM2CM;

	//RK
	StepRungeKutta(aStep, vx, vp);
	
	fMomentum = vp;
	fPosition = vx*CM2MM;
}
#endif
