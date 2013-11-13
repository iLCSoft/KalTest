#include <iostream>
#include <iomanip>
#include "TString.h"

#include "TTrackFrame.h"
#include "TMath.h"
#include "J4Timer.h"

//#define __J4TIMER__

#ifdef __FRAME_TIMER1__
double TTrackFrame::fTimeCtor = 0.;
#endif
#ifdef __FRAME_TIMER2__
double TTrackFrame::fTimeVec  = 0.;
#endif
#ifdef __FRAME_TIMER3__
double TTrackFrame::fTimeSv   = 0.;
#endif

ClassImp(TTrackFrame)

//-------------------------------------------------------
// Ctors and Dtor
//-------------------------------------------------------

TTrackFrame::TTrackFrame(const TTrackFrame &orig)
	           :fRotMat(orig.fRotMat),
			    fShift(orig.fShift),
				fDeltaRotMat(orig.fDeltaRotMat),
				fDeltaShift(orig.fDeltaShift)
{
}

TTrackFrame::TTrackFrame(const TTrackFrame& lastFrame, 
		                       const TVector3& v, 
							   const TVector3& b)
               :fRotMat(lastFrame.fRotMat), fShift(lastFrame.fShift), fDeltaShift(v)

{
#ifdef __FRAME_TIMER1__   
	clock_t t1, t2;
	t1=clock();
#endif

#ifdef __J4TIMER__
    static int timerid = -1;
    J4Timer timer(timerid, "TTrackFrame", "TTrackFrame_Ctor");
    timer.Start();
#endif

#ifdef __J4TIMER2__
    static int timerid2 = -1;
    J4Timer timer2(timerid2, "TTrackFrame2", "TTrackFrame_Ctor");
    timer2.Start();
#endif
	//Calculate local B field
	TVector3 localBField = fRotMat * b;

	//Calculate local transformation matrix and shift vector
	double theta = localBField.Theta();
	double phi = localBField.Phi();

	fDeltaRotMat.RotateZ(-phi);
	fDeltaRotMat.RotateY(-theta);
	fDeltaRotMat.RotateZ(phi);

	//Calculate fShift and fRotMat for the this frame
	TRotation invRotMat = fRotMat.Inverse();
	fShift  += invRotMat * fDeltaShift;
	fRotMat = fDeltaRotMat * fRotMat;

#ifdef __J4TIMER2__
	timer2.Stop();
#endif

#ifdef __FRAME_TIMER1__   
    t2 = clock();
	fTimeCtor += Double_t(t2 - t1)/CLOCKS_PER_SEC;
#endif

#ifdef __J4TIMER__
	timer.Stop();
#endif
}

//transform vector
TVector3 TTrackFrame::Transform(const TVector3 &v, TRType transformType) const
{
#ifdef __FRAME_TIMER2__   
	clock_t t1, t2;
	t1=clock();
#endif

#ifdef __J4TIMER__
    static int timerid = -1;
    J4Timer timer(timerid, "TTrackFrame", "TransformVector");
    timer.Start();
#endif

#ifdef __J4TIMER2__
    static int timerid2 = -1;
    J4Timer timer2(timerid2, "TTrackFrame2", "TransformVector");
    timer2.Start();
#endif

    TRotation invRotMat = fRotMat.Inverse();
    TVector3  transVector;

	switch (transformType)
	{
	    //LocalToLocal means from last local frame to this local frame
        case kLocalToLocal: 
			 transVector = fDeltaRotMat * (v - fDeltaShift);
			 break;

        case kLocalToGlobal:
             transVector = invRotMat * v + fShift;
			 break;

        case kGlobalToLocal:
			 transVector = fRotMat * (v - fShift);
			break;

		default:
			break;
	}

#ifdef __J4TIMER2__
	timer2.Stop();
#endif

#ifdef __J4TIMER__
	timer.Stop();
#endif

#ifdef __FRAME_TIMER2__   
	t2 = clock();
	fTimeVec += Double_t(t2 - t1)/CLOCKS_PER_SEC;
#endif

	return transVector;
}

//transform vector
TVector3 TTrackFrame::TransformBfield(const TVector3 &b, TRType transformType)
{
#ifdef __J4TIMER__
   static int timerid = -1;
   J4Timer timer(timerid, "TTrackFrame", "TransformBfield");
   timer.Start();
#endif
	
#ifdef __J4TIMER2__
   static int timerid2 = -1;
   J4Timer timer2(timerid2, "TTrackFrame2", "TransformBfield");
   timer2.Start();
#endif
    TRotation invRotMat = fRotMat.Inverse();
    TVector3  transVector;

	switch (transformType)
	{
	    //LocalToLocal means from last local frame to this local frame
        case kLocalToLocal: 
			 transVector = fDeltaRotMat * b;
			 break;

        case kLocalToGlobal:
             transVector = invRotMat * b;
			 break;

        case kGlobalToLocal:
			 transVector = fRotMat * b;
			break;

		default:
			break;
	}


#ifdef __J4TIMER2__
	timer2.Stop();
#endif

#ifdef __J4TIMER__
	timer.Stop();
#endif

	return transVector;
}

//transform matrix
void TTrackFrame::Transform(TKalMatrix* sv, TKalMatrix* Fr)
{
#ifdef __FRAME_TIMER3__   
	clock_t t1, t2;
	t1=clock();
#endif

#ifdef __J4TIMER__
   static int timerid = -1;
   J4Timer timer(timerid, "TTrackFrame", "Transform_sv");
   timer.Start();
#endif

#ifdef __J4TIMER2__
   static int timerid2 = -1;
   J4Timer timer2(timerid2, "TTrackFrame2", "Transform_sv");
   timer2.Start();
#endif

//#define __DEBUG__
#ifdef __DEBUG__
	std::cout << "TTrackFrame::Transform:" << std::endl;
#endif

	TKalMatrix localRot(fDeltaRotMat);

	//charge of particle
	Double_t chg = 1.;

	Double_t drhosign = TMath::Sign(1., (*sv)(0,0));
	Double_t cpasign  = TMath::Sign(1., (*sv)(2,0));

	Double_t drho  = (*sv)(0,0);
	Double_t phi0  = (*sv)(1,0);
	Double_t cpa   = (*sv)(2,0);
	Double_t dz0   = (*sv)(3,0);
	Double_t tanl  = (*sv)(4,0);

    const Double_t zeroDistLimit = 1.e-5;
    const Double_t PI2 = 2*TMath::Pi();

	Bool_t onPivot = false;

    if(fabs(drho)<zeroDistLimit && fabs(dz0)<zeroDistLimit)
	{
		onPivot = true;
	}
	
    TKalMatrix F(5,5);  //propagator from rotation
    
    if(onPivot==false)   //The helix is not on the pivot
    {
        //tansverse momentum
        Double_t pt = chg/fabs(cpa);
#ifdef __DEBUG__
		std::cout << "not on pivot, Pt = " << pt << std::endl;
#endif
      
        Double_t sinphi0 = sin(phi0);
        Double_t cosphi0 = cos(phi0);
      
		//momentum and positon vector
        TKalMatrix t(6,1);
        t(0,0) = -pt*sinphi0;
        t(1,0) =  pt*cosphi0;
        t(2,0) =  pt*tanl;
        t(3,0) = drho*cosphi0;
        t(4,0) = drho*sinphi0;
        t(5,0) = dz0;
  
        //DtDap
        TKalMatrix DtDap(6,5);
        CalcDtDap(DtDap, *sv, cpasign, 1);
      
        //calculate the roatated vector;
        TKalMatrix DtrDt(6,6);

        CalcDtrDt(DtrDt, localRot);
        TKalMatrix tr = DtrDt * t;
      
        Double_t px = tr(0,0);
        Double_t py = tr(1,0);
        Double_t pz = tr(2,0);
        Double_t dx = tr(3,0);
        Double_t dy = tr(4,0);
        Double_t dz = tr(5,0);
      
		//Momentum and position in new frame
        TVector3 P(px,py,pz);
        TVector3 dX(dx,dy,dz);
      
		//
        //Calcualte state vector in new frame
		//

		//drho
        (*sv)(0,0) = dX.Perp() * drhosign;
         
		//phi0
        //If ..., then change the phi0 
        if(drho!=0.&&dx!=0.&&dy!=0.)
		{ 
			(*sv)(1,0) = atan2(drhosign*dy, drhosign*dx);
		}

        (*sv)(2,0) = 1./P.Perp() * chg * cpasign;
        (*sv)(3,0) = dX.Z();
        (*sv)(4,0) = 1./P.Perp() * P.Z();
      
        while((*sv)(1,0)<0.0) (*sv)(1,0) += PI2;
        while((*sv)(1,0)>PI2) (*sv)(1,0) -= PI2;
    
        //Get DappDtr
        TKalMatrix DappDtr(5,6);
        CalcDappDtr(DappDtr, tr, drhosign, cpasign, 1);
    
        F = DappDtr * DtrDt * DtDap; 
	}
    else //The helix is on the pivot
    {
		//Just get momentum.
        TKalMatrix t(3,1); 

        Double_t pt = chg/fabs(cpa);
    
#ifdef __DEBUG__
		std::cout << "helix is on pivot, input Pt = " << pt << std::endl;
#endif

        Double_t sinphi0 = sin(phi0);
        Double_t cosphi0 = cos(phi0);
    
        t(0,0) = -pt*sinphi0;
        t(1,0) =  pt*cosphi0;
        t(2,0) =  pt*tanl;
    
        TKalMatrix DtDap(3,3);
        CalcDtDap(DtDap, (*sv), cpasign);
    
        //Calculate momentum in new frame
        TKalMatrix DtrDt(3,3);
		
		//localRot.DebugPrint(" rotMat = ", 3);
        TKalMatrix tr = localRot * t;
        TVector3 P = TKalMatrix::ToThreeVec(tr);
    
        //Calculate new state vector
        (*sv)(1,0) = atan2(-cpasign * P.Y(), -cpasign * P.X()) + cpasign*TMath::Pi()/2;
        //(*sv)(1,0) = atan2(-P.X(), P.Y());
        (*sv)(2,0) = 1./P.Perp() * chg * cpasign;
        (*sv)(4,0) = 1./P.Perp() * P.Z();
    
        while((*sv)(1,0)<0.0) (*sv)(1,0) += PI2;
        while((*sv)(1,0)>PI2) (*sv)(1,0) -= PI2;
    
#ifdef __DEBUG__
		std::cout << "output Pt = " << P.Perp() << std::endl;
#endif
        TKalMatrix DappDtr(3,3);
        CalcDappDtr(DappDtr, tr, drhosign, cpasign);
    
        TKalMatrix temp = DappDtr * localRot * DtDap; 
        //F.Zero();
    
        F(0,0) = 1; //
         
        F(1,1) = temp(0,0);
        F(1,2) = temp(0,1);
        F(1,4) = temp(0,2);
    
        F(2,1) = temp(1,0);
        F(2,2) = temp(1,1);
        F(2,4) = temp(1,2);
    
        F(3,3) = 1; //
    
        F(4,1) = temp(2,0);
        F(4,2) = temp(2,1);
        F(4,4) = temp(2,2); 
	} 
	
	//To eliminate the divergence of matrix elements
	//if(drho==0.&&(*sv)(0,0)==0.&&dz0!=0.)
	//FIXME
	if(drho<zeroDistLimit&&(*sv)(0,0)==0.&&fabs(dz0)>zeroDistLimit)
    { 
#if 1
		std::cout << "------Matrix rotation, UnitMatrix" << std::endl;
#endif
		F.UnitMatrix();
	}

	//F.Print();

    TKalMatrix Ft(TKalMatrix::kTransposed, F);
    
	if(Fr!=0) 
	{
		F.ResizeTo(*Fr); 
		Ft.ResizeTo(*Fr);

		*Fr = F;
	}


#ifdef __J4TIMER2__
    timer2.Stop();
#endif

#ifdef __FRAME_TIMER3__   
	t2 = clock();
	fTimeSv += Double_t(t2 - t1)/CLOCKS_PER_SEC;
#endif

#ifdef __J4TIMER__
    timer.Stop();
#endif

#ifdef __DEBUG__
		std::cout << "Transformation completed" << std::endl;
#endif
}

void TTrackFrame::CalcDtrDt(TKalMatrix& dtrdt, const TKalMatrix& R)
{
   dtrdt.Zero();

   dtrdt(0,0) = dtrdt(3,3) = R(0,0);
   dtrdt(0,1) = dtrdt(3,4) = R(0,1);
   dtrdt(0,2) = dtrdt(3,5) = R(0,2);

   dtrdt(1,0) = dtrdt(4,3) = R(1,0);
   dtrdt(1,1) = dtrdt(4,4) = R(1,1);
   dtrdt(1,2) = dtrdt(4,5) = R(1,2);

   dtrdt(2,0) = dtrdt(5,3) = R(2,0);
   dtrdt(2,1) = dtrdt(5,4) = R(2,1);
   dtrdt(2,2) = dtrdt(5,5) = R(2,2);
}

void TTrackFrame::CalcDappDtr(TKalMatrix& dappdtr, TKalMatrix& tr, Double_t drhosign, 
                                Double_t cpasign,    Int_t  mode)
{ 
   Double_t chg  = 1.;   
   Double_t px, py, pz, dx, dy;

   px = tr(0,0);
   py = tr(1,0);
   pz = tr(2,0);

   if(mode)
   {
     dx = tr(3,0);
     dy = tr(4,0);
   }

   //
   Double_t dr2, dr, pt, pt2, pt3;

   if(mode)
   {
     dr2 = dx*dx + dy*dy;
     dr  = sqrt(dr2);
   }

   pt  = sqrt(px*px + py*py);
   pt2 = pt*pt;
   pt3 = pt*pt2;

   dappdtr.Zero();

   if(mode)
   {
     if(dr!=0.)
     {
       dappdtr(0,3) = drhosign*dx/dr;
       dappdtr(0,4) = drhosign*dy/dr;
     }

     if(dr2!=0.)
     {
       dappdtr(1,3) = -dy/dr2;
       dappdtr(1,4) =  dx/dr2;
     }
     
     dappdtr(2,0) = -cpasign*chg*px/pt3;
     dappdtr(2,1) = -cpasign*chg*py/pt3;
     
     dappdtr(3,5) = 1.;
     
     dappdtr(4,0) = -px*pz/pt3;
     dappdtr(4,1) = -py*pz/pt3;
     dappdtr(4,2) = 1./pt;
   }
   else
   {
     dappdtr(0,0) = -py/pt2;
     dappdtr(0,1) =  px/pt2;

     dappdtr(1,0) = -cpasign*chg*px/pt3;
     dappdtr(1,1) = -cpasign*chg*py/pt3;

     dappdtr(2,0) = -px*pz/pt3;
     dappdtr(2,1) = -py*pz/pt3;
     dappdtr(2,2) = 1./pt;
   }
}

void TTrackFrame::CalcDtDap(TKalMatrix& dtdap, TKalMatrix& a, Double_t cpasign, Int_t mode)
{
	//mode 0: helix is on its pivot

    Double_t chg     = 1.;    //the absolute value of particle charge
 
    Double_t dr, phi0, cpa, tanl;
    dr = phi0 = cpa = tanl = 0.;
 
    if(mode) dr = a(0,0);
    phi0  = a(1,0);
    cpa   = a(2,0);
    tanl  = a(4,0);
 
    Double_t sinphi0 = sin(phi0);
    Double_t cosphi0 = cos(phi0);
    Double_t cpa2    = cpa*cpa;
    Double_t acpa    = fabs(cpa);
 
    dtdap.Zero();
 
    if(mode)
    { 
		dtdap(0,1) = -chg/acpa*cosphi0;
        dtdap(0,2) = chg*cpasign/cpa2*sinphi0;
   
        dtdap(1,1) = -chg/acpa*sinphi0;
        dtdap(1,2) = -chg*cpasign/cpa2*cosphi0;
   
        dtdap(2,2) = -chg*cpasign/cpa2*tanl;
        dtdap(2,4) = chg/acpa;
   
        dtdap(3,0) = cosphi0;
        dtdap(3,1) = -dr*sinphi0;

        dtdap(4,0) = sinphi0;
        dtdap(4,1) = dr*cosphi0; 
		
		dtdap(5,3) = 1;
   }
   else
   {
        dtdap(0,0) = -chg/acpa*cosphi0;
        dtdap(0,1) = chg*cpasign/cpa2*sinphi0;
   
        dtdap(1,0) = -chg/acpa*sinphi0;
        dtdap(1,1) = -chg*cpasign/cpa2*cosphi0;
   
        dtdap(2,1) = -chg*cpasign/cpa2*tanl;
        dtdap(2,2) = chg/acpa; 
   }
}
