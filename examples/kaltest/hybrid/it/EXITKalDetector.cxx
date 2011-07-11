
#include "EXITKalDetector.h"
#include "EXITMeasLayer.h"
#include "EXITFBMeasLayer.h"
#include "TRandom.h"
#include <sstream>

ClassImp(EXITKalDetector)

EXITKalDetector::EXITKalDetector(Int_t m)
               : EXVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("ITAir", "", A, Z, density, radlen, 0.);

   A       = 28.0855;
   Z       = 14.;
   density = 2.33;
   radlen  = 9.36;
   TMaterial &si = *new TMaterial("ITSi", "", A, Z, density, radlen, 0.);

   ////////// Barrel part of IT /////////////////////

   static const Int_t    nlayers  = 2;
   static const Double_t lhalfmin = 36.8;      // layer min half length [cm]
   static const Double_t rmin     = 17.92;     // layer radius min
   static const Double_t rstep    = 10.55;     // layer radius step
   static const Double_t lstep    = 31.;       // layer length step
   static const Double_t thick    = 0.05616;   // layer thick
   static const Double_t sigmax   = 1.e-3;
   static const Double_t sigmaz   = 1.e-3;
  
   ///////// Forward & Backward Part of IT //////////

   static const Double_t rinfb    = 1.6;      // dummy layer inner radius 
   
   InitFTDGeometry(); // get FTD geometry;
   
   Bool_t active = EXVMeasLayer::kActive;
   Bool_t dummy  = EXVMeasLayer::kDummy;

   // create measurement layers of inner tracker
   Double_t r     = rmin;
   Double_t len   = lhalfmin;
   Double_t lastz = 0; 

   for (Int_t layer = 0; layer < _nDiskHalf; layer++) {

     Double_t  npetal = _FTDgeo[layer].nPetal;
     Double_t       z = _FTDgeo[layer].zCoord;
     Double_t     rin = _FTDgeo[layer].rin;
     Double_t    rout = _FTDgeo[layer].rout;
     Double_t spolicy = _FTDgeo[layer].spolicy;
     Double_t    cosa = _FTDgeo[layer].cosAlpha;
     Double_t    sina = _FTDgeo[layer].sinAlpha;
     Double_t   dxMax = _FTDgeo[layer].dxMax;
     Double_t   dxMin = _FTDgeo[layer].dxMin;
     Double_t fbthick = _FTDgeo[layer].lyrthick;  
     
     
     static const  Double_t eps1 = 1e-6;
     static const  Double_t eps2 = 1e-4;
     
     for(Int_t pet=0; pet<npetal; pet++){

       std::stringstream ss;
       std::stringstream ssdm;
       
       Double_t cosphi = _FTDgeo[layer].cosphi[pet];
       Double_t sinphi = _FTDgeo[layer].sinphi[pet]; 
       
       TVector3 xc1(sinphi*(rout+rin)/2,cosphi*(rout+rin)/2,+z);         // for forward part  
       TVector3 xc2(xc1.X(),xc1.Y(),xc1.z()+fbthick); 
       TVector3 xc3(sinphi*(rout+rin)/2,cosphi*(rout+rin)/2,-z);         // for backward part 
       TVector3 xc4(xc3.X(),xc3.Y(),xc3.Z()-fbthick); 
       TVector3 normalF(cosphi*sina,-sinphi*sina,cosa);
       TVector3 normalB(-normalF);
       
       if(pet==0){           // Separate the first petal for numbering of overlapped part.
	 // Forward part
	 ss << "ITF" << layer << "_petalR" << std::ends;
	 ssdm << "ITFdm" << layer << "_petalR" << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc1, normalF, spolicy+2*pet*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, active, 1 ,ss.str().data()));	   
	 Add(new EXITFBMeasLayer(si, air, xc2, normalF, spolicy+(2*pet+1)*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,1,ssdm.str().data()));	   
	 ss.clear();
	 ssdm.clear();
	 ss << "ITF" << layer << "_petalL" << std::ends;
	 ssdm << "ITFdm" << layer << "_petalL" << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc1, normalF, spolicy+(2*npetal)*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, active, -1 ,ss.str().data()));	   
	 Add(new EXITFBMeasLayer(si, air, xc2, normalF, spolicy+(2*npetal+1)*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,-1,ssdm.str().data()));	   
	 
	 //  Backward part
	 ss.clear();
	 ssdm.clear();
	 ss << "ITB" << layer << "_petalR" << std::ends;
	 ssdm << "ITBdm" << layer << "_petalR" << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc3, normalB, spolicy+(2*npetal)*eps1+eps2, rin, rout,dxMax, dxMin, sigmax, sigmaz, active, 1 ,ss.str().data()));	   
	 Add(new EXITFBMeasLayer(si, air, xc4, normalB, spolicy+(2*npetal+1)*eps1+eps2, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,1,ssdm.str().data()));
	 ss.clear();
	 ssdm.clear();
	 ss << "ITB" << layer << "_petalL" << std::ends;
	 ssdm << "ITBdm" << layer << "_petalL" << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc3, normalB, spolicy+(2*pet)*eps1+eps2, rin, rout, dxMax, dxMin, sigmax, sigmaz, active, -1 ,ss.str().data()));	   
	 Add(new EXITFBMeasLayer(si, air, xc4, normalB, spolicy+(2*pet+1)*eps1+eps2, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,-1,ssdm.str().data()));
	 
       }else{
	 // Forward part
	 ss << "ITF" << layer << "_petal" << pet  << std::ends;
	 ssdm << "ITFdm" << layer << "_petal" << pet  << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc1, normalF, spolicy+2*pet*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, active,0,ss.str().data()));
	 Add(new EXITFBMeasLayer(si, air, xc2, normalF, spolicy+(2*pet+1)*eps1, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,0,ssdm.str().data()));
	 
	 //  Backward part
	 ss.clear();
	 ssdm.clear();
	 ss << "ITB" << layer << "_petal" << pet << std::ends;
	 ssdm << "ITBdm" << layer << "_petal" << pet << std::ends;
	 Add(new EXITFBMeasLayer(air, si, xc3, normalB, spolicy+2*(npetal-pet)*eps1+eps2, rin, rout, dxMax, dxMin, sigmax, sigmaz, active, 0 ,ss.str().data()));
	 Add(new EXITFBMeasLayer(si, air, xc4, normalB, spolicy+(2*(npetal-pet)+1)*eps1+eps2, rin, rout, dxMax, dxMin, sigmax, sigmaz, dummy,0,ssdm.str().data()));
       }
     }
     
     if(layer < _nDiskHalf && layer > 0){  // create completely dummy layers of FTD
       Double_t rifb = rinfb;
       Double_t rofb = rout;

       Double_t ndiv =5;
       Double_t interval = (z-lastz)/ndiv;
       
       for(Int_t div = 1; div<ndiv; div++){
	 for(Int_t dmpet=0; dmpet<npetal; dmpet++){
	   std::stringstream dm;
	   Double_t dmcosphi=_FTDgeo[0].cosphi[dmpet];
	   Double_t dmsinphi=_FTDgeo[0].sinphi[dmpet];
	   
	   Double_t dmz     = z - interval*div;
	   Double_t dmdxMax = 2*rofb*TMath::Tan(TMath::Pi()/npetal);
	   Double_t dmdxMin = 2*rifb*TMath::Tan(TMath::Pi()/npetal);
	   
	   TVector3 dmxc1(dmsinphi*(rofb+rifb)/2,dmcosphi*(rofb+rifb)/2,+dmz);
	   TVector3 dmxc3(dmsinphi*(rofb+rifb)/2,dmcosphi*(rofb+rifb)/2,-dmz);
	   TVector3 dmnormalF( 0, 0, 1);
	   TVector3 dmnormalB( 0, 0,-1);
	   
	   dm << "dmf" << layer << "_petal" << dmpet  << std::ends;
	   Add(new EXITFBMeasLayer(air, air, dmxc1, dmnormalF, spolicy-div*eps1-2*eps2, rifb, rofb, dmdxMax, dmdxMin, sigmax, sigmaz, dummy, 0.,dm.str().data()));       
	   dm.clear();
	   dm << "dmb" << layer << "_petal" << dmpet  << std::ends;
	   Add(new EXITFBMeasLayer(air, air, dmxc3, dmnormalB, spolicy-div*eps1-eps2, rifb, rofb, dmdxMax, dmdxMin, sigmax, sigmaz, dummy, 0.,dm.str().data()));	   
	 }
       }
     }
     
     
     if (layer < nlayers) { 
       std::stringstream ss;
       std::stringstream ssdm;
       ss   << "IT" << layer << std::ends;
       ssdm << "ITdm" << layer << std::ends;
       Add(new EXITMeasLayer(air, si, r, len, sigmax, sigmaz, active, ss.str().data()));
       Add(new EXITMeasLayer(si, air, r + thick, len, sigmax, sigmaz, dummy, ssdm.str().data()));

       len += lstep;
       r   += rstep;       
     }

     lastz = z;
   }
   
   SetOwner();
}

EXITKalDetector::~EXITKalDetector()
{
}
