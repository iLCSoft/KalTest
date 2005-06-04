//*************************************************************************
//* =====================
//*  TVKalDetector Class
//* =====================
//*
//* (Description)
//*   Base class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalDetector
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************

#include "TVKalDetector.h"
#include "TKalTrack.h"
#include "TKalTrackSite.h"

ClassImp(TVKalDetector)

void TVKalDetector::CalcEnergyLoss(const TMaterial      &mat,
                                         Double_t        df,
                                   const TKalTrackState &afrom,
                                         TKalTrackState &ato) const
{
   Double_t cpa    = afrom(2, 0);
   Double_t tnl    = afrom(4, 0); 
   Double_t tnl2   = tnl * tnl;
   Double_t tnl21  = 1. + tnl2;
   Double_t cslinv = TMath::Sqrt(1. + tnl2);
   Double_t mom2   = tnl21 / (cpa * cpa);
#if 0
   static const Double_t dedx = 1.5e-3;   // [GeV/g]   very temporary!!
#else
   // -----------------------------------------
   // Bethe-Bloch eq. (Physical Review D P195.)
   // -----------------------------------------
   static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
   static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
   static const Double_t kMpi = 0.13957018;      // pion mass [GeV]

   TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSite::GetKalSystemPtr());
   Double_t   mass = ktp ? ktp->GetMass() : kMpi;
   //Double_t   mass = ktp ? ktp->GetMass() : 105.6e-3;  // muon mass
   //Double_t   mass = ktp ? ktp->GetMass() : 938e-3;    // proton mass

   Double_t dnsty = mat.GetDensity();  // density
   Double_t A     = mat.GetA();        // atomic mass
   Double_t Z     = mat.GetZ();        // atomic number
   //Double_t I    = Z * 1.e-8;   // mean excitation energy [GeV]
   //Double_t I    = (2.4 +Z) * 1.e-8;   // mean excitation energy [GeV]
   Double_t I    = (9.76 * Z + 58.8 * TMath::Power(Z, -0.19)) * 1.e-9;
   Double_t hwp  = 28.816 * TMath::Sqrt(dnsty * Z/A) * 1.e-9;
   Double_t bg2  = mom2 / (mass * mass);
   Double_t gm2  = 1. + bg2;
   Double_t meM  = kMe / mass;
   Double_t x    = log10(TMath::Sqrt(bg2));
   Double_t C0   = - (2. * log(I/hwp) + 1.);
   Double_t a    = -C0/27.;
   Double_t del;
   if (x >= 3.)            del = 4.606 * x + C0;
   else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * TMath::Power(3.-x, 3.);
   else                    del = 0.;
   Double_t tmax = 2.*kMe*bg2 / (1. + meM*(2.*TMath::Sqrt(gm2) + meM)); 
   Double_t dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
                 - bg2/gm2 - del);
#endif
   Double_t path = static_cast<const TKalTrackSite &>(afrom.GetSite()).IsInB()
                 ? TMath::Abs(afrom.GetHelix().GetRho()*df)*cslinv
                 : TMath::Abs(df)*cslinv;
   Double_t edep = dedx * dnsty * path;
   Double_t cpaa = TMath::Sqrt(tnl21 / (mom2 + edep
                 * (edep + 2. * TMath::Sqrt(mom2 + mass * mass))));
   Double_t dcpa = TMath::Abs(cpa) - cpaa;
#if 0
   if (cpa > 0 && !ktp->GetFitDirection() || cpa < 0 && ktp->GetFitDirection()) {
      ato(2, 0) -= dcpa;
   } else {
      ato(2, 0) += dcpa;
   }
#else
   if (cpa > 0) ato(2, 0) -= dcpa;
   else         ato(2, 0) += dcpa;
#endif
}

TKalMatrix TVKalDetector::CalcSigmaMS0(const TMaterial      &mat, 
                                             Double_t        df, 
                                       const TKalTrackState &from) const
{
   Double_t cpa    = from(2, 0);
   Double_t tnl    = from(4, 0); 
   Double_t tnl2   = tnl * tnl;
   Double_t tnl21  = 1. + tnl2;
   Double_t cpatnl = cpa * tnl;
   Double_t cslinv = TMath::Sqrt(1. + tnl2);
   Double_t mom    = TMath::Abs(1. / cpa) * cslinv;

   static const Double_t kMpi = 0.13957018; // pion mass [GeV]
   TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSite::GetKalSystemPtr());
   Double_t   mass = ktp ? ktp->GetMass() : kMpi;
   //Double_t   mass = ktp ? ktp->GetMass() : 105.6e-3;
   //Double_t   mass = ktp ? ktp->GetMass() : 938e-3;
   Double_t   beta = mom / TMath::Sqrt(mom * mom + mass * mass);

   Double_t x0inv = 1. / mat.GetRadLength();  // radiation length inverse

   // *Calculate sigma_ms0 =============================================
   static const Double_t kMS1  = 0.0136;
   static const Double_t kMS12 = kMS1 * kMS1;
   static const Double_t kMS2  = 0.038;
   
   Double_t path = static_cast<const TKalTrackSite &>(from.GetSite()).IsInB()
                 ? TMath::Abs(from.GetHelix().GetRho()*df)*cslinv
                 : TMath::Abs(df)*cslinv;
   Double_t xl   = path * x0inv;
   // ------------------------------------------------------------------
   // Very Crude Treatment!!
   Double_t tmp = 1. + kMS2 * TMath::Log(TMath::Max(1.e-4, xl));
   tmp /= (mom * beta);
   Double_t sgms2 = kMS12 * xl * tmp * tmp;
   // ------------------------------------------------------------------

   Int_t p = from.GetDimension();
   TKalMatrix q(p,p);
   q(1,1) = sgms2 * tnl21;
   q(2,2) = sgms2 * cpatnl * cpatnl;
   q(2,4) = sgms2 * cpatnl * tnl21;
   q(4,2) = sgms2 * cpatnl * tnl21;
   q(4,4) = sgms2 * tnl21  * tnl21;

   return q;
}
