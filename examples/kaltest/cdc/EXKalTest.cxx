//#define __DEBUG__
#define __MS_ON__
//#define __T0DET__
//#define __HELIXFIT__
//#define __BACKWARD__
//#define __AXIALONLY__

#include <sstream>
#include "EXKalTest.h"
#include "EXHit.h"
#include "EXMeasLayer.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TKalTrack.h"
#include "TKalDetCradle.h"
#include "EXKalDetector.h"

#include "THelicalTrack.h"
#include "THype.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"

#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"

//-----------------------------------
// Default # Tracks to Generate
//-----------------------------------

#define __NTRKS__ 1

//-----------------------------------
// Detector Parameters
//-----------------------------------

#define __BFIELD__ 30.

#define __RT0DET__ 40.
#define __SIGMAXT0__ 0.0010
#define __SIGMAZT0__ 0.0010

#define __NSTEPS__ 50
#define __NWIRES__ 5
#define __RMIN__   45.
#define __RMAX__   155.
#define __RSTEP__  (__RMAX__-__RMIN__)/(__NSTEPS__-1)
#define __CELLHW__ 2.5
#define __LHALF__  150.
#ifdef __AXIALONLY__
#define __TANA__   0.
#else
#define __TANA__   0.05
#endif

#define __VDRIFT__ 0.7e-3
#define __SIGMAX__ 0.0085
//#define __SIGMAZ__ 0.3
#define __SIGMAZ__ 5.

//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __FI0__    0.
#define __CPA__   -0.01
#define __DZ__     0.
#define __TNLMIN__ -0.7
#define __TNLMAX__ 0.7
#define __T0__     14.
#define __X0__     0.
#define __Y0__     0.
#define __Z0__     __Y0__*__TNLMIN__
#define __Z0RANGE__ 0.



using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   //
   // -------------------------------------------------------------------
   // Read in job parameters
   //    Usage:
   //            ./EXKalTest [ntrks [t0 [pt]]
   // -------------------------------------------------------------------
   //

   Int_t    ntrks;
   Double_t t0in;
   Double_t cpain;
   
   switch (argc) {
      case 4: {
                 stringstream stin0;
                 stin0 << argv[1] << ends;
                 stin0 >> ntrks;
                 stringstream stin1;
                 stin1 << argv[2] << ends;
                 stin1 >> t0in;
                 stringstream stin2;
                 stin2 << argv[3] << ends;
                 stin2 >> cpain;
                 cpain = 1/cpain;
                 break;
              }
      case 3: {
                 stringstream stin0;
                 stin0 << argv[1] << ends;
                 stin0 >> ntrks;
                 stringstream stin1;
                 stin1 << argv[2] << ends;
                 stin1 >> t0in;
                 cpain = __CPA__;
                 break;
              }
      case 2: {
                 stringstream stin0;
                 stin0 << argv[1] << ends;
                 stin0 >> ntrks;
                 t0in  = __T0__;
                 cpain = __CPA__;
                 break;
              }
      case 1: {
                 ntrks = __NTRKS__;
                 t0in  = __T0__;
                 cpain = __CPA__;
                 break;
              }
      default:
              cerr << "Too many command line arguments!" << endl;
   }

   cerr << "---------------------------------------------" << endl
        << " " << ntrks << " events will be processed.   " << endl
        << " Input: t0 = " << t0in << " kappa = " << cpain << endl
        << "---------------------------------------------" << endl;

   Double_t t0 = t0in;

   //
   // -------------------------------------------------------------------
   // Define Hists and Plots
   // -------------------------------------------------------------------
   //

#ifdef __MSON__
   Double_t acpa  = TMath::Abs(cpain);
   Double_t cpamx = (acpa > 1. ? 0.05 : 0.005);
#else
   Double_t cpamx = 0.005;
#endif
   TFile hfile("h.root","RECREATE","KalTest");
 
   TH1D *hNDFPtr   = new TH1D("hNDF","NDF",200,0.,200.);
   TH1D *hChi2Ptr  = new TH1D("hChi2","Chi2",200,0.,200.);
   TH1D *hChi2sPtr = new TH1D("hChi2s","Chi2",200,0.,20.);
   TH1D *hChi2fPtr = new TH1D("hChi2f","Chi2",200,0.,20.);
   TH1D *hCLPtr    = new TH1D("hCL","CL",200,0.,1.);
   TH1D *hT0Ptr;
   if (kSdim == 6) hT0Ptr  = new TH1D("hT0","t0",100,-10,30.);
   TH1D *hDcpaPtr  = new TH1D("hDcpa","dKappa",200,-cpamx,cpamx);
#ifdef __HELIXFIT__
   TH1D *hNDFhPtr  = new TH1D("hNDFh","NDF h",200,0.,200.);
   TH1D *hChi2hPtr = new TH1D("hChi2h","Chi2 h",200,0.,200.);
   TH1D *hChi2hhPtr = new TH1D("hChi2hh","Chi2 h",200,0.,20.);
   TH1D *hCLhPtr   = new TH1D("hCLh","CL h",200,0.,1.);
   TH1D *hT0hPtr;
   if (kSdim == 6) hT0hPtr = new TH1D("hT0h","t0 h",100,-10,30.);
   TH1D *hDcpahPtr = new TH1D("hDcpah","dKappa h",200,-cpamx,cpamx);
#endif

   //
   // -------------------------------------------------------------------
   // Prepare a detector
   // -------------------------------------------------------------------
   //

   Int_t    nlayers = __NSTEPS__;
   Int_t    nwires  = __NWIRES__;
   Double_t lhalf   = __LHALF__;
   Double_t celhw   = __CELLHW__;
   Double_t tana    = __TANA__;
   Double_t rmin    = __RMIN__;
   Double_t rstep   = __RSTEP__;
   Double_t rt0det  = __RT0DET__;

#ifndef __T0DET__
   EXKalDetector detector(nlayers,nwires,lhalf,celhw,tana,rmin,rstep);
#else
   EXKalDetector detector(nlayers,nwires,lhalf,celhw,tana,rmin,rstep,rt0det);
#endif
   TKalDetCradle cradle;
   cradle.Install(detector);
#ifndef __MS_ON__
   cradle.SwitchOffMS();	// switch off multiple scattering
#endif

   //
   // -------------------------------------------------------------------
   // Start event loop
   // -------------------------------------------------------------------
   //

   for (Int_t itrk=0; itrk<ntrks; itrk++) {

   cerr << " -------------------------------------------------- " << endl;
   cerr << " -------------- " << endl;
   cerr << " track: " << itrk << endl;
   cerr << " -------------- " << endl;
 
   // ---------------------------
   //  Create helix track
   // ---------------------------

   Double_t dr  = __DR__;
   Double_t fi0 = __FI0__ + 2*TMath::Pi()*(gRandom->Uniform()-0.5);
   Double_t cpa = cpain;
   Double_t dz  = __DZ__;
   Double_t tnl = gRandom->Uniform(__TNLMIN__, __TNLMAX__);
   Double_t x0  = __X0__;
   Double_t y0  = __Y0__;
   Double_t z0  = __Z0__;

   Double_t b   = __BFIELD__;

   THelicalTrack heltrk(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);
   THelicalTrack hel1st = heltrk;

   // ---------------------------
   //  Create hits
   // ---------------------------

   TObjArray   hits(nlayers);
   hits.SetOwner();
   EXMeasLayer *msPtr;
   EXMeasLayer *msPtrs = 0;
   Int_t        nhits  = 0;
   Int_t        lyrmx  = 0;
   TVector3     x0first;
   Double_t     dfi  = -rt0det*cpa*0.03e-2*b;
   Bool_t       is1stloop = kTRUE;
   
   Int_t dlyr = 1;
   for (Int_t lyr = 0; lyr < nlayers && lyr >= 0; lyr += dlyr) {
       msPtr = (EXMeasLayer *)detector.At(lyr);
       EXMeasLayer &ms = *msPtr; // measurement layer
       TVector3 xx;
       if (lyr && lyrmx) dfi  = 0.;
       if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)) {
           if (is1stloop && hits.GetEntries() > 3) {
               is1stloop = kFALSE;
               dlyr  = -1;
               lyrmx = lyr;
               continue;
           } else break;
       }
       if (msPtrs) {       // Multiple scattering
          static const Double_t kMpi  = 0.13957018;
          static const Double_t kMpi2 = kMpi*kMpi;

          Double_t kpa    = heltrk.GetKappa();
          Double_t tanl   = heltrk.GetTanLambda();
          Double_t cslinv = TMath::Sqrt(1 + tanl*tanl);
          Double_t path   = TMath::Abs(heltrk.GetRho()*dfi)*cslinv;
          Bool_t   dir    = ((xx * TKalMatrix::ToThreeVec(heltrk.CalcDxDphi(dfi)))
                            * TMath::Sign(-1., kpa)) > 0;
          Double_t x0inv  = 1. / ms.GetMaterial(dir).GetRadLength();
          Double_t xl     = path * x0inv;
          Double_t mom    = TMath::Abs(1/kpa) * cslinv;
          Double_t beta   = mom/TMath::Sqrt(mom*mom + kMpi2); // pion assumed
          static const Double_t kMS1 = 0.01316;
          static const Double_t kMS2 = 0.038;
          Double_t tmp    = 1. + kMS2 * TMath::Log(TMath::Max(1.e-4, xl));
          tmp /= (mom * beta);
          Double_t sgms   = kMS1 * TMath::Sqrt(xl) * tmp;
          Double_t sgphi  = sgms*cslinv;
          Double_t sgtnl  = sgms*cslinv*cslinv;
          Double_t delphi = gRandom->Gaus(0.,sgphi);
          Double_t deltnl = gRandom->Gaus(0.,sgtnl);

          dfi *= 0.5;
          TVector3 x0ms = heltrk.CalcXAt(dfi);
          heltrk.MoveTo(x0ms,dfi);     // M.S. at mid point

          heltrk.ScatterBy(delphi,deltnl);
          dfi  = 0.;
          if (!ms.CalcXingPointWith(heltrk,xx,dfi)) break;// recalc exact hit
       }
       heltrk.MoveTo(xx,dfi);	// move pivot to current hit
       if (!nhits) x0first = heltrk.GetPivot();

       nhits++; // hit found

       Int_t lr     = 0;
       Int_t cellno = -99999999;
       TKalMatrix h = ms.XvToMv(xx, lr, cellno);
       Double_t   x = h(0,0);	            // drift length
       Double_t   z = h(1,0);	            // drift length

       Double_t vdrift;
       Double_t rdet  = ms.GetWireEnd().Perp();
       if (rdet > rt0det) {
          vdrift = __VDRIFT__;              // shift due to T0
       } else {
          vdrift = 0.;
       }
       x += vdrift*t0 * lr;

       Double_t dx;
       Double_t dz;
       if (rdet > rt0det) {
          dx    = __SIGMAX__;
          dz    = __SIGMAZ__;
       } else {
          dx    = __SIGMAXT0__;
          dz    = __SIGMAZT0__;
       }
       x += gRandom->Gaus(0.,dx);
       z += gRandom->Gaus(0.,dz);

       Double_t meas [2];
       Double_t dmeas[2];
       meas [0] = x;
       meas [1] = z;
       dmeas[0] = dx;
       dmeas[1] = dz;
       
       hits.Add(new EXHit(ms,meas,dmeas,lr,cellno,vdrift,xx,b));
#ifdef __DEBUG__ 
       cerr << " r= "    << setw(5) << setprecision(4) << ms.GetWireEnd().Y()
            << " ta= "   << setw(5) << setprecision(2) << ms.GetWireDir().X()
            << " "
            << " cel= "  << setw(2) << cellno
            << " "
            << " lr= "   << setw(2) << lr
            << " "
            << " x= "    << setw(7) << setprecision(4) << x 
            << " dx= "   << setw(4) << setprecision(2) << dx
            << " "
            << " z= "    << setw(5) << setprecision(4) << z 
            << " dz= "   << setw(2) << setprecision(2) << dz 
            << endl;
       cerr << setprecision(7);
#endif
       msPtrs = msPtr;
   }

   THelicalTrack hellast = heltrk;
#ifdef __BACKWARD__
   hel1st.MoveTo(heltrk.GetPivot(),dfi);// move pivot to 1st hit
#else
   hel1st.MoveTo(x0first,dfi);// move pivot to 1st hit
#endif

   // ---------------------------
   //  Create Kalman Filter
   // ---------------------------

   TKalTrack kalsys;			// a track is a kal system
   kalsys.SetOwner();			// kalsys owns sites

   // ---------------------------
   //  Prepare hit iterrator
   // ---------------------------

#ifdef __BACKWARD__
   TIter next(&hits,kIterBackward);
#else
   TIter next(&hits);
#endif

   // ---------------------------
   //  Create a dummy site: sited
   // ---------------------------

   if (!hits.GetEntries()) continue;

   EXHit       hitd        = *(EXHit *)next();
               hitd(0,1)   = 1.e4;      // give a huge error to d
               hitd(1,1)   = 1.e4;      // give a huge error to z
                                                                                
   next.Reset();                        // rewind iterator
                                                                                
   TKalTrackSite  &sited       = *(new TKalTrackSite(hitd));
                                                                                
   // sited.Lock();                     // dummy site should not be used
   sited.SetOwner();                    // site owns states
                                                                                
   // ---------------------------
   //  Set dummy state to sited
   // ---------------------------
                                                                                
   static TKalMatrix svd(kSdim,1);
   
   Int_t i1 = 0;
   Int_t i3 = nhits-2;
   Int_t i2 = (i3-i1)/2;
   EXHit &h1 = *(EXHit *)hits.At(i1);
   EXHit &h2 = *(EXHit *)hits.At(i2);
   EXHit &h3 = *(EXHit *)hits.At(i3);
   TVector3 x1 = h1.TVTrackHit::GetMeasLayer().HitToXv(h1);
   TVector3 x2 = h2.TVTrackHit::GetMeasLayer().HitToXv(h2);
   TVector3 x3 = h3.TVTrackHit::GetMeasLayer().HitToXv(h3);
   THelicalTrack helstart(x1, x2, x3, b);

   svd(0,0) = 0.;
   svd(1,0) = helstart.GetPhi0();
   svd(2,0) = helstart.GetKappa();
   svd(3,0) = 0.;
   svd(4,0) = helstart.GetTanLambda();
#if 0
   cerr << "----------------------------" << endl;
   cerr << "nhits  = " << nhits << endl;;
   cerr << "phi0t  = " << hel1st.GetPhi0() << endl;;
   cerr << "phi0s  = " << helstart.GetPhi0() << endl;;
   cerr << "kappat = " << hel1st.GetKappa() << endl;;
   cerr << "kappas = " << helstart.GetKappa() << endl;;
   cerr << "tanlt  = " << hel1st.GetTanLambda() << endl;;
   cerr << "tanls  = " << helstart.GetTanLambda() << endl;;
   cerr << "----------------------------" << endl;
#endif
   if (kSdim == 6) svd(5,0) = 0.;
                                                                                
   static TKalMatrix C(kSdim,kSdim);
   for (Int_t i=0; i<kSdim; i++) {
      C(i,i) = 1.e4;			// dummy error matrix
   }
                                                                                
   sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
   sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));
                                                                                
   // ---------------------------
   //  Add sited to the system
   // ---------------------------
                                                                                
   kalsys.Add(&sited);
                                                                                
#ifdef __DEBUG__ 
   cerr << " --------------------------- " << endl;
   cerr << " site: d" << endl;
   sited.DebugPrint();
   kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
#endif


   // ---------------------------
   //  Start Kalman Filter
   // ---------------------------

   Int_t  loop   = 1;
   Int_t  nsites = 1;
   EXHit *hitPtr;
   while ((hitPtr = (EXHit *)next())) {
#ifdef __DEBUG__ 
cerr << " --------------------------- " << endl;
cerr << " site: " << loop << endl;
#endif
      EXHit      &hit  = *hitPtr;
      TKalTrackSite  &site = *(new TKalTrackSite(hit));
      if (kalsys.AddAndFilter(site)) {
#ifdef __DEBUG__ 
cerr << " predicted:" << endl;
         kalsys.GetState(TVKalSite::kPredicted).DebugPrint();
         site.DebugPrint();
cerr << " filtered:" << endl;
         kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
         Int_t     ndf = kalsys.GetNDF();
         Double_t chi2 = kalsys.GetChi2();
         cerr << " ndf = "  << ndf
              << " chi2 = " << chi2 << endl;
         cerr << " cl = " << TMath::Prob(chi2,ndf) << endl;
#endif
         nsites++;
      } else {
         cerr << " site discarded!" << endl;
         delete &site;
      }
      loop++;
   } // end of Kalman Filter


   // ---------------------------
   // Moniter Fit Result
   // ---------------------------

#ifdef __DEBUG__ 
   kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
#endif

   TVector3 x0last = ((TKalTrackSite *)&kalsys.GetCurSite())->GetPivot();
   Double_t fid;
   THelicalTrack trktru = heltrk;
   trktru.MoveTo(x0last,fid,0,0);
#ifdef __DEBUG__
   dr   = trktru.GetDrho();
   fi0  = trktru.GetPhi0();
   cpa  = trktru.GetKappa();
   dz   = trktru.GetDz();
   tnl  = trktru.GetTanLambda();
   cerr << " (dr,fi0,cpa,dz,tnl)_tru = ("  
        << dr << ", " << fi0 << ", " << cpa << ", "
        << dz << ", " << tnl << ")"  << endl
        << " t0  = "  << t0 << endl
        << " X0  = (" << x0last.X() << ", "
                      << x0last.Y() << ", "
                      << x0last.Z() << ")" << endl;
#endif
   Int_t    ndf  = kalsys.GetNDF();
   Double_t chi2 = kalsys.GetChi2();
   Double_t t0f;
   if (kSdim == 6) t0f = kalsys.GetState(TVKalSite::kFiltered)(5,0);
   Double_t cl   = TMath::Prob(chi2,ndf);
   Double_t cpaf = kalsys.GetState(TVKalSite::kFiltered)(2,0);
#ifdef __DEBUG__ 
   Double_t drf  = kalsys.GetState(TVKalSite::kFiltered)(0,0);
   Double_t fi0f = kalsys.GetState(TVKalSite::kFiltered)(1,0);
   Double_t dzf  = kalsys.GetState(TVKalSite::kFiltered)(3,0);
   Double_t tnlf = kalsys.GetState(TVKalSite::kFiltered)(4,0);
   cerr << " (dr,fi0,cpa,dz,tnl)_fit = ("  
        << drf << ", " << fi0f << ", " << cpaf << ", "
        << dzf << ", " << tnlf << ")"  << endl;
   if (kSdim == 6) {
      cerr
        << " t0  = "  << t0f << endl;
   }
   cerr << " ndf = "  << ndf << " chi2 = " << chi2 << endl;
   cerr << " cl  = "  << cl  << endl;
#endif
   
#if 1
   TKalTrackState &afil = *(TKalTrackState *)&kalsys.GetState(TKalTrackSite::kFiltered);
   TKalMatrix atrufil(6,1);
   trktru.PutInto(atrufil);
   atrufil(5,0) = t0in;
   TKalMatrix delafil  = afil - atrufil;
   TKalMatrix delafilt = TKalMatrix(TKalMatrix::kTransposed,delafil);
   TKalMatrix cfilinv  = TKalMatrix(TKalMatrix::kInverted,afil.GetCovMat());
   Double_t chi2fil    = (delafilt * cfilinv * delafil)(0,0);
#ifdef __DEBUG__
   afil.GetCovMat().DebugPrint("cfil = ");
#endif
#endif

#ifdef __HELIXFIT__
   // ---------------------------
   // Helix fit if requested
   // ---------------------------

   Int_t ndfhel;
   TKalTrackState ahel
     = TKalTrackState((*(TKalTrackSite *)kalsys.At(nsites-1)).GetCurState(), 
                   *(TKalTrackSite *)kalsys.At(nsites-1));
   TKalMatrix Chel(kSdim,kSdim);
   (*(TKalTrackSite *)kalsys.At(0)).Lock();

   Double_t chi2hel  = kalsys.FitToHelix(ahel,Chel,ndfhel);
   
   TKalMatrix delah  = ahel - atrufil;
   TKalMatrix delaht = TKalMatrix(TKalMatrix::kTransposed,delah);
   TKalMatrix chinv  = TKalMatrix(TKalMatrix::kInverted,Chel);
   Double_t chi2helh = (delaht * chinv * delah)(0,0);
   Double_t clhel   = TMath::Prob(chi2hel,ndfhel);
   Double_t cpafhel = ahel(2,0);
   Double_t t0fhel;
   if (kSdim == 6) t0fhel = ahel(5,0);
#ifdef __DEBUG__ 
   cerr << "chi2helh = " << chi2helh << endl;
   delah.DebugPrint("delah = ");
   Chel.DebugPrint("Chel = ");
   chinv.DebugPrint("chinv = ");

   Double_t drfhel  = ahel(0,0);
   Double_t fi0fhel = ahel(1,0);
   Double_t dzfhel  = ahel(3,0);
   Double_t tnlfhel = ahel(4,0);
   cerr << " (dr,fi0,cpa,dz,tnl)_hel = ("  
        << drfhel << ", " << fi0fhel << ", " << cpafhel << ", "
        << dzfhel << ", " << tnlfhel << ")"  << endl;
   if (kSdim == 6) {
      cerr
        << " t0  = "  << t0fhel << endl;
   }

   cerr << " ndfhel = "  << ndfhel << " chi2hel = " << chi2hel << endl;
   cerr << " clhel  = "  << clhel  << endl;
#endif
#endif

   // ---------------------------
   //  Smoothing
   // ---------------------------

   kalsys.SmoothBackTo(3);

#ifdef __DEBUG__
   TVector3 x03rd  = ((TKalTrackSite *)kalsys.At(3))->GetPivot();
   trktru.MoveTo(x03rd,fid,0,0);
   dr   = trktru.GetDrho();
   fi0  = trktru.GetPhi0();
   cpa  = trktru.GetKappa();
   dz   = trktru.GetDz();
   tnl  = trktru.GetTanLambda();
   cerr << " (dr,fi0,cpa,dz,tnl)_3rd = ("  
        << dr << ", " << fi0 << ", " << cpa << ", "
        << dz << ", " << tnl << ")"  << endl
        << " t0  = "  << t0 << endl
        << " X0  = (" << x03rd.X() << ", "
                      << x03rd.Y() << ", "
                      << x03rd.Z() << ")" << endl;

   cerr << "Smoothing-----------" << endl;
   kalsys.GetState(TKalTrackSite::kSmoothed).DebugPrint();
#endif

   TKalTrackState &asmo = *(TKalTrackState *)&kalsys.GetState(TKalTrackSite::kSmoothed);
   TKalMatrix      atru(6,1);
   heltrk.PutInto(atru);
   atru(5,0) = t0in;
   TKalMatrix dela  = asmo - atru;
   TKalMatrix delat = TKalMatrix(TKalMatrix::kTransposed,dela);
   TKalMatrix cinv  = TKalMatrix(TKalMatrix::kInverted,asmo.GetCovMat());
   Double_t chi2smo = (delat * cinv * dela)(0,0);

#ifdef __HELIXFIT__
   hNDFhPtr ->Fill(ndfhel,1.);
   hChi2hPtr->Fill(chi2hel,1.);
   hChi2hhPtr->Fill(chi2helh,1.);
   hCLhPtr  ->Fill(clhel,1.);
   if (kSdim == 6) hT0hPtr  ->Fill(t0fhel,1.);
   hDcpahPtr->Fill(cpafhel-cpa,1.);
#endif

   hNDFPtr ->Fill(ndf,1.);
   hChi2Ptr->Fill(chi2,1.);
   hChi2sPtr->Fill(chi2smo,1.);
   hChi2fPtr->Fill(chi2fil,1.);
   hCLPtr  ->Fill(cl,1.);
   if (kSdim == 6) hT0Ptr  ->Fill(t0f,1.);
   hDcpaPtr->Fill(cpaf-cpa,1.);
   }

   hfile.Write();

   return 0;
}
