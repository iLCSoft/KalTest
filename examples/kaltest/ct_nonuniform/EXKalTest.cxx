#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "EXHYBTrack.h"        // from KalTrackLib

#include "EXKalTest.h"
#include "EXKalDetector.h"
#include "EXRKEventGen.h"
#include "EXRungeKuttaEventGen.h"
#include "EXHit.h"

#include "TNtupleD.h"         // from ROOT
#include "TFile.h"            // from ROOT
#include "TRandom.h"
#include "TCanvas.h"
#include "TView.h"
#include "TRotMatrix.h"
#include "TString.h"

#include "TH1F.h"

#include "TBField.h"

#include <iostream>

//#define __CLOCK__

static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   /////////////////////////////////////////////
   // input parameters
   // -b: batch mode
   // nevent: event number
   // tp: total momentum
   // theta: theta angle
   // coeff: the coefficient for magnetic field
   // useRK: use Runge-Kutta track
   /////////////////////////////////////////////

   Int_t offset = 0;

   if (argc > 1 && TString(argv[1]) == "-b") 
   {
	   offset = 1;
	   gROOT->SetBatch();   // batch mode without event display
   }
   else
   {
	   gROOT->SetBatch(kFALSE);
   }

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "tant:tanl:ndf:chi2:cl:cpa:pinv");

#ifdef __CLOCK__   
   clock_t t1, t2, total=0;
#endif

   Int_t    nevents = 1;
   Double_t tp      = 10.;
   Double_t theta   = 0.3;
   Double_t coeff   = 3.0;
   Bool_t   useRK   = kFALSE;

   switch (argc-offset) {
      case 6: 
         nevents = atoi(argv[1+offset]);
         tp      = atof(argv[2+offset]);
         theta   = atof(argv[3+offset]);
         coeff   = atof(argv[4+offset]);
		 useRK   = atof(argv[5+offset]) == 1 ? kTRUE : kFALSE;
         break;

      case 5: 
         nevents = atoi(argv[1+offset]);
         tp      = atof(argv[2+offset]);
         theta   = atof(argv[3+offset]);
         coeff   = atof(argv[4+offset]);
         break;

      case 4: 
         nevents = atoi(argv[1+offset]);
         tp      = atof(argv[2+offset]);
         theta   = atof(argv[3+offset]);
         break;

      case 3:
         nevents = atoi(argv[1+offset]);
         tp      = atof(argv[2+offset]);
         break;

      case 2:
         nevents = atoi(argv[1+offset]);
         break;

      case 1:
		 break;

      default:
         cerr << "Too many command line arguments!" << argc << " " << offset << endl;
         abort();
   }

   // ===================================================================
   //  Create TApplication
   // ===================================================================
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   // ===================================================================
   //  Prepare a detector
   // ===================================================================

   TObjArray     kalhits;    // hit buffer
   TKalDetCradle cradle;     // detctor system
   EXKalDetector detector;   // CT detector

   cradle.Install(detector); // install detector into its cradle

   cradle.SwitchOnMS();       // switch on multiple scattering
   cradle.SwitchOnDEDX();     // switch on energy loss

   // ===================================================================
   //  Prepare a Event Generator
   // ===================================================================

   kalhits.SetOwner();

   //EXRKEventGen gen(cradle, kalhits);
   EXRungeKuttaEventGen gen(cradle, kalhits);

   // ===================================================================
   //  Event loop
   // ===================================================================

   Double_t phi  = gRandom->Uniform(0, 2*TMath::Pi());
   
   //There is -1 difference between the sign of Runge-Kutta generator
   //in ROOT and that of KalTest.
   Double_t chg = 1;
   
   for (Int_t eventno = 0; eventno < nevents; eventno++) { 
      cerr << "------ Event " << eventno << " ------" << endl;

      // ---------------------------
      //  Reset hit data
      // ---------------------------

      kalhits.Delete();

      // ============================================================
      //  Generate a partcle
      // ============================================================

      phi  = gRandom->Uniform(0, 2*TMath::Pi());
      theta = gRandom->Uniform(0, 0.5);
	  Double_t costh = cos(theta);

	  TVector3 xstart;

	  Double_t pt = tp * costh;
	  TVector3 p(pt*cos(phi), pt*sin(phi), tp*sqrt(1-costh*costh));
	  
      TRungeKuttaTrack hel = gen.GenerateRKTrack(chg, xstart, p);

      // ============================================================
      //  Swim the particle in detector
      // ============================================================
      
	  TBField::SetBfieldCoeff(coeff);
      TBField::SetUseUniformBfield(kFALSE);
	  TKalDetCradle::SetUseRungeKuttaTrack(useRK);
      
	  gen.Swim(hel);

      // ============================================================
      //  Do Kalman Filter
      // ============================================================

#ifdef __CLOCK__   
      t1=clock();
#endif

	  //TBField::SetUseUniformBfield(kTRUE);

	  if(kalhits.GetEntries()<3)
	  {
		  cout << "Shortage of Hits!" << endl;
		  continue;
	  }

      Int_t i1, i2, i3;
      if (gkDir == kIterBackward) {
         i3 = 0;
         i1 = kalhits.GetEntries() - 1;
         i2 = i1 / 2;
      } else {
         i1 = 0;
         i3 = kalhits.GetEntries() - 1;
         i2 = i3 / 2;
      }

      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------

      EXHit hitd = *dynamic_cast<EXHit *>(kalhits.At(i1));
      hitd(0,1) = 1.e6;   // give a huge error to d
      hitd(1,1) = 1.e6;   // give a huge error to z

      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetOwner();   // site owns states

      // ---------------------------
      // Create initial helix
      // ---------------------------

      EXHit   &h1 = *dynamic_cast<EXHit *>(kalhits.At(i1));   // first hit
      EXHit   &h2 = *dynamic_cast<EXHit *>(kalhits.At(i2));   // last hit
      EXHit   &h3 = *dynamic_cast<EXHit *>(kalhits.At(i3));   // middle hit
      TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
      TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
      TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------

      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0.;
      svd(1,0) = 3.;  //3.
      svd(2,0) = 1.1; //1.
      svd(3,0) = 0.;  //0.
      svd(4,0) = 0.;  //0.

      svd(1,0) = helstart.GetPhi0();
      svd(2,0) = helstart.GetKappa();
      svd(4,0) = helstart.GetTanLambda();
      if (kSdim == 6) svd(5,0) = 0.;

      static TKalMatrix C(kSdim,kSdim);
      for (Int_t i=0; i<kSdim; i++) {
         C(i,i) = 1.e4;   // dummy error matrix
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      EXHYBTrack kaltrack;    // a track is a kal system
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track

      // ---------------------------
      //  Prepare hit iterrator
      // ---------------------------

      TIter next(&kalhits, gkDir);   // come in to IP

      // ---------------------------
      //  Start Kalman Filter
      // ---------------------------
      EXHit *hitp = 0;

      while ((hitp = dynamic_cast<EXHit *>(next()))) {     // loop over hits
          TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
         if (!kaltrack.AddAndFilter(site)) {               // add and filter this site
            cerr << " site discarded!" << endl;           
            delete &site;                                  // delete this site, if failed
         }
      }

      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t cpa  = kaltrack.GetCurSite().GetCurState()(2, 0);
	  Double_t tanl = kaltrack.GetCurSite().GetCurState()(4, 0);
	  Double_t pinv    = cpa / sqrt(1 + tanl*tanl);
      hTrackMonitor->Fill(tan(theta), tanl, ndf, chi2, cl, cpa, pinv);

#ifdef __CLOCK__
      t2=clock();
#endif

#ifdef __CLOCK__
	  total += t2-t1;
#endif

      // ============================================================
      //  Very Primitive Event Display
      // ============================================================
      static TCanvas *cvp = 0;
      if (!gROOT->IsBatch()) 
	  {
         if (!cvp) 
		 {
            cvp = new TCanvas("OED", "Event Display", 400, 10, 610, 610);
         } 
		 else 
		 {
            cvp->cd();
            cvp->Clear();
         }

         TView   *vwp = TView::CreateView(1,0,0);

         vwp->SetRange(-1200.,-1200.,-1200.,1200.,1200.,1200.);
         Int_t ierr;
         vwp->SetView(10.,80.,80.,ierr);

         detector.Draw(40);
         kaltrack.Draw(kRed,"");

         cout << "Next? [yes/no/edit/quit] " << flush;
         static const Int_t kMaxLen = 1024;
         Char_t temp[kMaxLen];
         cin.getline(temp,kMaxLen);
         TString opts(temp);
         opts.ToLower();
         
         if (!opts.Length()) 
		 {
            continue;
         } 
		 else if (opts[0] == 'n' || opts[0] == 'q') 
		 {
            break;
         } 
		 else if (opts[0] == 'e') 
		 {
            cout << "Select \"Quit ROOT\" from \"File\" to display next" << endl;
            cout << "\"CNTRL+C\" to really quit" << endl;
            app.Run(kTRUE);
         }
      } // end of event display
   }

   hfile.Write();

#ifdef __CLOCK__
   float ttime = float(total)/CLOCKS_PER_SEC;
   cout << "Time consumption of filtering: " << ttime << endl;
#endif

   return 0;
}
