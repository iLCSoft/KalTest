#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib

#include "EXKalTest.h"
#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"

#include "TNtupleD.h"         // from ROOT
#include "TFile.h"            // from ROOT

#include <iostream>

static const Bool_t gkDir = kIterBackward;
//static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:cpa:t0");

   Double_t pt      = -100.;
   Double_t t0in    =   14.;
   Int_t    nevents = 1;
   switch (argc) {
      case 4: 
         nevents = atoi(argv[1]);
         pt      = atof(argv[2]);
         t0in    = atof(argv[3]);
         break;
      case 3: 
         nevents = atoi(argv[1]);
         pt      = atof(argv[2]);
         break;
      case 2:
         nevents = atoi(argv[1]);
         break;
      case 1:
         break;
      default:
         cerr << "Too many command line arguments!" << endl;
         abort();
   }

   // ===================================================================
   //  Prepare a detector
   // ===================================================================

   TObjArray     kalhits;    // hit buffer
   TKalDetCradle cradle;     // detctor system
   EXKalDetector detector;   // CDC detector

   cradle.Install(detector); // install detector into its cradle
#ifdef __MS_OFF__
   cradle.SwitchOffMS();     // switch off multiple scattering
#endif

   // ===================================================================
   //  Prepare a Event Generator
   // ===================================================================

   EXEventGen gen(cradle, kalhits);
   gen.SetT0(t0in);

   // ===================================================================
   //  Event loop
   // ===================================================================

   for (Int_t eventno = 0; eventno < nevents; eventno++) { 
      cerr << "------ Event " << eventno << " ------" << endl;

      // ---------------------------
      //  Reset hit data
      // ---------------------------

      kalhits.Delete();

      // ============================================================
      //  Generate a partcle
      // ============================================================

      THelicalTrack hel = gen.GenerateHelix(pt);

      // ============================================================
      //  Swim the particle in detector
      // ============================================================

      gen.Swim(hel);

      // ============================================================
      //  Do Kalman Filter
      // ============================================================

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

      EXHit  hitd = *dynamic_cast<EXHit *>(kalhits.At(i1));
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
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield(),gkDir); // initial helix 

      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------

      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0.;
      svd(1,0) = helstart.GetPhi0();
      svd(2,0) = helstart.GetKappa();
      svd(3,0) = 0.;
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

      TKalTrack kaltrack;    // a track is a kal system
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track

      // ---------------------------
      //  Prepare hit iterrator
      // ---------------------------

      TIter next(&kalhits, gkDir);     // come in to IP

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
#if 1
      TVKalSite::EStType stype = TVKalSite::kFiltered;
#else
      TVKalSite::EStType stype = TVKalSite::kSmoothed;
      kaltrack.SmoothBackTo(1);                          // smooth back.
#endif

      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t cpa  = kaltrack.GetState(stype)(2, 0);
      Double_t t0   = kaltrack.GetState(stype)(5, 0);
      hTrackMonitor->Fill(ndf, chi2, cl, cpa, t0);
   }

   hfile.Write();

   return 0;
}
