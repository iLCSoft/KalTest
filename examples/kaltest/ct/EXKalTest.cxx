#include "TNtupleD.h"
#include "TFile.h"
#include "TKalDetCradle.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TKalTrack.h"
#include "EXKalTest.h"
#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"
#include "EXToyGLD.h"

#include <iostream>

//static const Bool_t gkDir = kIterBackward;
static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:cpa");

   Double_t pt      = 1.;
   Int_t    nevents = 1;
   switch (argc) {
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
   // Prepare a detector
   // ===================================================================

   EXToyGLD toygld;
   EXKalDetector detector;

   toygld.Install(detector); // install detector into its toygld
#ifdef __MS_OFF__
   toygld.SwitchOffMS();     // switch off multiple scattering
#endif

   // ===================================================================
   // Prepare a Event Generator
   // ===================================================================

   EXEventGen gen(toygld);

   // ===================================================================
   // Event loop
   // ===================================================================

   for (Int_t eventno = 0; eventno < nevents; eventno++) { 
      cerr << "------ Event " << eventno << " ------" << endl;

      // ---------------------------
      //  Reset hit data
      // ---------------------------

      toygld.Reset();

      // ============================================================
      // Generate a partcle
      // ============================================================

      THelicalTrack hel = gen.GenerateHelix(pt);

      // ============================================================
      // Swim the particle in detector
      // ============================================================

      gen.Swim(hel);

      // ============================================================
      // Create Kalman Filter
      // ============================================================

      // ---------------------------
      //  Prepare hit iterrator
      // ---------------------------
      TObjArray &kalhits = toygld.GetHits();
      TIter next(&kalhits, gkDir);   // come in to IP

      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------

      EXHit &hitd = *dynamic_cast<EXHit *>(next());
      hitd(0,1) = 1.e6;   // give a huge error to d
      hitd(1,1) = 1.e6;   // give a huge error to z
      next.Reset();       // rewind iterator

      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      // sited.Lock();    // dummy site should not be used
      sited.SetOwner();   // site owns states

      // ---------------------------
      // Create initial helix
      // ---------------------------

      Int_t i1, i2, i3;
      if (gkDir == kIterBackward) {
         i3 = 0;
         i1 = toygld.GetHits().GetEntries() - 1;
         i2 = i1 / 2;
      } else {
         i1 = 0;
         i3 = toygld.GetHits().GetEntries() - 1;
         i2 = i3 / 2;
      }
      EXHit   &h1 = *dynamic_cast<EXHit *>(kalhits.At(i1));   // first hit
      EXHit   &h2 = *dynamic_cast<EXHit *>(kalhits.At(i2));   // last hit
      EXHit   &h3 = *dynamic_cast<EXHit *>(kalhits.At(i3));   // middle hit
      TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
      TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
      TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield());   // initial helix 

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
      kaltrack.Add(&sited);

      // ---------------------------
      //  Start Kalman Filter
      // ---------------------------

      EXHit *hitp = 0;
      while ((hitp = dynamic_cast<EXHit *>(next()))) {
         TKalTrackSite  &site = *new TKalTrackSite(*hitp);
         if (!kaltrack.AddAndFilter(site)) {
            cerr << " site discarded!" << endl;
            delete &site;
         }
      }
      //kaltrack.SmoothBackTo(3);

      // ============================================================
      // Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t cpa  = kaltrack.GetState(TVKalSite::kFiltered)(2, 0);
      hTrackMonitor->Fill(ndf, chi2, cl, cpa);
   }

   hfile.Write();

   return 0;
}
