#include "TNtupleD.h"
#include "TFile.h"
#include "TKalDetCradle.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TKalTrack.h"
#include "EXKalTest.h"
#include "EXTPCKalDetector.h"
#include "EXITKalDetector.h"
#include "EXVTXKalDetector.h"
#include "EXVTXHit.h"
#include "EXITHit.h"
#include "EXITFBHit.h"
#include "EXTPCHit.h"
#include "EXEventGen.h"
#include "TVTrackHit.h"

#include <iostream>

//#define __MS_OFF__

static const Bool_t gkDir = kIterBackward;
//static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:fi0:cpa:cs:t0");

   Double_t pt      =  1.;
   Double_t t0in    = 14.;
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

   TObjArray kalhits;
   TKalDetCradle    toygld;
   EXVTXKalDetector vtxdet;
   EXITKalDetector  itdet;
   EXTPCKalDetector tpcdet;

   toygld.Install(vtxdet); // install vtx into its toygld
   toygld.Install(itdet);  // install it into its toygld
   toygld.Install(tpcdet); // install tpc into its toygld
#ifdef __MS_OFF__
   toygld.SwitchOffMS();     // switch off multiple scattering
#endif
   toygld.Sort();   // temporary treatment

   // ===================================================================
   //  Prepare a Event Generator
   // ===================================================================

   EXEventGen gen(toygld, kalhits);
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

      if (kalhits.GetEntries() < 2) {
         cerr << "<<<<<< Shortage of Hits! nhits = " 
              << kalhits.GetEntries() << " >>>>>>>" << endl;
         continue;
      }

      Int_t i1, i2, i3;
      if (gkDir == kIterBackward) {
         i3 = 1;
         i1 = kalhits.GetEntries() - 1;
         i2 = i1 / 2;
      } else {
         i1 = 1;
         i3 = kalhits.GetEntries() - 1;
         i2 = i3 / 2;
      }

      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------

      TVTrackHit *ht1p = dynamic_cast<TVTrackHit *>(kalhits.At(i1));
      TVTrackHit *htdp = 0;
      if (dynamic_cast<EXVTXHit *>(ht1p)) {
         htdp = new EXVTXHit(*dynamic_cast<EXVTXHit *>(ht1p));
      } else if (dynamic_cast<EXITHit *>(ht1p)) {
         htdp = new EXITHit(*dynamic_cast<EXITHit *>(ht1p));
      } else if (dynamic_cast<EXITFBHit *>(ht1p)) {
         htdp = new EXITFBHit(*dynamic_cast<EXITFBHit *>(ht1p));
      } else if (dynamic_cast<EXTPCHit *>(ht1p)) {
         htdp = new EXTPCHit(*dynamic_cast<EXTPCHit *>(ht1p));
      }
      TVTrackHit &hitd = *htdp;

      hitd(0,1) = 1.e6;   // give a huge error to d
      hitd(1,1) = 1.e6;   // give a huge error to z

      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      // sited.Lock();    // dummy site should not be used
      sited.SetHitOwner();// site owns hit
      sited.SetOwner();   // site owns states

      // ---------------------------
      //  Create initial helix
      // ---------------------------

      TVTrackHit &h1 = *dynamic_cast<TVTrackHit *>(kalhits.At(i1)); // first hit
      TVTrackHit &h2 = *dynamic_cast<TVTrackHit *>(kalhits.At(i2)); // middle hit
      TVTrackHit &h3 = *dynamic_cast<TVTrackHit *>(kalhits.At(i3)); // last hit
      TVector3    x1 = h1.GetMeasLayer().HitToXv(h1);
      TVector3    x2 = h2.GetMeasLayer().HitToXv(h2);
      TVector3    x3 = h3.GetMeasLayer().HitToXv(h3);
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

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
      //  Prepare hit iterrator
      // ---------------------------

      TIter next(&kalhits, gkDir);   // come in to IP

      // ---------------------------
      //  Start Kalman Filter
      // ---------------------------

      TVTrackHit *hitp = 0;
      while ((hitp = dynamic_cast<TVTrackHit *>(next()))) {
         TKalTrackSite  &site = *new TKalTrackSite(*hitp);
         if (!kaltrack.AddAndFilter(site)) {
            cerr << " site discarded!" << endl;
            delete &site;
         }
      }
      //kaltrack.SmoothBackTo(3);

      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t fi0  = kaltrack.GetState(TVKalSite::kFiltered)(1, 0);
      Double_t cpa  = kaltrack.GetState(TVKalSite::kFiltered)(2, 0);
      Double_t tnl  = kaltrack.GetState(TVKalSite::kFiltered)(4, 0);
      Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl);
      Double_t t0   = kaltrack.GetState(TVKalSite::kFiltered)(5, 0);
      hTrackMonitor->Fill(ndf, chi2, cl, fi0, cpa, cs, t0);
   }

   hfile.Write();

   return 0;
}
