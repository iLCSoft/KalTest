//#define __MS_OFF__
//#define __DEDX_OFF__

#include "TNtupleD.h"
#include "TFile.h"
#include "TKalDetCradle.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TVTrackHit.h"
#include "EXKalTest.h"
#include "EXTPCKalDetector.h"
#include "EXITKalDetector.h"
#include "EXVTXKalDetector.h"
#include "EXVTXHit.h"
#include "EXITHit.h"
#include "EXITFBHit.h"
#include "EXTPCHit.h"
#include "EXEventGen.h"
#include "EXHYBTrack.h"

#include "TCanvas.h"
#include "TView.h"
#include "TRotMatrix.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TString.h"

#include <iostream>

static const Bool_t gkDir = kIterBackward;
//static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   // ===================================================================
   //  Get job parameters from command line arguments, if any
   // ===================================================================

   Int_t    offset  = 0;
   if (argc > 1 && TString(argv[1]) == "-b") {
      offset = 1;
      gROOT->SetBatch();   // batch mode without event display
   }

   Double_t pt      =  1.; // default Pt [GeV]
   Double_t t0in    = 14.; // default tp [nsec]
   Int_t    nevents = 1;   // default number of events to generate
   switch (argc-offset) {
      case 4: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         break;
      case 3: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         break;
      case 2:
         nevents = atoi(argv[1+offset]);
         break;
      case 1:
         break;
      default:
         cerr << "Too many command line arguments!" << endl;
         abort();
   }

   // ===================================================================
   //  Create TApplication and an ouput ntuple file
   // ===================================================================

   static TCanvas     *cvp    = 0;
   static TTUBE       *itubep = 0;
   static TTUBE       *otubep = 0;
   static TRotMatrix  *orotp  = 0;
   static TNode       *onodep = 0;
   static TNode       *inodep = 0;

   TApplication app("EXKalTest", &argc, argv, 0, 0);

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:fi0:cpa:cs:t0");

   // ===================================================================
   //  Prepare a detector
   // ===================================================================

   TKalDetCradle    toygld; // toy GLD detector
   EXVTXKalDetector vtxdet; // vertex detector (vtx)
   EXITKalDetector  itdet;  // intermediate tracker (it)
   EXTPCKalDetector tpcdet; // TPC (tpc)

   toygld.Install(vtxdet);  // install vtx into its toygld
   toygld.Install(itdet);   // install it into its toygld
   toygld.Install(tpcdet);  // install tpc into its toygld
   toygld.Sort();           // sort meas. layers from inside to outside

#ifdef __MS_OFF__
   toygld.SwitchOffMS();    // switch off multiple scattering
#endif
#ifdef __DEDX_OFF__
   toygld.SwitchOffDEDX();  // switch off enery loss
#endif

   // ===================================================================
   //  Prepare an Event Generator
   // ===================================================================

   TObjArray kalhits;                // array to store hits
   EXEventGen gen(toygld, kalhits);  // create event generator
   gen.SetT0(t0in);                  // set bunch crossing timing (t0)

   // ===================================================================
   //  Loop over events
   // ===================================================================

   for (Int_t eventno = 0; eventno < nevents; eventno++) { 
      cerr << "------ Event " << eventno << " ------" << endl;

      kalhits.Delete(); // clear hit data

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

      Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
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
      svd(0,0) = 0.;                        // dr
      svd(1,0) = helstart.GetPhi0();        // phi0
      svd(2,0) = helstart.GetKappa();       // kappa
      svd(3,0) = 0.;                        // dz
      svd(4,0) = helstart.GetTanLambda();   // tan(lambda)
      if (kSdim == 6) svd(5,0) = 0.;        // t0

      static TKalMatrix C(kSdim,kSdim);
      for (Int_t i=0; i<kSdim; i++) {
         C(i,i) = 1.e4;   // dummy error matrix
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      EXHYBTrack kaltrack;   // a track is a kal system
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to this track

      // ---------------------------
      //  Prepare hit iterrator
      // ---------------------------

      TIter next(&kalhits, gkDir); // come in to IP, if gkDir = kIterBackward

      // ---------------------------
      //  Start Kalman Filter
      // ---------------------------

      TVTrackHit *hitp = 0;
      while ((hitp = dynamic_cast<TVTrackHit *>(next()))) {
         TKalTrackSite  &site = *new TKalTrackSite(*hitp); // new site
         if (!kaltrack.AddAndFilter(site)) {               // filter it
            cerr << " site discarded!" << endl;
            delete &site;                        // delete it if failed
         }
      } // end of Kalman filter

      // ---------------------------
      //  Smooth the track
      // ---------------------------
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

      // ============================================================
      //  Very Primitive Event Display
      // ============================================================

      if (!gROOT->IsBatch()) {
         if (!cvp) {
            cvp = new TCanvas("OED", "Event Display", 400, 10, 610, 610);
         } else {
            cvp->cd();
            cvp->Clear();
         }

         TView   *vwp = new TView(1);
         vwp->SetRange(-260.,-260.,-260.,+260.,+260.,+260.);
         Int_t ierr;
         vwp->SetView(10.,80.,80.,ierr);

         if (!otubep) otubep = new TTUBE("TPCO","TPCO","void",210.,210.,260.);
         if (!itubep) itubep = new TTUBE("TPCI","TPCI","void", 44., 44.,255.);
         if (!orotp)  orotp  = new TRotMatrix("orot","orot",10.,80.,10.,80.,10.,80.);
         if (!onodep) {
            onodep = new TNode("ONODE","ONODE","TPCO",0.,0.,0.,"orot");
            onodep->cd();
            inodep = new TNode("INODE","INODE","TPCI");
         }
         onodep->Draw("pad same");
         kaltrack.Draw(2,"");         

         cout << "Next? [yes/no/edit/quit] " << flush;
         static const Int_t kMaxLen = 1024;
         Char_t temp[kMaxLen];
         cin.getline(temp,kMaxLen);
         TString opts(temp);
         opts.ToLower();
         if (!opts.Length()) {
            continue;
         } else if (opts[0] == 'n' || opts[0] == 'q') {
            break;
         } else if (opts[0] == 'e') {
            cout << "Select \"Quit ROOT\" from \"File\" to display next" << endl;
            cout << "\"CNTRL+C\" to really quit" << endl;
            app.Run(kTRUE);
         }
      } // endo fo event display
   } // end of event loop

   // ===================================================================
   //  Write results to file.
   // ===================================================================

   hfile.Write();

   return 0;
}
