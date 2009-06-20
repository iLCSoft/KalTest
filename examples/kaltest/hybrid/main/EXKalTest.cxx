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
#include <sstream>

#define SAVE_RESIDUAL

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

   Double_t pt      =  1.;   // default Pt [GeV]
   Double_t t0in    = 14.;   // default tp [nsec]
   Double_t cosmin  = -0.97; // default cos minimum [deg]
   Double_t cosmax  = 0.97;  // default cos maximum [deg]
   Int_t    nevents = 1;     // default number of events to generate
   switch (argc-offset) {
      case 6: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         cosmin  = atof(argv[4+offset]);
         cosmax  = atof(argv[5+offset]);
         break;
      case 5: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         cosmin = cosmax = atof(argv[4+offset]);
         break;
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
   //  Create TApplication
   // ===================================================================

   TApplication app("EXKalTest", &argc, argv, 0, 0);

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
   toygld.Close();          // close the cradle
   toygld.Sort();           // sort meas. layers from inside to outside

   //vtxdet.PowerOff();       // power off vtx not to process hit
   //itdet.PowerOff();        // power off it not to process hit
   //toygld.SwitchOffMS();    // switch off multiple scattering
   //toygld.SwitchOffDEDX();  // switch off enery loss

   // ===================================================================
   //  Prepare an output n-tuple
   // ===================================================================

   TFile hfile("h.root","RECREATE","KalTest");

   stringstream sout;
   sout << "ndf:chi2:cl:fi0:cpa:cs:t0"; // 7 items
#ifdef SAVE_RESIDUAL
   Int_t nitems  = 7;
   Int_t itemID[2000][2];
   TIter nextlayer(&toygld);
   TVMeasLayer *mlp;
   while ((mlp = dynamic_cast<TVMeasLayer *>(nextlayer()))) {
      TVMeasLayer &ml = *mlp;
      if (ml.IsActive()) {
         Int_t index = ml.GetIndex();
         sout << ":dxin" << setw(3) << setfill('0') << index;
	 itemID[index][0] = nitems++;
         sout << ":dxot" << setw(3) << setfill('0') << index;
	 itemID[index][1] = nitems++;
      }
   }
   Double_t *data = new Double_t [nitems];
#endif
   sout << ends;
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", sout.str().data());

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

      THelicalTrack hel = gen.GenerateHelix(pt, cosmin, cosmax);

      // ============================================================
      //  Swim the particle in detector
      // ============================================================

      gen.Swim(hel);

      // ============================================================
      //  Do Kalman Filter
      // ============================================================

      if (kalhits.GetEntries() < 3) {
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
#ifndef SAVE_RESIDUAL
      TVKalSite &cursite = kaltrack.GetCurSite();
#else
      Int_t isite = 1;
      kaltrack.SmoothBackTo(isite);
      TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);
#endif

      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t fi0  = cursite.GetCurState()(1, 0); 
      Double_t cpa  = cursite.GetCurState()(2, 0); 
      Double_t tnl  = cursite.GetCurState()(4, 0); 
      Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl);
      Double_t t0   = cursite.GetCurState()(5, 0); 
#ifndef SAVE_RESIDUAL
      hTrackMonitor->Fill(ndf, chi2, cl, fi0, cpa, cs, t0);
#else
      data[0] = ndf;
      data[1] = chi2;
      data[2] = cl;
      data[3] = fi0;
      data[4] = cpa;
      data[5] = cs;
      data[6] = t0;
      for (Int_t i=7; i<nitems; i++) data[i] = 9999999.;
      TIter nextsite(&kaltrack);
            nextsite(); // skip dummy site
      TKalTrackSite *sitep;
      while ((sitep = static_cast<TKalTrackSite *>(nextsite()))) {
         TKalTrackSite &site = *sitep;
         Int_t index = site.GetHit().GetMeasLayer().GetIndex();
	 data[itemID[index][0]] = site.GetResVec(TVKalSite::kSmoothed)(0,0);
	 if (site.InvFilter()) {
	    data[itemID[index][1]] = site.GetResVec(TVKalSite::kInvFiltered)(0,0);
	 }
      }
      hTrackMonitor->Fill(data);
#endif

      // ============================================================
      //  Very Primitive Event Display
      // ============================================================

      static TCanvas     *cvp    = 0;
      if (!gROOT->IsBatch()) {
         if (!cvp) {
            cvp = new TCanvas("OED", "Event Display", 400, 10, 610, 610);
         } else {
            cvp->cd();
            cvp->Clear();
         }

#if ROOT_VERSION_CODE < ROOT_VERSION(5,16,0)
         TView   *vwp = new TView(1);
#else
         TView   *vwp = TView::CreateView(1,0,0);
#endif
         vwp->SetRange(-260.,-260.,-260.,+260.,+260.,+260.);
         Int_t ierr;
         vwp->SetView(10.,80.,80.,ierr);

         vtxdet.Draw(40);
         itdet.Draw(40);
         tpcdet.Draw(40);
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
