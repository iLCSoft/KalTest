//#define __DEBUG__

#include "EXKalTest.h"
#include "EXHit.h"
#include "EXKalState.h"
#include "EXKalSite.h"
#include "EXKalSystem.h"
#include "TFile.h"
#include "TNtupleD.h"

#define __TSTEP__     1.
#define __SLOPE__     0.2
#define __INTERCEPT__ 2.
#define __SIGMAX__    0.01

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   // -------------------------------------------------------------------
   //  Define Hists and Plots
   // -------------------------------------------------------------------

   TFile hfile("h.root","RECREATE","KalTest");
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl");

   // -------------------------------------------------------------------
   //  Start event loop
   // -------------------------------------------------------------------

   Int_t ntrks = 1000;
   for (Int_t itrk=0; itrk<ntrks; itrk++) {
                                                                                
   cerr << " -------------------------------------------------- " << endl;
   cerr << " -------------- " << endl;
   cerr << " track: " << itrk << endl;
   cerr << " -------------- " << endl;
                                                                                
   // ---------------------------
   //  Create hits
   // ---------------------------

   Int_t nhits = 50;
   TObjArray hits(nhits);
   
   for (Int_t ihit=0; ihit<nhits; ihit++) {
       Double_t t  = __TSTEP__*ihit;
       Double_t dx = __SIGMAX__;
       Double_t x  = __SLOPE__*t + __INTERCEPT__;
       x += gRandom->Gaus(0.,dx);
       hits.Add(new EXHit(t,&x,&dx));
       //cerr << "t = " << t << " x = " << x << " dx = " << dx << endl;
   }

   // ---------------------------
   //  Create Kalman Filter
   // ---------------------------

   EXKalSystem kalsys;

   // ---------------------------
   // Prepare hit iterrator 
   // ---------------------------

   TIter next(&hits);

   // ---------------------------
   //  Create a dummy site: sited 
   // ---------------------------

   if (!hits.GetEntries()) continue;
                                                                                
   EXHit       hitd        = *(EXHit *)next();
               hitd(0,1)   = 1.e4;      // give a huge error to x
                                                                                
   next.Reset();                        // rewind iterator
                                                                                
   EXKalSite  &sited       = *(new EXKalSite(hitd));
                                                                                
   // sited.Lock();                     // dummy site should not be used
   sited.SetOwner();                    // site owns states

   // ---------------------------
   //  Set dummy state to sited
   // ---------------------------

   Int_t i1st  = 0;
   Int_t ilst  = hits.GetEntries() - 1;
   EXHit &h1st = *(EXHit *)hits.At(i1st);
   EXHit &hlst = *(EXHit *)hits.At(ilst);

   TKalMatrix  svd(2,1);
   svd(0,0) = (hlst.GetX(0) - h1st.GetX(0))/(hlst.GetT() - h1st.GetT());
   svd(1,0) = h1st.GetX(0);

   TKalMatrix C(2,2);
   C(0,0) = 1.e8;
   C(1,1) = 1.e8;

   sited.Add(new EXKalState(svd,C,TVKalSite::kPredicted));
   sited.Add(new EXKalState(svd,C,TVKalSite::kFiltered));

   // ---------------------------
   //  Add sited to the system
   // ---------------------------

   kalsys.Add(&sited);

   // ---------------------------
   //  Start Kalman Filter
   // ---------------------------

   Int_t  loop = 1;
   EXHit *hitPtr;
   Int_t    ndf  = 0;
   Double_t chi2 = 0.;
   Double_t cl   = 0.;
   while ((hitPtr = (EXHit *)next())) {
#ifdef __DEBUG__ 
      cerr << "----------------------" << endl
           << "loop = " << loop        << endl;
#endif
      EXHit &hit = *hitPtr;
      EXKalSite & site = *(new EXKalSite(hit));
      if (kalsys.AddAndFilter(site)) {
#ifdef __DEBUG__
         site.DebugPrint();
         kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
#endif
         ndf  = kalsys.GetNDF();
         chi2 = kalsys.GetChi2();
         cl   = TMath::Prob(chi2,ndf);
#ifdef __DEBUG__ 
         cerr << " ndf = "  << ndf
              << " chi2 = " << chi2 << endl;
         cerr << " cl = " << cl << endl;
#endif
      } else {
         cerr << " site discarded!" << endl;
         delete &site;
      }
      loop++;
   }

   kalsys.SmoothBackTo(2);
#ifdef __DEBUG__ 
   cerr << " Smoothed-------------------" << endl;
   kalsys.GetState(TVKalSite::kSmoothed).DebugPrint();
#endif   

   hTrackMonitor->Fill(ndf,chi2,cl);

   }
   hfile.Write();
   return 0;
}
