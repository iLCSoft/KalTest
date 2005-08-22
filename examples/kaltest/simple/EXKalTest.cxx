#define __DEBUG__

#include "EXKalTest.h"
#include "EXHit.h"
#include "EXKalState.h"
#include "EXKalSite.h"
#include "EXKalSystem.h"
#include "TFile.h"
#include "TH1D.h"


#define __TSTEP__     1.
#define __SLOPE__     0.2
#define __INTERCEPT__ 2.
#define __SIGMAX__    0.01

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   gROOT->SetBatch();
   TApplication app("EXKalTest", &argc, argv, 0, 0);

   //
   // -------------------------------------------------------------------
   // Define Hists and Plots
   // -------------------------------------------------------------------
   //
                                                                                
   TFile hfile("h.root","RECREATE","KalTest");
   TH1D *hNDFPtr   = new TH1D("hNDF","NDF",200,0.,200.);
   TH1D *hChi2Ptr  = new TH1D("hChi2","Chi2",200,0.,200.);
   TH1D *hCLPtr    = new TH1D("hCL","CL",200,0.,1.);

   //
   // -------------------------------------------------------------------
   // Start event loop
   // -------------------------------------------------------------------
   //
                                                                                
   Int_t ntrks = 1;
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

#if 0  // no dummy site
   // ---------------------------
   //  Create 1st site: site1
   // ---------------------------

   EXHit      &hit1   = *(EXHit *)next();
   EXKalSite  &site1  = *(new EXKalSite(hit1));

   // ---------------------------
   //  Set initial state to site1
   // ---------------------------

   TKalMatrix  sv1(2,1);
   sv1(0,0) = 0.;
   sv1(1,0) = hit1(0,0);

   TKalMatrix C(2,2);
#if 0
   C(0,0) = 1.e8;
   C(1,1) = 1.e8;
#else
   Double_t big = 1.e8;
   C(0,0) = 1 + big;
   C(0,1) = - big * hit1.GetT();
   C(1,0) = - big * hit1.GetT();
   C(1,1) = 1 + big * hit1.GetT() * hit1.GetT();
   C *= hit1(0,1) * hit1(0,1) / (hit1.GetT() * hit1.GetT() + 1);
#endif

   site1.Add(new EXKalState(sv1,C,TVKalSite::kPredicted));
   site1.Add(new EXKalState(sv1,C,TVKalSite::kFiltered));

   // ---------------------------
   //  Add site1 to the system
   // ---------------------------

   kalsys.Add(&site1);

#else // with dummy site
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

   TKalMatrix  svd(2,1);
   svd(0,0) = 0.;
   svd(1,0) = 0.;

   TKalMatrix C(2,2);
   C(0,0) = 1.e8;
   C(1,1) = 1.e8;

   sited.Add(new EXKalState(svd,C,TVKalSite::kPredicted));
   sited.Add(new EXKalState(svd,C,TVKalSite::kFiltered));

   // ---------------------------
   //  Add sited to the system
   // ---------------------------

   kalsys.Add(&sited);

#endif

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

   hNDFPtr ->Fill(ndf,1.);
   hChi2Ptr->Fill(chi2,1.);
   hCLPtr  ->Fill(cl,1.);

   }
   hfile.Write();
   return 0;
}
