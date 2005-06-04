//*************************************************************************
//* ======================
//*  TKalTrackState Class
//* ======================
//*
//* (Description)
//*   Track state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class TKalTrackState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  Y.Yamaguchi       Improved CalcProcessNoise().
//*
//*************************************************************************
//
#include <iostream>
#include <memory>

#include "TKalDetCradle.h"
#include "TVKalDetector.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TKalTrack.h"
#include "TObjNum.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//
// --------------------------------------------------------------------
// Ctors and Dtor
//
TKalTrackState::TKalTrackState(Int_t p) 
           : TVKalState(p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, Int_t type, Int_t p) 
           : TVKalState(sv,type,p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                                     Int_t type, Int_t p) 
           : TVKalState(sv,c,type,p), fX0()
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TVKalSite &site,
                                     Int_t type, Int_t p) 
           : TVKalState(sv,site,type,p), 
             fX0(((TKalTrackSite *)&site)->GetPivot())
{
}

TKalTrackState::TKalTrackState(const TKalMatrix &sv, const TKalMatrix &c,
                               const TVKalSite &site, Int_t type, Int_t p) 
           : TVKalState(sv,c,site,type,p),
             fX0(((TKalTrackSite *)&site)->GetPivot())
{
}

//
// --------------------------------------------------------------------
// Implementation of base-class pure virtuals
//

TKalTrackState * TKalTrackState::MoveTo(const TVKalSite  &to,
                                              TKalMatrix &F,
                                              TKalMatrix *QPtr) const
{
   TKalTrackSite &siteto = *(TKalTrackSite *)&to;
   TVector3       x0to   = siteto.GetPivot();

   auto_ptr<TVTrack> helto(&CreateTrack());

   TKalMatrix  av(5,1), Fto(5,5);
   Double_t    fid;
   helto->MoveTo(x0to,fid,&Fto);
   helto->PutInto(av);

   Int_t sdim = GetDimension();
   TKalMatrix sv(sdim,1);
   for (Int_t i=0; i<5; i++) {
      sv(i,0) = av(i,0);
      for (Int_t j=0; j<5; j++) {
         F(i,j) = Fto(i,j);
      }
   }
   if (sdim == 6) {
      sv(5,0) = (*this)(5,0);
      F (5,5) = 1.;
   }

   if (QPtr) {
      TKalTrackState *atop
                   = new TKalTrackState(sv, siteto, TVKalSite::kPredicted);
       //  fid = | phi_to - phi_from |
       //  assuming that pivots are hits
      *QPtr = CalcProcessNoise(siteto, *atop, *helto, fid);
      return atop;
   } else {
      return 0;
   }
}

TKalTrackState & TKalTrackState::MoveTo(const TVKalSite  &to,
                                              TKalMatrix &F,
                                              TKalMatrix &Q) const
{
   return *MoveTo(to, F, &Q);
}

TKalMatrix TKalTrackState::CalcProcessNoise(const TKalTrackSite  &to,
                                                  TKalTrackState &ato,
                                            const TVTrack        &tto,
                                                  Double_t        dfi) const
{
   const TKalTrackSite &from = static_cast<const TKalTrackSite &>(GetSite());
   TKalDetCradle       &det  = const_cast<TKalDetCradle &>
                              (static_cast<const TKalDetCradle &>
                              (from.GetHit().GetMeasLayer().GetParent()));
   Int_t p = GetDimension();
   TKalMatrix Q(p,p);
   if (!det.IsMSOn()) return Q;

   det.CalcTable(from, to);
   const TObjArray &mlt = det.GetMeasLayerTable();
   const TObjArray &dft = det.GetDPhiTable();
   Int_t dir = det.GetDir();

   Int_t  nel   = mlt.GetEntries();
   if (!nel) return Q;

   Double_t delfi = dfi;
   Double_t dr    = 0.;
   Double_t drp   = ato(0, 0);

   for (Int_t i = 0; i < nel; i++) {
      TKalMatrix D(p, p);
      tto.CalcDapDa(delfi, dr, drp, D);
      if (p == 6) D(5, 5) = 1.;
      TVMeasLayer         &ml  = *dynamic_cast<TVMeasLayer *>(mlt[i]);
      const TVKalDetector &dtr = dynamic_cast<const TVKalDetector &>
                                (ml.GetParent(kFALSE));
      Double_t   df  = dynamic_cast<TObjNum *>(dft[i])->GetNum();
      TKalMatrix Qi  = dtr.CalcSigmaMS0(ml.GetMaterial(dir), df, *this);
      //TKalMatrix Qi  = TKalMatrix(p, p);
      TKalMatrix Dt  = TKalMatrix(TMatrixD::kTransposed, D);
      Q += D * Qi * Dt;
      delfi -= df;
      dtr.CalcEnergyLoss(ml.GetMaterial(dir), df, *this, ato);
   }
   return Q;
}

void TKalTrackState::DebugPrint() const
{
   cerr << "          +-     -+   " << "+-" <<  endl
        << "          | drho  |   " << "| " << (*this)(0,0) << endl
        << "          | phi0  |   " << "| " << (*this)(1,0) << endl
        << " a      = | kappa | = " << "| " << (*this)(2,0) << endl
        << "          | dz    |   " << "| " << (*this)(3,0) << endl
        << "          | tanl  |   " << "| " << (*this)(4,0) << endl;
   if (GetDimension() == 6) {
      cerr 
        << "          | t0    |   " << "| " << (*this)(5,0) << endl;
   }
   cerr << "          +-     -+   " << "+-" << endl;
   cerr << "          +-" << endl 
        << " X0     = | " << fX0.X() << endl
        << "          | " << fX0.Y() << endl
        << "          | " << fX0.Z() << endl
        << "          +-" << endl;
   GetCovMat().DebugPrint(" covMat = ", 6);
}

//
// --------------------------------------------------------------------
// Derived class methods
//

THelicalTrack TKalTrackState::GetHelix() const
{
   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   return THelicalTrack(a,fX0,((TKalTrackSite *)&GetSite())->GetBfield());
}

TStraightTrack TKalTrackState::GetLine() const
{
   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   return TStraightTrack(a,fX0,((TKalTrackSite *)&GetSite())->GetBfield());
}

TVTrack &TKalTrackState::CreateTrack() const
{
   TVTrack *tkp = 0;

   TKalMatrix a(5,1);
   for (Int_t i=0; i<5; i++) a(i,0) = (*this)(i,0);
   Double_t bfield = static_cast<const TKalTrackSite &>(GetSite()).GetBfield();

   if (bfield == 0.) tkp = new TStraightTrack(a,fX0);
   else              tkp = new THelicalTrack(a,fX0, bfield);

   return *tkp;
}

