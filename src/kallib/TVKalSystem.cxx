//*************************************************************************
//* ====================
//*  TVKalSystem Class
//* ====================
//*
//* (Description)
//*   This is the base class of Kalman filter.
//* (Requires)
//* 	TObject
//* (Provides)
//* 	class TVKalSystem
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*
//*************************************************************************
//
#include <iostream>
#include "TVKalSystem.h"
#include "TVKalState.h"

//_____________________________________________________________________
//  ------------------------------
//  Base Class for measurement vector used by Kalman filter
//  ------------------------------
//
ClassImp(TVKalSystem)

TVKalSystem *TVKalSystem::fgCurInstancePtr = 0;

TVKalSystem::TVKalSystem(Int_t n) 
            :TObjArray(n),
             fCurSitePtr(0),
             fChi2(0.)
{
}

TVKalSystem::~TVKalSystem() 
{
}

//-------------------------------------------------------
// AddAndFilter 
//-------------------------------------------------------

Bool_t TVKalSystem::AddAndFilter(TVKalSite &next)
{
   //
   // Propagate current state to the next site
   //

   GetState(TVKalSite::kFiltered).Propagate(next);
   
   //
   // Calculate new pull and gain matrix
   //

   if (next.Filter()) {
      //
      // Add this to the system if accepted.
      //

      Add(&next);
      fChi2 += next.GetDeltaChi2();
      return kTRUE;
   } else {
      return kFALSE; 
   }
}

//-------------------------------------------------------
// GetNDF
//-------------------------------------------------------

Int_t TVKalSystem::GetNDF(Bool_t self)
{
   Int_t ndf    = 0;
   Int_t nsites = GetEntries();
   for (Int_t isite=1; isite<nsites; isite++) {
       TVKalSite &site = *(TVKalSite *)At(isite);
       if (!site.IsLocked()) ndf += site.GetDimension();
   }
   if (self) ndf -= GetCurSite().GetCurState().GetDimension();
   return ndf;
}

//-------------------------------------------------------
// SmoothBackTo 
//-------------------------------------------------------

void TVKalSystem::SmoothBackTo(Int_t k)
{

   TIter previous(this,kIterBackward);
   TIter cur     (this,kIterBackward);

   TVKalSite  *prePtr;
   TVKalSite  *curPtr = (TVKalSite *)cur();
   TVKalState &cura   = curPtr->GetState(TVKalSite::kFiltered);
   curPtr->Add(&curPtr->CreateState(cura, cura.GetCovMat(),
                                    TVKalSite::kSmoothed));

   while((curPtr = (TVKalSite *)cur()) && (prePtr = (TVKalSite *)previous())){
      curPtr->Smooth(*prePtr);
      fCurSitePtr = curPtr;
      if (IndexOf(curPtr) == k) break;
   }

}

//-------------------------------------------------------
// SmoothAll 
//-------------------------------------------------------

void TVKalSystem::SmoothAll()
{
   SmoothBackTo(0);
}
