//*************************************************************************
//* ===================
//*  EXKalSystem Class
//* ===================
//*
//* (Description)
//*   Sample steering class for Kalman filter
//* (Requires)
//*     TVKalSystem
//* (Provides)
//*     class EXKalSystem
//* (Update Recored)
//*   2003/09/30  Y.Nakashima   Original version.
//*
//*************************************************************************
                                                                                
#include "EXKalState.h"
#include "EXKalSite.h"
#include "EXKalSystem.h"

//_____________________________________________________________________
//  ------------------------------
//  Base Class for measurement vector used by Kalman filter
//  ------------------------------
//
ClassImp(EXKalSystem)
                                                                                
EXKalSystem::EXKalSystem(Int_t n)
            :TVKalSystem(n)
{
}
