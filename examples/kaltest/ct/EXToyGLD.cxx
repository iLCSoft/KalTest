//*************************************************************************
//* =====================
//*  EXToyGLD Class
//* =====================
//*
//* (Description)
//*   A singleton class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//*     TObjArray
//*     TVKalDetector
//* (Provides)
//*     class EXToyGLD
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi	Original version.
//*
//*************************************************************************

#include "EXToyGLD.h"

ClassImp(EXToyGLD)

EXToyGLD::EXToyGLD(Int_t n)
             : TKalDetCradle(n)
{
}

EXToyGLD::~EXToyGLD()
{
}

