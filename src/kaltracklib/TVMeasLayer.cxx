//*************************************************************************
//* ====================
//*  TVMeasLayer Class
//* ====================
//*
//* (Description)
//*   Measurement layer interface class.
//* (Requires)
//* (Provides)
//*     class TVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2005/02/23  A.Yamaguchi       Added new data members, fFwdX0Inv,
//*                                 fBwdX0Inv and fIndex, and their
//*                                 corresponding getters and setters.
//*                                 Added a new method, GetX0Inv().
//*
//*************************************************************************
//

#include "TVMeasLayer.h"

ClassImp(TVMeasLayer)

TVMeasLayer::TVMeasLayer(TMaterial     &matIn,
                         TMaterial     &matOut)
           : fMaterialInPtr(&matIn),
             fMaterialOutPtr(&matOut),
             fIndex(0)
{
}

