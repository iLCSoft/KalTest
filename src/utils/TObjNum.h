#ifndef __TOBJNUM__
#define __TOBJNUM__

//*************************************************************************
//* -------
//* TObjNum
//* -------
//* (Description)
//*     A simple class to hold an integer.
//* (Update Record)
//*	2000/05/16  K.Fujii	Original version.
//*************************************************************************


#include "TObject.h"
#include <iostream>

class TObjNum : public TObject {
   private:
     Double_t  num; //! data
   public:
     TObjNum(Double_t i = 0) : num(i) { }
     ~TObjNum() { }
     void    SetNum(Double_t i) { num = i; }
     Double_t GetNum() { return num; }
     void    Print(Option_t *) { std::cerr << "TObjNum = " << num << std::endl; }
     Bool_t   IsEqual(TObject *obj) { return num == ((TObjNum*)obj)->num; }
     Bool_t   IsSortable() const { return kTRUE; }
     Int_t    Compare(const TObject *obj) const { 
	if (num > ((TObjNum*)obj)->num) 	return 1;
        else if (num < ((TObjNum*)obj)->num) 	return -1;
        else 					return 0; 
     }
   ClassDef(TObjNum, 1) // Double_t object
};

#endif
