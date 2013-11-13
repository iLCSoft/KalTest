#ifndef J4TIMER_H
#define J4TIMER_H
//*************************************************************************
//* -----------------------
//* J4Timer
//* -----------------------
//* (Description)
//*   Timer class for JUPITER.
//* (Usage)
//*   
//*     
//* (Update Record)
//*	2001/03/27  K.Hoshina	Original version.
//*************************************************************************

#include "TStopwatch.h"
#include <string>
#include <iostream>

using namespace std;

// typedefs ============================================================

// =====================================================================

class J4Timer : public TStopwatch
{
 
 public:

   J4Timer(int          &timerid, 
           const string &classname, 
           const string &timername);
          
   virtual ~J4Timer();
  
   static void    ResetAllTimers();
   static void    PrintAllAccumulatedTimes();
  
   inline virtual void     Start();
   inline virtual void     Stop();
   //inline virtual bool   IsValid() const;
   inline virtual double GetRealElapsed() ;
   inline virtual double GetSystemElapsed() ;
   inline virtual double GetUserElapsed() ;
   
   inline virtual void     ResetTimer(int id = -1);
   inline virtual string GetClassName(int id = -1) const;
   inline virtual string GetTimerName(int id = -1) const;
   inline virtual double GetAccumulatedRealElapsed(int id = -1) const;
   inline virtual double GetAccumulatedSystemElapsed(int id = -1) const;
   inline virtual double GetAccumulatedUserElapsed(int id = -1) const;
   inline virtual int    GetNCalls(int id = -1) const;

 private:
   // class AccumulatedTime
   class AccumulatedTime 
   {
    public:
       AccumulatedTime(int          &id,
                       const string &classname,
                       const string &timername)
                       :fID(id), fClassName(classname), 
                        fTimerName(timername),
                        fAccumulatedRealElapsed(0),
                        fAccumulatedSystemElapsed(0),
                        fAccumulatedUserElapsed(0),
	                    fNCalls(0) {}
       virtual ~AccumulatedTime(){}
       
       inline void ResetTimes() 
                   { 
                      fAccumulatedRealElapsed = 0;
                      fAccumulatedSystemElapsed = 0;
                      fAccumulatedUserElapsed = 0;
					  fNCalls = 0;
                   }
       inline void AccumulateRealElapsed(double t) 
                             { fAccumulatedRealElapsed += t; }
       inline void AccumulateSystemElapsed(double t) 
                             { fAccumulatedSystemElapsed += t; }
       inline void AccumulateUserElapsed(double t) 
                             { fAccumulatedUserElapsed += t; }
       inline void AccumulateNCalls() 
                             { fNCalls++; }
       
       inline int    GetID() const { return fID; }
       inline string GetClassName() const { return fClassName; }
       inline string GetTimerName() const { return fTimerName; }
       inline double GetAccumulatedRealElapsed() const 
                             { return fAccumulatedRealElapsed; }
       inline double GetAccumulatedSystemElapsed() const 
                             { return fAccumulatedSystemElapsed; }
       inline double GetAccumulatedUserElapsed() const 
                             { return fAccumulatedUserElapsed; }
	   inline int    GetNCalls() const
		                     { return fNCalls; }
       
    private:
       int      fID;
       string   fClassName;
       string   fTimerName;
       double   fAccumulatedRealElapsed;
       double   fAccumulatedSystemElapsed;
       double   fAccumulatedUserElapsed;
	   int      fNCalls;
   };
   // class AccumulatedTime end
                     
typedef  vector<J4Timer::AccumulatedTime*> J4TimerArray;
 
private:
  static int           fgNtimers;
  static J4TimerArray  fgTimers;
  AccumulatedTime      *fCurrentTimer;

  ClassDef(J4Timer,1)
    
};

//=====================================================================
//---------------------
// inline functions
//---------------------

void J4Timer::Start()
{
   this->TStopwatch::Start();
}

void J4Timer::Stop()
{
   this->TStopwatch::Stop();
   fCurrentTimer->AccumulateUserElapsed(GetUserElapsed());
   fCurrentTimer->AccumulateSystemElapsed(GetSystemElapsed());
   fCurrentTimer->AccumulateRealElapsed(GetRealElapsed());
   fCurrentTimer->AccumulateNCalls();
}

//G4bool J4Timer::IsValid() const
//{
//   return this->G4Timer::IsValid();
//}

double J4Timer::GetRealElapsed() 
{
   return this->TStopwatch::RealTime();
}

double J4Timer::GetSystemElapsed() 
{
   return this->TStopwatch::CpuTime();
}

double J4Timer::GetUserElapsed() 
{
   return this->TStopwatch::CpuTime();
}

void J4Timer::ResetTimer(int id) 
{
   if (id == -1) {
      fCurrentTimer->ResetTimes();
   } else if (id < fgNtimers) {
      fgTimers[id]->ResetTimes();
   } else {
      std::cerr << "J4Timer::ResetTimer: your id exceeds current fgNtimers."
             << " abort. your id = " << id << std::endl;
      abort();
   }
}

string J4Timer::GetClassName(int id) const 
{
   if (id == -1) {
      return fCurrentTimer->GetClassName();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetClassName();
   } else {
      std::cerr << "J4Timer::GetClassName: your id exceeds current fgNtimers."
             << " abort. your id = " << id << std::endl;
      abort();
   }
}

string J4Timer::GetTimerName(int id) const 
{
   if (id == -1) {
      return fCurrentTimer->GetTimerName();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetTimerName();
   } else {
      std::cerr << "J4Timer::GetTimerName: your id exceeds current fgNtimers."
             << " abort. your id = " << id << std::endl;
      abort();
   }
}

double J4Timer::GetAccumulatedRealElapsed(int id) const 
{
   if (id == -1) {
      return fCurrentTimer->GetAccumulatedRealElapsed();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetAccumulatedRealElapsed();
   } else {
      std::cerr << "J4Timer::GetAccumulatedRealElapsed: your id exceeds"
             << " current fgNtimers. abort. your id = " << id << std::endl;
      abort();
   }
}

double J4Timer::GetAccumulatedSystemElapsed(int id) const 
{
   if (id == -1) {
      return fCurrentTimer->GetAccumulatedSystemElapsed();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetAccumulatedSystemElapsed();
   } else {
      std::cerr << "J4Timer::GetAccumulatedSystemElapsed: your id exceeds"
             << " current fgNtimers. abort. your id = " << id << std::endl;
      abort();
   }
}

double J4Timer::GetAccumulatedUserElapsed(int id) const 
{
   if (id == -1) {
      return fCurrentTimer->GetAccumulatedUserElapsed();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetAccumulatedUserElapsed();
   } else {
      std::cerr << "J4Timer::GetAccumulatedUserElapsed: your id exceeds"
             << " current fgNtimers. abort. your id = " << id << std::endl;
      abort();
   }
}

int J4Timer::GetNCalls(int id) const
{
   if (id == -1) {
      return fCurrentTimer->GetNCalls();
   } else if (id < fgNtimers) {
      return fgTimers[id]->GetNCalls();
   } else {
      std::cerr << "J4Timer::GetNCalls: your id exceeds"
             << " current fgNtimers. abort. your id = " << id << std::endl;
      abort();
   }
}

#endif
