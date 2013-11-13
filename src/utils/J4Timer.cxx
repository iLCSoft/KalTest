// $Id: J4Timer.cc,v 1.6 2004/02/03 15:01:59 fujiik Exp $
//*************************************************************************
//* --------------------
//* J4Timer
//* --------------------
//* (Description)
//* 	Class for describing his/her detector compornents.
//*     
//* (Update Record)
//*	2000/12/08  K.Hoshina	Original version.
//*************************************************************************

#include <sstream>
#include <iomanip>

#include "J4Timer.h"


// ====================================================================
//--------------------------------
// constants (detector parameters)
//--------------------------------

J4Timer::J4TimerArray J4Timer::fgTimers;
int                 J4Timer::fgNtimers = 0;

ClassImp(J4Timer)
//=====================================================================
//---------------------
// Class Description
//---------------------

// ====================================================================
//* constructor -------------------------------------------------------

J4Timer::J4Timer(int          &timerid, 
                 const string &classname, 
                 const string &timername)
        : fCurrentTimer(0)
{
          
   if (fgNtimers > (int)fgTimers.size()) {
      abort();
   } 

   if (timerid == -1) {
      AccumulatedTime *timer = new AccumulatedTime(fgNtimers,
                                                   classname,
                                                   timername);
      fgTimers.push_back(timer);
      timerid = fgNtimers;
      fgNtimers ++;


      std::cerr << "J4Timer::New timer is created! timerID, name = "
      << timerid << " " << classname << " " << timername << std::endl;
   } 
   
   fCurrentTimer = fgTimers[timerid];
   
}

// ====================================================================
//* destructor --------------------------------------------------------
J4Timer::~J4Timer()
{	
}

// ====================================================================
//* ResetAllTimers ----------------------------------------------------
void J4Timer::ResetAllTimers()
{
   for (int i=0; i<fgNtimers; i++) {
      if (fgTimers[i]) {
         fgTimers[i]->ResetTimes();
      }
   }
}

// ====================================================================
//* PrintAllAccumulatedTimes ------------------------------------------
void J4Timer::PrintAllAccumulatedTimes()
{
   std::cerr.precision(6);
   std::cerr << " *********************************************************************************" << std::endl;
   std::cerr << " * Output of Accumulated Time ****************************************************" << std::endl;
   std::cerr << " * ---------+---------+---------+---------+---------+---------+---------+---------" << std::endl;
   std::cerr << " * Timer Name                                   Real[s]   System[s]     User[s]   " 
	   << "#Calls"  << std::endl;
   std::cerr << " * ---------+---------+---------+---------+---------+---------+---------+---------" << std::endl;
   
   for (int i=0; i<fgNtimers; i++) {
      if (fgTimers[i]) {
         AccumulatedTime *timer = fgTimers[i];
         std::stringstream name;
         name << timer->GetClassName() << ":" << timer->GetTimerName();
         std::cerr << " * " << std::setw(40) << name.str()
                << std::setw(12) <<  timer->GetAccumulatedRealElapsed()
                << std::setw(12) <<  timer->GetAccumulatedSystemElapsed()
                << std::setw(12) <<  timer->GetAccumulatedUserElapsed() 
                << std::setw(12) <<  timer->GetNCalls() << std::endl;
      }
   }
   std::cerr << " *********************************************************************************" << std::endl;
}

