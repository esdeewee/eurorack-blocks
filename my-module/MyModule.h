// MyModule.h

#pragma once

#include "artifacts/MyModuleUi.h"
//#include "artifacts/MyModuleData.h" Uncomment if you use some 'data' resources
#include "erb/erb.h"


struct MyModule
{
   // The UI elements defined in MyModule.erbui are in 'ui'
   MyModuleUi ui;

   // The data resources defined in MyModule.erbb are in 'data'
   // Uncomment if you use some 'data' resources
   //MyModuleData data;

   void  init ();
   void  process ();

   // Put here your DSP objects
};
