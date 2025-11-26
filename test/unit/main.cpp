/*****************************************************************************

      main.cpp
      Copyright (c) 2020 Raphael DINGE

*Tab=3***********************************************************************/

#include "test.h"
#include <iostream>
#include <string>

// Include existing legacy tests if we want to run them
// We can manually register them or just run them.
#include "TestAnimation.h"
#include "TestDebounce.h"
#include "TestSdramPtr.h"

int main (int argc, char* argv[])
{
   std::string filter = "";
   if (argc > 1) {
       filter = argv[1];
   }
   
   if (filter.empty()) {
       std::cout << "Running Legacy Tests..." << std::endl;
       { std::cout << "TestDebounce: "; TestDebounce t; t.run(); std::cout << "OK" << std::endl; }
       { std::cout << "TestAnimation: "; TestAnimation t; t.run(); std::cout << "OK" << std::endl; }
       { std::cout << "TestSdramPtr: "; TestSdramPtr t; t.run(); std::cout << "OK" << std::endl; }
   }

   std::cout << "Running Registry Tests..." << std::endl;
   TestRegistry::instance().runFiltered(filter);

   std::cout << "All Tests Finished." << std::endl;
   return 0;
}
