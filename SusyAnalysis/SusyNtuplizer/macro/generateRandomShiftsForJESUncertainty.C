#include <fstream>
#include <iostream>
#include "TRandom3.h"

void generateRandomShiftsForJESUncertainty(const unsigned int nToys, const string& outFileName)
{
  //open output file to store 1 random shift per toy
  ofstream out(outFileName.c_str());
  if (!out.is_open()) {
    cerr << "Error opening file " << outFileName << ".  Quitting\n";
    return;
  }

  //write random shifts to file
  TRandom3 random(0);
  for (unsigned int iToy = 0; iToy < nToys; ++iToy) out << random.Gaus(0.0, 1.0) << endl;

  //close
  out.close();
}
