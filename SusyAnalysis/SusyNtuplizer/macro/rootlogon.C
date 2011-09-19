{
  cout << "Loading libSusy.so...\n";
  gSystem->Load("libSusy.so");
  cout << "Loading libJetMETObjects.so...\n";
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  cout << "Loading libFilters.so...\n";
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");

}
