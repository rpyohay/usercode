#ifndef Skims_SkimTypes_CutScanKey_h
#define Skims_SkimTypes_CutScanKey_h

#include <string>
#include <vector>

class CutScanKey {

 public:

  //default constructor
  CutScanKey();

  //constructor from the cut name, bit number, and cut values
  explicit CutScanKey(const std::string&, const std::vector<double>&, const std::vector<double>&);

  //LCG dictionary doesn't like copy constructors
  //copy constructor
  //explicit CutScanKey(const CutScanKey&);

  //destructor
  ~CutScanKey();

  //assignment operator
  CutScanKey& operator=(const CutScanKey&);

  //setters
  void setCutName(const std::string&);
  void setCutVal1(const std::vector<double>&);
  void setCutVal2(const std::vector<double>&);

  //getters
  std::string getCutName() const;
  std::vector<double> getCutVal1() const;
  std::vector<double> getCutVal2() const;

 private:

  std::string cutName_;
  std::vector<double> cutVal1_;
  std::vector<double> cutVal2_;

};

typedef std::vector<CutScanKey> CutScanKeyCollection;

#endif
