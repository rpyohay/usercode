#ifndef BoostedTauAnalysis_Common_interface_Common_h
#define BoostedTauAnalysis_Common_interface_Common_h

#include <vector>
#include <algorithm>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"

class Common {

 public:

/*   static void sortByMinDR2(double, std::vector<double>&); */

  template<typename T, typename U>
    static std::vector<double> sortByProximity(const std::vector<T*>& objectsToSort, 
					       const U& referenceObject)
    {
      std::vector<double> sortedDR2;
      for (typename std::vector<T*>::const_iterator iObject = objectsToSort.begin(); 
	   iObject != objectsToSort.end(); ++iObject) {
	sortedDR2.push_back(reco::deltaR2(referenceObject, **iObject));
      }
      std::sort(sortedDR2.begin(), sortedDR2.end());
      return sortedDR2;
    }

  static void sortByPT(std::vector<reco::Candidate*>&);

 private:

  static bool comparePT(reco::Candidate*, reco::Candidate*);

};

#endif
