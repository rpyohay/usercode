#include "BoostedTauAnalysis/Common/interface/Common.h"

// void Common::sortByMinDR2(double dr2, std::vector<double>& minDR2Vec)
// {
//   bool foundSpot = false;
//   unsigned int i = 0;
//   while ((i < minDR2Vec.size()) && (!foundSpot)) {
//     if (dr2 < minDR2Vec[i]) {
//       for (unsigned int j = minDR2Vec.size() - 1; j >= (i + 1); --j) {
// 	minDR2Vec[j] = minDR2Vec[j - 1];
//       }
//       minDR2Vec[i] = dr2;
//       foundSpot = true;
//     }
//     ++i;
//   }
// }

void Common::sortByPT(std::vector<reco::Candidate*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePT);
}

bool Common::comparePT(reco::Candidate* object1, reco::Candidate* object2)
{
  return (object1->pt() < object2->pt());
}
