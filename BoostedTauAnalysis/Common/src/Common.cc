#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

void minDR::setCandidate(const reco::Candidate* cand)
{
  cand_ = const_cast<reco::Candidate*>(cand);
}

void minDR::deleteCandidate() { cand_ = NULL; }

reco::Candidate* minDR::getCandidate() const { return cand_; }

void Common::sortByPT(std::vector<reco::Candidate*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePT);
}

bool Common::comparePT(reco::Candidate* object1, reco::Candidate* object2)
{
  return (object1->pt() < object2->pt());
}

reco::Vertex* Common::getPrimaryVertex(edm::Handle<reco::VertexCollection>& pVertices)
{
  reco::VertexCollection::const_iterator iVtx = pVertices->begin();
  reco::Vertex* pPV = NULL;
  while ((iVtx != pVertices->end()) && (pPV == NULL)) {
    if (!iVtx->isFake() && 
	(iVtx->ndof() > 4) && 
	(fabs(iVtx->x()) <= 24.0/*cm*/) && 
	(fabs(iVtx->position().Rho()) <= 2.0/*cm*/)) pPV = const_cast<reco::Vertex*>(&*iVtx);
    ++iVtx;
  }
  return pPV;
}

std::vector<reco::Muon*> Common::getTightRecoMuons(edm::Handle<reco::MuonCollection>& pMuons, 
						   const reco::Vertex* pPV)
{
  std::vector<reco::Muon*> tightMuons;
  for (reco::MuonCollection::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    if (muon::isTightMuon(*iMuon, *pPV) && 
	iMuon->isPFMuon() && 
	(fabs(iMuon->innerTrack()->dz(pPV->position())) < 0.5) && 
	(iMuon->track()->hitPattern().trackerLayersWithMeasurement() > 5)) {
      tightMuons.push_back(const_cast<reco::Muon*>(&*iMuon));
    }
  }
  return tightMuons;
}
