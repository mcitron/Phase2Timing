#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"


JetTimingTools::JetTimingTools(edm::ConsumesCollector && cc):
    caloGeometryToken_(cc.esConsumes()),
    ecalCellEnergyThresh_(0.5),
    ecalCellTimeThresh_(12.5),
    ecalCellTimeErrorThresh_(100),
    matchingRadius2_(0.16)
{
}

void  JetTimingTools::setMatchingRadius(double matchingRadius) {matchingRadius2_ = matchingRadius*matchingRadius; }

void  JetTimingTools::setEcalCellEnergyThreshold(double ecalCellEnergyThresh){ ecalCellEnergyThresh_ = ecalCellEnergyThresh;}

void  JetTimingTools::setEcalCellTimeThreshold(double ecalCellTimeThresh){ ecalCellTimeThresh_ = ecalCellTimeThresh; }

void  JetTimingTools::setEcalCellTimeErrorThreshold(double ecalCellTimeErrorThresh){ ecalCellTimeErrorThresh_ = ecalCellTimeErrorThresh;}

void JetTimingTools::init(const edm::EventSetup &es){
  caloGeometry_ = es.getHandle(caloGeometryToken_);
}
//calculate jet time
void JetTimingTools::jetTimeFromEcalCells(
    const reco::Jet& jet,
    const edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>& ecalRecHits,
    float& weightedTimeCell,
    float& totalEmEnergyCell,
    uint& nCells) {
  for (auto const& ecalRH : ecalRecHits) {
    if (ecalRH.checkFlag(EcalRecHit::kSaturated) || ecalRH.checkFlag(EcalRecHit::kLeadingEdgeRecovered) ||
        ecalRH.checkFlag(EcalRecHit::kPoorReco) || ecalRH.checkFlag(EcalRecHit::kWeird) ||
        ecalRH.checkFlag(EcalRecHit::kDiWeird))
      continue;
    if (ecalRH.energy() < ecalCellEnergyThresh_)
      continue;
    if (ecalRH.timeError() <= 0. || ecalRH.timeError() > ecalCellTimeErrorThresh_)
      continue;
    if (fabs(ecalRH.time()) > ecalCellTimeThresh_)
      continue;
    auto const pos = caloGeometry_->getPosition(ecalRH.detid());
    if (reco::deltaR2(jet, pos) > matchingRadius2_)
      continue;
    weightedTimeCell += ecalRH.time() * ecalRH.energy() * sin(pos.theta());
    totalEmEnergyCell += ecalRH.energy() * sin(pos.theta());
    nCells++;
  }
  if (totalEmEnergyCell > 0) {
    weightedTimeCell /= totalEmEnergyCell;
  }
}
