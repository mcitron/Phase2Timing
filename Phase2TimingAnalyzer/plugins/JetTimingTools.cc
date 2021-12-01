#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"


JetTimingTools::JetTimingTools(edm::ConsumesCollector && cc):
    caloGeometryToken_(cc.esConsumes()),
    mtdGeometryToken_(cc.esConsumes()),
    mtdTopologyToken_(cc.esConsumes()),
    ecalCellEnergyThresh_(0.5),
    ecalCellTimeThresh_(12.5),
    ecalCellTimeErrorThresh_(100),
    matchingRadius2_(0.16),
    mtdCellEnergyThresh_(0.5),
    mtdCellTimeThresh_(25.),
    mtdCellTimeErrorThresh_(100)
{

}

void  JetTimingTools::setMatchingRadius(double matchingRadius) {matchingRadius2_ = matchingRadius*matchingRadius; }

void  JetTimingTools::setEcalCellEnergyThreshold(double ecalCellEnergyThresh){ ecalCellEnergyThresh_ = ecalCellEnergyThresh;}

void  JetTimingTools::setEcalCellTimeThreshold(double ecalCellTimeThresh){ ecalCellTimeThresh_ = ecalCellTimeThresh; }

void  JetTimingTools::setEcalCellTimeErrorThreshold(double ecalCellTimeErrorThresh){ ecalCellTimeErrorThresh_ = ecalCellTimeErrorThresh;}

void  JetTimingTools::setMTDCellEnergyThreshold(double mtdCellEnergyThresh){ mtdCellEnergyThresh_ = mtdCellEnergyThresh;}

void  JetTimingTools::setMTDCellTimeThreshold(double mtdCellTimeThresh){ mtdCellTimeThresh_ = mtdCellTimeThresh; }

void  JetTimingTools::setMTDCellTimeErrorThreshold(double mtdCellTimeErrorThresh){ mtdCellTimeErrorThresh_ = mtdCellTimeErrorThresh;}

void JetTimingTools::init(const edm::EventSetup &es){
  caloGeometry_ = es.getHandle(caloGeometryToken_);
  mtdGeometry_ = es.getHandle(mtdGeometryToken_);
  mtdTopology_ = es.getHandle(mtdTopologyToken_);
}

//calculate jet time from ecal cells
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



//calculate jet time from ecal cells
void JetTimingTools::jetTimeFromMTDCells(
    const reco::Jet& jet,
    const edm::SortedCollection<FTLRecHit,edm::StrictWeakOrdering<FTLRecHit> >& mtdRecHits,
    float& weightedTimeCell,
    float& totalEmEnergyCell,
    uint& nCells) {
  
  for (auto const& mtdRH : mtdRecHits) {
    
    if (mtdRH.energy() < mtdCellEnergyThresh_)
      continue;
    if (mtdRH.timeError() <= 0. || mtdRH.timeError() > mtdCellTimeErrorThresh_)
      continue;
    if (fabs(mtdRH.time()) > mtdCellTimeThresh_)
      continue;

    BTLDetId detId = mtdRH.id();
    DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(mtdTopology_->getMTDTopologyMode()));
    const MTDGeomDet* thedet = mtdGeometry_->idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());      
    Local3DPoint local_point(0., 0., 0.);
    local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
    const auto& global_point = thedet->toGlobal(local_point);


    if (reco::deltaR2(jet, global_point) > matchingRadius2_)
    continue;
    //    std::cout<< "MTD RECHIT: eta: "<<global_point.eta()<<" phi: "<<global_point.phi()<<" pt: "<<mtdRH.energy()* sin(global_point.theta())<<" time: "<<mtdRH.time()<<std::endl;
    weightedTimeCell += mtdRH.time() * mtdRH.energy() * sin(global_point.theta());
    totalEmEnergyCell += mtdRH.energy() * sin(global_point.theta());
    nCells++;
  }
  if (totalEmEnergyCell > 0) {
    weightedTimeCell /= totalEmEnergyCell;
  }
}
