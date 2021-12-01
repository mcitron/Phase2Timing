#ifndef JetTimingTools_hh
#define JetTimingTools_hh
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/JetReco/interface/Jet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
}  // namespace edm

class JetTimingTools {
public:
  explicit JetTimingTools(edm::ConsumesCollector &&);
  // virtual ~JetTimingTools(){};
  void init(const edm::EventSetup &es);
  void jetTimeFromEcalCells(const reco::Jet&,
                            const edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>&,
                            float&,
                            float&,
                            uint&);
  void jetTimeFromMTDCells(const reco::Jet&,
                            const edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>&,
                            float&,
                            float&,
                            uint&);
  void setMatchingRadius(double );
  void setEcalCellEnergyThreshold(double);
  void setEcalCellTimeThreshold(double );
  void setEcalCellTimeErrorThreshold(double );
  void setMTDCellEnergyThreshold(double);
  void setMTDCellTimeThreshold(double );
  void setMTDCellTimeErrorThreshold(double );

private:
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeometry_;
    
    edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdGeometryToken_;
    edm::ESHandle<MTDGeometry> mtdGeometry_;

    edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdTopologyToken_;
    edm::ESHandle<MTDTopology> mtdTopology_;

    
    // Configurables for timing definition
    double ecalCellEnergyThresh_;
    double ecalCellTimeThresh_;
    double ecalCellTimeErrorThresh_;
    double matchingRadius2_;
    double mtdCellEnergyThresh_;
    double mtdCellTimeThresh_;
    double mtdCellTimeErrorThresh_;

};

#endif //JetTimingTools_hh
