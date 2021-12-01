#ifndef JetTimingTools_hh
#define JetTimingTools_hh
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/JetReco/interface/Jet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"

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
                            const edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>&,
                            float&,
                            float&,
                            uint&);
  void setMatchingRadius(double );
  void setEcalCellEnergyThreshold(double);
  void setEcalCellTimeThreshold(double );
  void setEcalCellTimeErrorThreshold(double );
private:
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeometry_;
    // Configurables for timing definition
    double ecalCellEnergyThresh_;
    double ecalCellTimeThresh_;
    double ecalCellTimeErrorThresh_;
    double matchingRadius2_;
};

#endif //JetTimingTools_hh
