#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"


JetTimingTools::JetTimingTools(edm::ConsumesCollector && cc):
    caloGeometryToken_(cc.esConsumes()),
    ecalCellEnergyThresh_(0.5),
    ecalCellTimeThresh_(12.5),
    ecalCellTimeErrorThresh_(100),
    hgcalTracksterEnergyThresh_(0.5),
    hgcalTracksterTimeThresh_(12.5),
    hgcalTracksterTimeErrorThresh_(100),
    matchingRadius2_(0.16)
{

}

void  JetTimingTools::setMatchingRadius(double matchingRadius) {matchingRadius2_ = matchingRadius*matchingRadius; }

void  JetTimingTools::setEcalCellEnergyThreshold(double ecalCellEnergyThresh){ ecalCellEnergyThresh_ = ecalCellEnergyThresh;}

void  JetTimingTools::setEcalCellTimeThreshold(double ecalCellTimeThresh){ ecalCellTimeThresh_ = ecalCellTimeThresh; }

void  JetTimingTools::setEcalCellTimeErrorThreshold(double ecalCellTimeErrorThresh){ ecalCellTimeErrorThresh_ = ecalCellTimeErrorThresh;}

void  JetTimingTools::setHgcalTracksterEnergyThreshold(double hgcalTracksterEnergyThresh){ hgcalTracksterEnergyThresh_ = hgcalTracksterEnergyThresh;}

void  JetTimingTools::setHgcalTracksterTimeThreshold(double hgcalTracksterTimeThresh){ hgcalTracksterTimeThresh_ = hgcalTracksterTimeThresh; }

void  JetTimingTools::setHgcalTracksterTimeErrorThreshold(double hgcalTracksterTimeErrorThresh){ hgcalTracksterTimeErrorThresh_ = hgcalTracksterTimeErrorThresh;}

void JetTimingTools::init(const edm::EventSetup &es){
  caloGeometry_ = es.getHandle(caloGeometryToken_);
}

//Calculate gen delay
std::vector<double> JetTimingTools::surfaceIntersection(const reco::GenParticle& genParticle, reco::GenParticle& genParticleMother, double radius){
    //all in cm
    // double radius = 130;
    double lightSpeed = 29979245800;

    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> particlePos = genParticle.vertex();
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> parentPos = genParticleMother.vertex();
    double genParticleBeta = genParticleMother.p()/genParticleMother.energy();
    double parentTime = TMath::Sqrt((particlePos-parentPos).Mag2())/(lightSpeed*genParticleBeta);
    double a = (genParticle.px()*genParticle.px() + genParticle.py()*genParticle.py())
    *lightSpeed*lightSpeed/(genParticle.energy()*genParticle.energy());
    double b = 2*(genParticle.vx()*genParticle.px() + genParticle.vy()*genParticle.py())*lightSpeed/genParticle.energy();
    double c = (genParticle.vx()*genParticle.vx() + genParticle.vy()*genParticle.vy()) - radius*radius;
    if (c > 0) return {1000.,1000.,1000.,1000.};
    double sqrt_disc = TMath::Sqrt(b*b - 4*a*c);
    double t1 = (-b + sqrt_disc)/(2*a);
    double xPos = t1*(genParticle.px()/genParticle.energy())*lightSpeed + genParticle.vx();
    double yPos = t1*(genParticle.py()/genParticle.energy())*lightSpeed + genParticle.vy();
    double zPos = t1*(genParticle.pz()/genParticle.energy())*lightSpeed + genParticle.vz();
    double phi = TMath::ATan2(yPos,xPos);
    double theta = TMath::ACos(zPos/TMath::Sqrt(xPos*xPos+yPos*yPos+zPos*zPos));
    double eta = -TMath::Log(TMath::Tan(theta/2.));
    // double t1Light = TMath::Sqrt(zPos*zPos+radius*radius)/lightSpeed;
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> intersectionXYZ(xPos,yPos,zPos);
    double t1Light = TMath::Sqrt((intersectionXYZ-parentPos).Mag2())/lightSpeed;
    std::vector<double> outVec = {eta,phi,t1Light,t1+parentTime-t1Light,t1};
    return outVec;
}
std::vector<double> JetTimingTools::endCapIntersection(const reco::GenParticle& genParticle, reco::GenParticle& genParticleMother, double zMin, double zMax){
    //all in cm
    // double radius = 130;
    double lightSpeed = 29979245800;

    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> particlePos = genParticle.vertex();
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> parentPos = genParticleMother.vertex();
    double genParticleBeta = genParticleMother.p()/genParticleMother.energy();
    double parentTime = TMath::Sqrt((particlePos-parentPos).Mag2())/(lightSpeed*genParticleBeta);
    double xPos;
    double yPos;
    double zPos;
    double t1;
    if (fabs(parentPos.z()) > zMax) return {1000.,1000.,1000.,1000.};
    else if (fabs(particlePos.z()) > zMin and fabs(particlePos.z()) < zMax)
    {
        xPos = particlePos.x();
        yPos = particlePos.y();
        zPos = particlePos.z();
	t1 = 0;
    }
    else {
	if (genParticle.pz() > 0) t1 = (zMin-genParticle.vz())/((genParticle.pz()/genParticle.energy())*lightSpeed);
	else t1 = (-zMin-genParticle.vz())/((genParticle.pz()/genParticle.energy())*lightSpeed);
	xPos = t1*(genParticle.px()/genParticle.energy())*lightSpeed + genParticle.vx();
	yPos = t1*(genParticle.py()/genParticle.energy())*lightSpeed + genParticle.vy();
	zPos = t1*(genParticle.pz()/genParticle.energy())*lightSpeed + genParticle.vz();
    }
    double phi = TMath::ATan2(yPos,xPos);
    double theta = TMath::ACos(zPos/TMath::Sqrt(xPos*xPos+yPos*yPos+zPos*zPos));
    double eta = -TMath::Log(TMath::Tan(theta/2.));
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> intersectionXYZ(xPos,yPos,zPos);
    double t1Light = TMath::Sqrt((intersectionXYZ-parentPos).Mag2())/lightSpeed;
    std::vector<double> outVec = {eta,phi,t1Light,t1+parentTime-t1Light,t1};
    return outVec;
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
void JetTimingTools::jetTimeFromHgcalTracksters(
    const reco::Jet& jet,
    const std::vector<ticl::Trackster>& tracksters,
    float& weightedTimeTrackster,
    float& totalEnergyTrackster,
    uint& nTracksters) {
    for (auto const& trackster : tracksters) {
    if (trackster.regressed_energy() < hgcalTracksterEnergyThresh_)
      continue;
    if (trackster.timeError() <= 0. || trackster.timeError() > hgcalTracksterTimeErrorThresh_)
      continue;
    if (fabs(trackster.time()) > hgcalTracksterTimeThresh_)
      continue;
    auto const pos = trackster.barycenter();
    if (reco::deltaR2(jet, pos) > matchingRadius2_)
      continue;
    weightedTimeTrackster += trackster.time() * trackster.regressed_energy() * sin(pos.theta());
    totalEnergyTrackster += trackster.regressed_energy() * sin(pos.theta());
    nTracksters++;
  }
  if (totalEnergyTrackster > 0) {
    weightedTimeTrackster /= totalEnergyTrackster;
    }
}
