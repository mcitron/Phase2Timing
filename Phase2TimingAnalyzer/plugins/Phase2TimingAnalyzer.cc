// -*- C++ -*-
//
// Package:    Phase2TimingAnalyzer/Phase2TimingAnalyzer
// Class:      Phase2TimingAnalyzer
//
/**\class Phase2TimingAnalyzer Phase2TimingAnalyzer.cc Phase2TimingAnalyzer/Phase2TimingAnalyzer/plugins/Phase2TimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Tue, 30 Nov 2021 21:39:01 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"

//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "TMath.h"
#include "TTree.h"


//
// class declaration
//
struct tree_struc_{
  float                         q1_eta;
  float                         q2_eta;
  float                         q3_eta;
  float                         q4_eta;
  float                         q1_phi;
  float                         q2_phi;
  float                         q3_phi;
  float                         q4_phi;
  float                         q1_pt;
  float                         q2_pt;
  float                         q3_pt;
  float                         q4_pt;
  float                         q1_vx;
  float                         q2_vx;
  float                         q3_vx;
  float                         q4_vx;
  float                         q1_vy;
  float                         q2_vy;
  float                         q3_vy;
  float                         q4_vy;
  float                         q1_vz;
  float                         q2_vz;
  float                         q3_vz;
  float                         q4_vz;
  int                           nrecojets;
  std::vector<float>            recojet_pt;
  std::vector<float>            recojet_eta;
  std::vector<float>            recojet_phi;
  std::vector<float>            recojet_e;
  std::vector<float>            recojet_ECALtime;
  std::vector<float>            recojet_ECALenergy;
  std::vector<float>            recojet_ECALnCells;
  std::vector<float>            recojet_MTDtime;
  std::vector<float>            recojet_MTDenergy;
  std::vector<float>            recojet_MTDnCells;
  std::vector<float>            recojet_HGCALtime;
  std::vector<float>            recojet_HGCALenergy;
  std::vector<float>            recojet_HGCALnCells;
};

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class Phase2TimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Phase2TimingAnalyzer(const edm::ParameterSet&); 
  ~Phase2TimingAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void initTreeStructure();
  void clearVectors();
  JetTimingTools _jetTimingTools;

  // ---------- member data -------------------- // 
  edm::Service<TFileService> fs;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > _genParticles; 
  edm::Handle< edm::View<reco::GenParticle> > _genParticlesH;
  const edm::EDGetTokenT< edm::View<reco::PFJet> > _recoak4PFJets; 
  edm::Handle< edm::View<reco::PFJet> > _recoak4PFJetsH;
  const edm::EDGetTokenT<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEBToken_;
  edm::Handle< edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>> > _ecalRecHitsEBH;

  // setup tree;                                                                                                                                             
  TTree* tree;
  tree_struc_ tree_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Phase2TimingAnalyzer::Phase2TimingAnalyzer(const edm::ParameterSet& iConfig):
  _jetTimingTools(consumesCollector()),
  _genParticles(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  _genParticlesH(),
  _recoak4PFJets(consumes< edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("recoak4PFJets"))),
  _recoak4PFJetsH(),
  ecalRecHitsEBToken_{consumes<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>>(												       iConfig.getParameter<edm::InputTag>("ebRecHitsColl"))},
  _ecalRecHitsEBH()
{

}

Phase2TimingAnalyzer::~Phase2TimingAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void Phase2TimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif

  //  _jetTimingTools.init(iSetup);

  iEvent.getByToken(_genParticles, _genParticlesH);
  iEvent.getByToken(_recoak4PFJets, _recoak4PFJetsH);
  iEvent.getByToken(ecalRecHitsEBToken_, _ecalRecHitsEBH);


  //  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);

  //variable declaration
  float q1_pt  = -10000;
  float q2_pt  = -10000;
  float q3_pt  = -10000;
  float q4_pt  = -10000;
  float q1_phi  = -10000;
  float q2_phi  = -10000;
  float q3_phi  = -10000;
  float q4_phi  = -10000;
  float q1_eta  = -10000;
  float q2_eta  = -10000;
  float q3_eta  = -10000;
  float q4_eta  = -10000;
  float q1_vx   = -10000;
  float q1_vy   = -10000;
  float q1_vz   = -10000;
  float q2_vx   = -10000;
  float q2_vy   = -10000;
  float q2_vz   = -10000;
  float q3_vx   = -10000;
  float q3_vy   = -10000;
  float q3_vz   = -10000;
  float q4_vx   = -10000;
  float q4_vy   = -10000;
  float q4_vz   = -10000;
  int interestingquark = 0;
  int nrecojets = 0;
  std::vector<float>    recojet_pt;
  std::vector<float>    recojet_eta;
  std::vector<float>    recojet_phi;
  std::vector<float>    recojet_e;
  std::vector<float>    recojet_ECALtime;
  std::vector<float>    recojet_ECALenergy;
  std::vector<float>    recojet_ECALnCells;
  std::vector<float>    recojet_MTDtime;
  std::vector<float>    recojet_MTDenergy;
  std::vector<float>    recojet_MTDnCells;
  std::vector<float>    recojet_HGCALtime;
  std::vector<float>    recojet_HGCALenergy;
  std::vector<float>    recojet_HGCALnCells;

  
  bool debug=1;


  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON GENPARTICLES --------------------------------------"<<std::endl; 
  for (const auto & genpar_iter : *_genParticlesH){

    if (genpar_iter.mother(0) == NULL)continue;
    if(abs(genpar_iter.pdgId()) !=5  || genpar_iter.status() != 23)continue;
    float vx = genpar_iter.vertex().x();
    float vy = genpar_iter.vertex().y();
    float vz = genpar_iter.vertex().z();
    interestingquark++;
    if (interestingquark==1){ q1_pt = genpar_iter.pt(); q1_eta = genpar_iter.eta(); q1_phi = genpar_iter.phi(); q1_vx = vx; q1_vy = vy; q1_vz = vz;}
    if (interestingquark==2){ q2_pt = genpar_iter.pt(); q2_eta = genpar_iter.eta(); q2_phi = genpar_iter.phi(); q2_vx = vx; q2_vy = vy; q2_vz = vz;}
    if (interestingquark==3){ q3_pt = genpar_iter.pt(); q3_eta = genpar_iter.eta(); q3_phi = genpar_iter.phi(); q3_vx = vx; q3_vy = vy; q3_vz = vz;}
    if (interestingquark==4){ q4_pt = genpar_iter.pt(); q4_eta = genpar_iter.eta(); q4_phi = genpar_iter.phi(); q4_vx = vx; q4_vy = vy; q4_vz = vz;}

  }

  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);

  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON RECO JETS --------------------------------------"<<std::endl; 
  for (const auto & recojet_iter : *_recoak4PFJetsH){

    if(recojet_iter.pt()<10)continue;
    if(fabs(recojet_iter.eta())>3)continue;
    nrecojets++;
    recojet_pt.push_back(recojet_iter.pt());
    recojet_eta.push_back(recojet_iter.eta());
    recojet_phi.push_back(recojet_iter.phi());
    recojet_e.push_back(recojet_iter.energy());
    


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM ECAL --------------------------------------"<<std::endl; 
    float weightedECALTimeCell = 0;
    float totalECALEnergyCell = 0;
    unsigned int ECALnCells = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromEcalCells(recojet_iter, ecalRecHitsEB, weightedECALTimeCell, totalECALEnergyCell, ECALnCells);
    recojet_ECALenergy.push_back(totalECALEnergyCell);
    recojet_ECALnCells.push_back(ECALnCells);
    if(ECALnCells>0)
      recojet_ECALtime.push_back(weightedECALTimeCell);
    else
      recojet_ECALtime.push_back(-50);



    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD --------------------------------------"<<std::endl; 
    float weightedMTDTimeCell = 0;
    float totalMTDEnergyCell = 0;
    unsigned int MTDnCells = 0;
    /// insert here jet timing tools function
    recojet_MTDenergy.push_back(totalMTDEnergyCell);
    recojet_MTDnCells.push_back(MTDnCells);
    if(MTDnCells>0)
      recojet_MTDtime.push_back(weightedMTDTimeCell);
    else
      recojet_MTDtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HGCAL --------------------------------------"<<std::endl; 
    float weightedHGCALTimeCell = 0;
    float totalHGCALEnergyCell = 0;
    unsigned int HGCALnCells = 0;
    /// insert here jet timing tools function
    recojet_HGCALenergy.push_back(totalHGCALEnergyCell);
    recojet_HGCALnCells.push_back(HGCALnCells);
    if(HGCALnCells>0)
      recojet_HGCALtime.push_back(weightedHGCALTimeCell);
    else
      recojet_HGCALtime.push_back(-50);

    
  }
  
  

  // --- setup tree values                                                                                                                                   
  initTreeStructure();
  clearVectors();
  tree_.q1_eta = q1_eta;
  tree_.q2_eta = q2_eta;
  tree_.q3_eta = q3_eta;
  tree_.q4_eta = q4_eta;
  tree_.q1_phi = q1_phi;
  tree_.q2_phi = q2_phi;
  tree_.q3_phi = q3_phi;
  tree_.q4_phi = q4_phi;
  tree_.q1_pt = q1_pt;
  tree_.q2_pt = q2_pt;
  tree_.q3_pt = q3_pt;
  tree_.q4_pt = q4_pt;
  tree_.q1_vx = q1_vx;
  tree_.q2_vx = q2_vx;
  tree_.q3_vx = q3_vx;
  tree_.q4_vx = q4_vx;
  tree_.q1_vy = q1_vy;
  tree_.q2_vy = q2_vy;
  tree_.q3_vy = q3_vy;
  tree_.q4_vy = q4_vy;
  tree_.q1_vz = q1_vz;
  tree_.q2_vz = q2_vz;
  tree_.q3_vz = q3_vz;
  tree_.q4_vz = q4_vz;

  tree_.nrecojets        = nrecojets;
  for (int ij = 0; ij < nrecojets; ij++){
    tree_.recojet_pt.push_back(recojet_pt[ij]);
    tree_.recojet_eta.push_back(recojet_eta[ij]);
    tree_.recojet_phi.push_back(recojet_phi[ij]); 
    tree_.recojet_e.push_back(recojet_e[ij]);
    tree_.recojet_ECALtime.push_back(recojet_ECALtime[ij]);
    tree_.recojet_ECALenergy.push_back(recojet_ECALenergy[ij]);
    tree_.recojet_ECALnCells.push_back(recojet_ECALnCells[ij]);
    tree_.recojet_MTDtime.push_back(recojet_MTDtime[ij]);
    tree_.recojet_MTDenergy.push_back(recojet_MTDenergy[ij]);
    tree_.recojet_MTDnCells.push_back(recojet_MTDnCells[ij]);
    tree_.recojet_HGCALtime.push_back(recojet_HGCALtime[ij]);
    tree_.recojet_HGCALenergy.push_back(recojet_HGCALenergy[ij]);
    tree_.recojet_HGCALnCells.push_back(recojet_HGCALnCells[ij]);
  }


  // --- fill tree
  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void Phase2TimingAnalyzer::beginJob() {

  bool verbose_ = true;
  if (verbose_) std::cout << "Starting job" << std::endl;
  // --- set up output tree                                                                                                                                  
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("q1_eta", &tree_.q1_eta, "q1_eta/F");
  tree->Branch("q2_eta", &tree_.q2_eta, "q2_eta/F");
  tree->Branch("q3_eta", &tree_.q3_eta, "q3_eta/F");
  tree->Branch("q4_eta", &tree_.q4_eta, "q4_eta/F");
  tree->Branch("q1_phi", &tree_.q1_phi, "q1_phi/F");
  tree->Branch("q2_phi", &tree_.q2_phi, "q2_phi/F");
  tree->Branch("q3_phi", &tree_.q3_phi, "q3_phi/F");
  tree->Branch("q4_phi", &tree_.q4_phi, "q4_phi/F");
  tree->Branch("q1_pt", &tree_.q1_pt, "q1_pt/F");
  tree->Branch("q2_pt", &tree_.q2_pt, "q2_pt/F");
  tree->Branch("q3_pt", &tree_.q3_pt, "q3_pt/F");
  tree->Branch("q4_pt", &tree_.q4_pt, "q4_pt/F");
  tree->Branch("q1_vx", &tree_.q1_vx, "q1_vx/F");
  tree->Branch("q2_vx", &tree_.q2_vx, "q2_vx/F");
  tree->Branch("q3_vx", &tree_.q3_vx, "q3_vx/F");
  tree->Branch("q4_vx", &tree_.q4_vx, "q4_vx/F");
  tree->Branch("q1_vy", &tree_.q1_vy, "q1_vy/F");
  tree->Branch("q2_vy", &tree_.q2_vy, "q2_vy/F");
  tree->Branch("q3_vy", &tree_.q3_vy, "q3_vy/F");
  tree->Branch("q4_vy", &tree_.q4_vy, "q4_vy/F");
  tree->Branch("q1_vz", &tree_.q1_vz, "q1_vz/F");
  tree->Branch("q2_vz", &tree_.q2_vz, "q2_vz/F");
  tree->Branch("q3_vz", &tree_.q3_vz, "q3_vz/F");
  tree->Branch("q4_vz", &tree_.q4_vz, "q4_vz/F");

  tree->Branch("nrecojets",              &tree_.nrecojets,                "nrecojets/I");
  tree->Branch("recoJet_pt",             &tree_.recojet_pt);
  tree->Branch("recoJet_eta",             &tree_.recojet_eta);
  tree->Branch("recoJet_phi",             &tree_.recojet_phi);
  tree->Branch("recoJet_e",             &tree_.recojet_e);
  tree->Branch("recoJet_ECALtime",             &tree_.recojet_ECALtime);
  tree->Branch("recoJet_ECALenergy",             &tree_.recojet_ECALenergy);
  tree->Branch("recoJet_ECALnCells",             &tree_.recojet_ECALnCells);
  tree->Branch("recoJet_MTDtime",             &tree_.recojet_MTDtime);
  tree->Branch("recoJet_MTDenergy",             &tree_.recojet_MTDenergy);
  tree->Branch("recoJet_MTDnCells",             &tree_.recojet_MTDnCells);
  tree->Branch("recoJet_HGCALtime",             &tree_.recojet_HGCALtime);
  tree->Branch("recoJet_HGCALenergy",             &tree_.recojet_HGCALenergy);
  tree->Branch("recoJet_HGCALnCells",             &tree_.recojet_HGCALnCells);
}


// ------------ initialize trees ------------                                                                                                                 
void Phase2TimingAnalyzer::initTreeStructure()
{
  tree_.q1_pt = -999;
  tree_.q2_pt = -999;
  tree_.q3_pt = -999;
  tree_.q4_pt = -999;
  tree_.q1_eta = -999;
  tree_.q2_eta = -999;
  tree_.q3_eta = -999;
  tree_.q4_eta = -999;
  tree_.q1_phi = -999;
  tree_.q2_phi = -999;
  tree_.q3_phi = -999;
  tree_.q4_phi = -999;
  tree_.q1_vx = -999;
  tree_.q2_vx = -999;
  tree_.q3_vx = -999;
  tree_.q4_vx = -999;
  tree_.q1_vy = -999;
  tree_.q2_vy = -999;
  tree_.q3_vy = -999;
  tree_.q4_vy = -999;
  tree_.q1_vz = -999;
  tree_.q2_vz = -999;
  tree_.q3_vz = -999;
  tree_.q4_vz = -999;

}


void Phase2TimingAnalyzer::clearVectors()
{

  tree_.recojet_pt.clear();
  tree_.recojet_eta.clear();
  tree_.recojet_phi.clear();
  tree_.recojet_e.clear();
  tree_.recojet_ECALtime.clear();
  tree_.recojet_ECALenergy.clear();
  tree_.recojet_ECALnCells.clear();
  tree_.recojet_MTDtime.clear();
  tree_.recojet_MTDenergy.clear();
  tree_.recojet_MTDnCells.clear();
  tree_.recojet_HGCALtime.clear();
  tree_.recojet_HGCALenergy.clear();
  tree_.recojet_HGCALnCells.clear();
}

// ------------ method called once each job just after ending the event loop  ------------
void Phase2TimingAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2TimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2TimingAnalyzer);
