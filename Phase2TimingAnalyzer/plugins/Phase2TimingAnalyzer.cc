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
#include "TLorentzVector.h"

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TMath.h"
#include "TTree.h"
#include "TGeoPolygon.h"


//
// class declaration
//
struct tree_struc_{
    // std::vector<float> tracksterEM_time;
    // std::vector<float> tracksterHAD_time;
    // std::vector<float> tracksterTrkEM_time;
    std::vector<float> tracksterMerge_time;
    std::vector<float> tracksterMerge_ticlIteration;
    // std::vector<float> tracksterEM_timeError;
    // std::vector<float> tracksterHAD_timeError;
    // std::vector<float> tracksterTrkEM_timeError;
    std::vector<float> tracksterMerge_timeError;
    // std::vector<float> tracksterEM_e;
    // std::vector<float> tracksterHAD_e;
    // std::vector<float> tracksterTrkEM_e;
    std::vector<float> tracksterMerge_e;
    // std::vector<float> tracksterEM_eta;
    // std::vector<float> tracksterHAD_eta;
    // std::vector<float> tracksterTrkEM_eta;
    std::vector<float> tracksterMerge_eta;
    // std::vector<float> tracksterEM_phi;
    // std::vector<float> tracksterHAD_phi;
    // std::vector<float> tracksterTrkEM_phi;
    std::vector<float> tracksterMerge_phi;
    std::vector<uint> tracksterMerge_iJ;
    std::vector<float> pvtrack_eta;
    std::vector<float> pvtrack_phi;
    std::vector<float> pvtrack_pt;
    std::vector<float> pfCand_eta;
    std::vector<float> pfCand_phi;
    std::vector<int> pfCand_iJ;
    std::vector<float> pfCand_pt;
    std::vector<float> pfCand_time;
    std::vector<float> pfCand_timeError;
    std::vector<float> q_eta;
    std::vector<float> q_phi;
    std::vector<float> q_pt;
    std::vector<int> q_pdgId;
    std::vector<float> q_ctau;
    std::vector<float> q_vx;
    std::vector<float> q_vy;
    std::vector<float> q_vz;
    std::vector<float> q_ebeta;
    std::vector<float> q_ebphi;
    std::vector<float> q_ebdelay;
    std::vector<float> q_etleta;
    std::vector<float> q_etlphi;
    std::vector<float> q_etldelay;
    std::vector<float> q_hgeta;
    std::vector<float> q_pathdelay;
    std::vector<bool> q_decayInHGCAL;
    std::vector<float> q_hgphi;
    std::vector<float> q_hgdelay;
  int                           ngen;
  int 				npvtrack;
  int 				npfCand;
  int                           nrecojets;
  std::vector<float>            recojet_pt;
  std::vector<float>            recojet_eta;
  std::vector<float>            recojet_phi;
  std::vector<float>            recojet_e;
  std::vector<float>            recojet_ECALtime;
  std::vector<float>            recojet_ECALenergy;
  std::vector<float>            recojet_ECALnCells;
  std::vector<float>            recojet_HCALtime;
  std::vector<float>            recojet_HCALenergy;
  std::vector<float>            recojet_HCALnCells;
  std::vector<float>            recojet_HCALTDCtime;
  std::vector<float>            recojet_HCALTDCenergy;
  std::vector<float>            recojet_HCALTDCnCells;
  std::vector<float>            recojet_MTDtime;
  std::vector<float>            recojet_MTDenergy;
  std::vector<float>            recojet_MTDnCells;
  std::vector<float>            recojet_MTDClutime;
  std::vector<float>            recojet_MTDCluenergy;
  std::vector<float>            recojet_MTDnClus;
  std::vector<float>            recojet_HGCALtime;
  std::vector<float>            recojet_HGCALenergy;
  std::vector<float>            recojet_HGCALnTracksters;
  std::vector<float>  recojet_closestGenIndex;
  std::vector<float>  recojet_closestGenR;
  std::vector<float>  recojet_closestEbGenIndex;
  std::vector<float>  recojet_closestEbGenR;
  std::vector<float>  recojet_closestEtlGenIndex;
  std::vector<float>  recojet_closestEtlGenR;
  std::vector<float>  recojet_closestHgGenIndex;
  std::vector<float>  recojet_closestHgGenR;
  int                           nmuons;
  std::vector<float>            muon_pt;
  std::vector<bool> 		muon_tightId;
  std::vector<bool> 		muon_looseId;
  std::vector<bool> 		muon_tightIso;
  std::vector<bool> 		muon_looseIso;
  std::vector<float>            muon_eta;
  std::vector<float>            muon_phi;
  std::vector<float>            muon_e;
  std::vector<float>            muon_ECALtime;
  std::vector<float>            muon_ECALenergy;
  std::vector<float>            muon_ECALnCells;
  std::vector<float>            muon_HCALtime;
  std::vector<float>            muon_HCALenergy;
  std::vector<float>            muon_HCALnCells;
  std::vector<float>            muon_HCALTDCtime;
  std::vector<float>            muon_HCALTDCenergy;
  std::vector<float>            muon_HCALTDCnCells;
  std::vector<float>            muon_MTDtime;
  std::vector<float>            muon_MTDenergy;
  std::vector<float>            muon_MTDnCells;
  std::vector<float>            muon_MTDClutime;
  std::vector<float>            muon_MTDCluenergy;
  std::vector<float>            muon_MTDnClus;
  std::vector<float>            muon_HGCALtime;
  std::vector<float>            muon_HGCALenergy;
  std::vector<float>            muon_HGCALnTracksters;
  std::vector<float>  muon_closestGenIndex;
  std::vector<float>  muon_closestGenR;
  std::vector<float>  muon_closestEbGenIndex;
  std::vector<float>  muon_closestEbGenR;
  std::vector<float>  muon_closestEtlGenIndex;
  std::vector<float>  muon_closestEtlGenR;
  std::vector<float>  muon_closestHgGenIndex;
  std::vector<float>  muon_closestHgGenR;
  int                           nelectrons;
  std::vector<float>            electron_pt;
  std::vector<float>            electron_eta;
  std::vector<float>            electron_phi;
  std::vector<float>            electron_e;
  std::vector<float>            electron_ECALtime;
  std::vector<float>            electron_ECALenergy;
  std::vector<float>            electron_ECALnCells;
  std::vector<float>            electron_HCALtime;
  std::vector<float>            electron_HCALenergy;
  std::vector<float>            electron_HCALnCells;
  std::vector<float>            electron_HCALTDCtime;
  std::vector<float>            electron_HCALTDCenergy;
  std::vector<float>            electron_HCALTDCnCells;
  std::vector<float>            electron_MTDtime;
  std::vector<float>            electron_MTDenergy;
  std::vector<float>            electron_MTDnCells;
  std::vector<float>            electron_MTDClutime;
  std::vector<float>            electron_MTDCluenergy;
  std::vector<float>            electron_MTDnClus;
  std::vector<float>            electron_HGCALtime;
  std::vector<float>            electron_HGCALenergy;
  std::vector<float>            electron_HGCALnTracksters;
  std::vector<float>  electron_closestGenIndex;
  std::vector<float>  electron_closestGenR;
  std::vector<float>  electron_closestEbGenIndex;
  std::vector<float>  electron_closestEbGenR;
  std::vector<float>  electron_closestEtlGenIndex;
  std::vector<float>  electron_closestEtlGenR;
  std::vector<float>  electron_closestHgGenIndex;
  std::vector<float>  electron_closestHgGenR;
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
  edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronInputToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonInputToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexCollectionToken_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidatesToken_;
  edm::Service<TFileService> fs;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > _genParticles; 
  edm::Handle< edm::View<reco::GenParticle> > _genParticlesH;
  const edm::EDGetTokenT< edm::View<reco::PFJet> > _recoak4PFJets; 
  edm::Handle< edm::View<reco::PFJet> > _recoak4PFJetsH;
  const edm::EDGetTokenT<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEBToken_;
  edm::Handle< edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>> > _ecalRecHitsEBH;
  const edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHitsToken_;
  edm::Handle< edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>> > _hcalRecHitsH;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterEMToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterMergeToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterHADToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterTrkEMToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterTrkToken;
  const edm::EDGetTokenT<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>> mtdRecHitsBTLToken_;
  edm::Handle< edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>> > _mtdRecHitsBTLH;
  const edm::EDGetTokenT<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>> mtdRecHitsETLToken_;
  edm::Handle< edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>> > _mtdRecHitsETLH;
  const edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::Handle<FTLClusterCollection> _btlRecCluH;
  const edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;
  edm::Handle<FTLClusterCollection> _etlRecCluH;
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
   electronInputToken_(consumes<std::vector<reco::GsfElectron>>(edm::InputTag("gedGsfElectrons"))),
   muonInputToken_(consumes<std::vector<reco::Muon>>(edm::InputTag("muons"))),
   vertexCollectionToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"))),
   pfCandidatesToken_(consumes<std::vector<reco::PFCandidate>>(edm::InputTag("pfCandidates"))),
  _genParticles(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  _genParticlesH(),
  _recoak4PFJets(consumes< edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("recoak4PFJets"))),
  _recoak4PFJetsH(),
  ecalRecHitsEBToken_{consumes<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>>(												       iConfig.getParameter<edm::InputTag>("ebRecHitsColl"))},
  _ecalRecHitsEBH(),
  hcalRecHitsToken_{consumes<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>>(												       iConfig.getParameter<edm::InputTag>("hcalRecHitsColl"))},
  _hcalRecHitsH(),
  _tracksterEMToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersEM"))),
  _tracksterMergeToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersMerge"))),
  _tracksterHADToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersHAD"))),
  _tracksterTrkEMToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersTrkEM"))),
  _tracksterTrkToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersTrk"))),
  mtdRecHitsBTLToken_{consumes<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>>(												       iConfig.getParameter<edm::InputTag>("mtdBTLRecHitsColl"))},
  _mtdRecHitsBTLH(),
  mtdRecHitsETLToken_{consumes<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>>(												       iConfig.getParameter<edm::InputTag>("mtdETLRecHitsColl"))},
  _mtdRecHitsETLH(),
  btlRecCluToken_(consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recBTLCluTag"))),
  _btlRecCluH(),
  etlRecCluToken_(consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recETLCluTag"))),
  _etlRecCluH()
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

   _jetTimingTools.init(iSetup);
  Handle< std::vector<reco::Vertex> > vertexCollectionH;
  Handle<std::vector<reco::PFCandidate>> pfCandidatesH;
  Handle<View<ticl::Trackster>> tracksterEMH;
  Handle<View<ticl::Trackster>> tracksterMergeH;
  Handle<View<ticl::Trackster>> tracksterHADH;
  Handle<View<ticl::Trackster>> tracksterTrkEMH;
  Handle<View<ticl::Trackster>> tracksterTrkH;
  Handle<std::vector<reco::Muon>> muonsH;
  Handle<std::vector<reco::GsfElectron>> electronsH;

  iEvent.getByToken(vertexCollectionToken_,vertexCollectionH);
  iEvent.getByToken(electronInputToken_, electronsH); 
  iEvent.getByToken(muonInputToken_, muonsH); 
  iEvent.getByToken(pfCandidatesToken_,pfCandidatesH);
  iEvent.getByToken(_genParticles, _genParticlesH);
  iEvent.getByToken(_recoak4PFJets, _recoak4PFJetsH);
  iEvent.getByToken(ecalRecHitsEBToken_, _ecalRecHitsEBH);
  iEvent.getByToken(_tracksterEMToken, tracksterEMH);
  iEvent.getByToken(_tracksterMergeToken, tracksterMergeH);
  iEvent.getByToken(_tracksterHADToken, tracksterHADH);
  iEvent.getByToken(_tracksterTrkEMToken, tracksterTrkEMH);
  iEvent.getByToken(_tracksterTrkToken, tracksterTrkH);
  iEvent.getByToken(mtdRecHitsBTLToken_, _mtdRecHitsBTLH);
  iEvent.getByToken(mtdRecHitsETLToken_, _mtdRecHitsETLH);
  
  iEvent.getByToken(btlRecCluToken_,_btlRecCluH);
  iEvent.getByToken(etlRecCluToken_,_etlRecCluH);
  //  auto _btlRecCluH = makeValid(iEvent.getHandle(btlRecCluToken_));
  // auto _etlRecCluH = makeValid(iEvent.getHandle(etlRecCluToken_));

  //  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);

  //variable declaration
  int interestingquark = 0;
  int nrecojets = 0;
  int nmuons = 0;
  int nelectrons = 0;
  int ngen = 0;
  int npvtrack = 0;
  int npfCand = 0;
  // std::vector<float> tracksterEM_time;
  // std::vector<float> tracksterHAD_time;
  // std::vector<float> tracksterTrkEM_time;
  std::vector<float> tracksterMerge_time;
  std::vector<int> tracksterMerge_ticlIteration;
  // std::vector<float> tracksterEM_timeError;
  // std::vector<float> tracksterHAD_timeError;
  // std::vector<float> tracksterTrkEM_timeError;
  std::vector<float> tracksterMerge_timeError;
  // std::vector<float> tracksterEM_e;
  // std::vector<float> tracksterHAD_e;
  // std::vector<float> tracksterTrkEM_e;
  std::vector<float> tracksterMerge_e;
  // std::vector<float> tracksterEM_eta;
  // std::vector<float> tracksterHAD_eta;
  // std::vector<float> tracksterTrkEM_eta;
  std::vector<float> tracksterMerge_eta;
  // std::vector<float> tracksterEM_phi;
  // std::vector<float> tracksterHAD_phi;
  // std::vector<float> tracksterTrkEM_phi;
  std::vector<float> tracksterMerge_phi;
  std::vector<uint> tracksterMerge_iJ;
  std::vector<float> pfCand_eta;
  std::vector<float> pfCand_phi;
  std::vector<int> pfCand_iJ;
  std::vector<float> pfCand_pt;
  std::vector<float> pfCand_time;
  std::vector<float> pfCand_timeError;
  std::vector<float> pvtrack_eta;
  std::vector<float> pvtrack_phi;
  std::vector<float> pvtrack_pt;
  std::vector<float> q_eta;
  std::vector<float> q_phi;
  std::vector<float> q_pt;
  std::vector<int> q_pdgId;
  std::vector<float> q_ctau;
  std::vector<float> q_vx;
  std::vector<float> q_vy;
  std::vector<float> q_vz;
  std::vector<float> q_ebeta;
  std::vector<float> q_ebphi;
  std::vector<float> q_ebdelay;
  std::vector<float> q_etleta;
  std::vector<float> q_etlphi;
  std::vector<float> q_etldelay;
  std::vector<float> q_hgeta;
std::vector<float> q_pathdelay;
std::vector<bool> q_decayInHGCAL;
  std::vector<float> q_hgphi;
  std::vector<float> q_hgdelay;
  std::vector<float>    genjet_pt;
  std::vector<float>    genjet_eta;
  std::vector<float>    genjet_phi;
  std::vector<float>    recojet_pt;
  std::vector<float>    recojet_eta;
  std::vector<float>    recojet_phi;
  std::vector<float>    recojet_e;
  std::vector<float>    recojet_ECALtime;
  std::vector<float>    recojet_ECALenergy;
  std::vector<float>    recojet_ECALnCells;
  std::vector<float>    recojet_HCALtime;
  std::vector<float>    recojet_HCALenergy;
  std::vector<float>    recojet_HCALnCells;
  std::vector<float>    recojet_HCALTDCtime;
  std::vector<float>    recojet_HCALTDCenergy;
  std::vector<float>    recojet_HCALTDCnCells;
  std::vector<float>    recojet_MTDtime;
  std::vector<float>    recojet_MTDenergy;
  std::vector<float>    recojet_MTDnCells;
  std::vector<float>            recojet_MTDClutime;
  std::vector<float>            recojet_MTDCluenergy;
  std::vector<float>            recojet_MTDnClus;
  std::vector<float>    recojet_HGCALtime;
  std::vector<float>    recojet_HGCALenergy;
  std::vector<float>    recojet_HGCALnTracksters;
  std::vector<float>  recojet_closestGenIndex;
  std::vector<float>  recojet_closestGenR;
  std::vector<float>  recojet_closestEbGenIndex;
  std::vector<float>  recojet_closestEbGenR;
  std::vector<float>  recojet_closestEtlGenIndex;
  std::vector<float>  recojet_closestEtlGenR;
  std::vector<float>  recojet_closestHgGenIndex;
  std::vector<float>  recojet_closestHgGenR;
  std::vector<float>    muon_pt;
    std::vector<bool> muon_tightId;
    std::vector<bool> muon_looseId;
    std::vector<bool> muon_tightIso;
    std::vector<bool> muon_looseIso;
  std::vector<float>    muon_eta;
  std::vector<float>    muon_phi;
  std::vector<float>    muon_e;
  std::vector<float>    muon_ECALtime;
  std::vector<float>    muon_ECALenergy;
  std::vector<float>    muon_ECALnCells;
  std::vector<float>    muon_HCALtime;
  std::vector<float>    muon_HCALenergy;
  std::vector<float>    muon_HCALnCells;
  std::vector<float>    muon_HCALTDCtime;
  std::vector<float>    muon_HCALTDCenergy;
  std::vector<float>    muon_HCALTDCnCells;
  std::vector<float>    muon_MTDtime;
  std::vector<float>    muon_MTDenergy;
  std::vector<float>    muon_MTDnCells;
  std::vector<float>            muon_MTDClutime;
  std::vector<float>            muon_MTDCluenergy;
  std::vector<float>            muon_MTDnClus;
  std::vector<float>    muon_HGCALtime;
  std::vector<float>    muon_HGCALenergy;
  std::vector<float>    muon_HGCALnTracksters;
  std::vector<float>  muon_closestGenIndex;
  std::vector<float>  muon_closestGenR;
  std::vector<float>  muon_closestEbGenIndex;
  std::vector<float>  muon_closestEbGenR;
  std::vector<float>  muon_closestEtlGenIndex;
  std::vector<float>  muon_closestEtlGenR;
  std::vector<float>  muon_closestHgGenIndex;
  std::vector<float>  muon_closestHgGenR;
  std::vector<float>    electron_pt;
  std::vector<float>    electron_eta;
  std::vector<float>    electron_phi;
  std::vector<float>    electron_e;
  std::vector<float>    electron_ECALtime;
  std::vector<float>    electron_ECALenergy;
  std::vector<float>    electron_ECALnCells;
  std::vector<float>    electron_HCALtime;
  std::vector<float>    electron_HCALenergy;
  std::vector<float>    electron_HCALnCells;
  std::vector<float>    electron_HCALTDCtime;
  std::vector<float>    electron_HCALTDCenergy;
  std::vector<float>    electron_HCALTDCnCells;
  std::vector<float>    electron_MTDtime;
  std::vector<float>    electron_MTDenergy;
  std::vector<float>    electron_MTDnCells;
  std::vector<float>            electron_MTDClutime;
  std::vector<float>            electron_MTDCluenergy;
  std::vector<float>            electron_MTDnClus;
  std::vector<float>    electron_HGCALtime;
  std::vector<float>    electron_HGCALenergy;
  std::vector<float>    electron_HGCALnTracksters;
  std::vector<float>  electron_closestGenIndex;
  std::vector<float>  electron_closestGenR;
  std::vector<float>  electron_closestEbGenIndex;
  std::vector<float>  electron_closestEbGenR;
  std::vector<float>  electron_closestEtlGenIndex;
  std::vector<float>  electron_closestEtlGenR;
  std::vector<float>  electron_closestHgGenIndex;
  std::vector<float>  electron_closestHgGenR;

  
  bool debug=0;


  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON GENPARTICLES --------------------------------------"<<std::endl; 
  double xArr[11] = {
    302.27485338939545,
318.02128494261274,
317.9398066651239,
382.31971673724706,
452.5378974580781,
524.0255298602788,
524.8161709233189,
519.1267741643437,
520.6084717290485,
318.5222254634702,
302.2919537686215
};
    double yArr[11] = {
     130.99494029955773,
     140.8085461659938,
    159.38559343345372,
     180.76609698932697,
    271.0208926397953,
    271.84070493798606,
    91.57454256485633,
    88.75700361119829,
    50.92995885850013,
    26.59410741049942,
    27.096053836017063
    };

  TGeoPolygon hgcal(11);
  hgcal.SetXY(xArr,yArr);
  hgcal.FinishPolygon();
  for (const auto & genpar_iter : *_genParticlesH){

    if (genpar_iter.mother(0) == NULL)continue;
    // if (genpar_iter.mother(0)->pdgId() == 6000113){
    // std::cout << genpar_iter.pdgId() << " " << genpar_iter.status() << " " << genpar_iter.mother(0)->pdgId() << std::endl;
    // }
    if(abs(genpar_iter.pdgId()) > 16  || (genpar_iter.status() != 23 && genpar_iter.status() != 1) || genpar_iter.mother(0)->pdgId() != 6000113)continue;
    float vx = genpar_iter.vertex().x();
    float vy = genpar_iter.vertex().y();
    float vz = genpar_iter.vertex().z();
    double vertexPos[2] = {fabs(genpar_iter.vertex().z()),TMath::Sqrt(genpar_iter.vertex().x()*genpar_iter.vertex().x()+genpar_iter.vertex().y()*genpar_iter.vertex().y())};
    bool decayInHGCAL = hgcal.Contains(vertexPos);
    interestingquark++;
    reco::GenParticle * genParticleMother = (reco::GenParticle *) genpar_iter.mother();
    std::vector<double> ecalIntersection = _jetTimingTools.surfaceIntersection(genpar_iter,*genParticleMother,130);
    std::vector<double> etlIntersection = _jetTimingTools.endCapIntersection(genpar_iter,*genParticleMother,300,300);
    std::vector<double> hgcalIntersection = _jetTimingTools.endCapIntersection(genpar_iter,*genParticleMother,300,520);
    ngen++;
    q_pt.push_back(genpar_iter.pt());
    q_eta.push_back(genpar_iter.eta());
    q_phi.push_back(genpar_iter.phi());
    double displacement = TMath::Sqrt((genpar_iter.vertex()-genParticleMother->vertex()).Mag2());
    double genParticleBeta = genParticleMother->p()/genParticleMother->energy();
    double genParticleGamma = 1./TMath::Sqrt(1.-genParticleBeta*genParticleBeta);
    double ctau = displacement*10 / (genParticleBeta*genParticleGamma);
    q_ctau.push_back(ctau);
    q_vx.push_back(vx);
    q_vy.push_back(vy);
    q_vz.push_back(vz);
    q_pdgId.push_back(genpar_iter.pdgId());
    q_ebphi.push_back(ecalIntersection[1]);
    q_ebeta.push_back(ecalIntersection[0]);
    q_ebdelay.push_back(ecalIntersection[3]*1E9);
    q_etlphi.push_back(etlIntersection[1]);
    q_etleta.push_back(etlIntersection[0]);
    q_etldelay.push_back(etlIntersection[3]*1E9);
    q_hgphi.push_back(hgcalIntersection[1]);
    q_hgeta.push_back(hgcalIntersection[0]);
    q_pathdelay.push_back(hgcalIntersection[5]*1E9);
    q_hgdelay.push_back(hgcalIntersection[3]*1E9);
    q_decayInHGCAL.push_back(decayInHGCAL);
  }

  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);
  auto const& hcalRecHits = iEvent.get(hcalRecHitsToken_);
  auto const& mtdRecHitsBTL = iEvent.get(mtdRecHitsBTLToken_);
  auto const& mtdRecHitsETL = iEvent.get(mtdRecHitsETLToken_);
  //  auto const& mtdClusBTL = iEvent.get(btlRecCluToken_);
  //  auto const& mtdClusETL = iEvent.get(btlRecCluToken_);
  //
  reco::Vertex primaryVertex = vertexCollectionH->at(0);
  for(auto pvTrack=primaryVertex.tracks_begin(); pvTrack!=primaryVertex.tracks_end(); pvTrack++){
      pvtrack_pt.push_back((*pvTrack)->pt());
      pvtrack_eta.push_back((*pvTrack)->eta());
      pvtrack_phi.push_back((*pvTrack)->phi());
      npvtrack++;
  }

  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON RECO JETS --------------------------------------"<<std::endl; 
  std::vector<ticl::Trackster> tracksters;
  tracksters.reserve(tracksters.size() + distance(tracksterMergeH->begin(),tracksterMergeH->end()));
  tracksters.insert(tracksters.end(),tracksterMergeH->begin(),tracksterMergeH->end());
    for (auto const & trackster : *tracksterMergeH){
      if (trackster.regressed_energy() < 1.) continue;
      tracksterMerge_time.push_back(trackster.time());
      tracksterMerge_ticlIteration.push_back(trackster.ticlIteration());
      tracksterMerge_timeError.push_back(trackster.timeError());
      tracksterMerge_e.push_back(trackster.regressed_energy());
      tracksterMerge_eta.push_back(trackster.barycenter().eta());
      tracksterMerge_phi.push_back(trackster.barycenter().phi());
      uint matchI = 999;
      double minMatch = 999;
      uint iJ = -1;
      for (const auto & recojet_iter : *_recoak4PFJetsH){

	if(recojet_iter.pt()<10)continue;
	if(fabs(recojet_iter.eta())>3)continue;
	iJ +=1;
	auto const pos = trackster.barycenter();
	double matchR = reco::deltaR2(recojet_iter, pos);
	if (matchR > _jetTimingTools.getMatchingRadius()) continue;
        if (matchR < minMatch) {matchI = iJ;minMatch=matchR;}
      }
      tracksterMerge_iJ.push_back(matchI);
  }
  for (const auto & electron_iter : *electronsH){

    if(electron_iter.pt()<10)continue;
    if(fabs(electron_iter.eta())>3)continue;
    
    nelectrons++;
    electron_pt.push_back(electron_iter.pt());
    electron_eta.push_back(electron_iter.eta());
    electron_phi.push_back(electron_iter.phi());
    electron_e.push_back(electron_iter.energy());
    int closestGenIndex = -1;
    float closestGenR = 999;
    int closestEbGenIndex = -1;
    float closestEbGenR = 999;
    int closestEtlGenIndex = -1;
    float closestEtlGenR = 999;
    int closestHgGenIndex = -1;
    float closestHgGenR = 999;
    for (int ig = 0; ig < ngen; ig++){
	TLorentzVector jetVec;
	jetVec.SetPtEtaPhiM(electron_iter.pt(),electron_iter.eta(),electron_iter.phi(),0);
	TLorentzVector pos;
	pos.SetPtEtaPhiM(1,q_eta[ig],q_phi[ig],0);
	TLorentzVector ebPos;
	ebPos.SetPtEtaPhiM(1,q_ebeta[ig],q_ebphi[ig],0);
	TLorentzVector hgPos;
	hgPos.SetPtEtaPhiM(1,q_hgeta[ig],q_hgphi[ig],0);
	TLorentzVector etlPos;
	etlPos.SetPtEtaPhiM(1,q_etleta[ig],q_etlphi[ig],0);
	
	float ebDeltaR = jetVec.DeltaR(ebPos);
	if (ebDeltaR < 0.4 and closestEbGenR > ebDeltaR){
	    closestEbGenR = ebDeltaR;
	    closestEbGenIndex = ig;
	}
	float etlDeltaR = jetVec.DeltaR(etlPos);
	if (etlDeltaR < 0.4 and closestEtlGenR > etlDeltaR){
	    closestEtlGenR = etlDeltaR;
	    closestEtlGenIndex = ig;
	}
	float hgDeltaR = jetVec.DeltaR(hgPos);
	if (hgDeltaR < 0.4 and closestHgGenR > hgDeltaR){
	    closestHgGenR = hgDeltaR;
	    closestHgGenIndex = ig;
	}
	float genDeltaR = jetVec.DeltaR(pos);
	if (genDeltaR < 0.4 and closestGenR > genDeltaR){
	    closestGenR = genDeltaR;
	    closestGenIndex = ig;
	}
    }
    electron_closestGenIndex.push_back(closestGenIndex);
    electron_closestGenR.push_back(closestGenR);
    electron_closestEbGenIndex.push_back(closestEbGenIndex);
    electron_closestEbGenR.push_back(closestEbGenR);
    electron_closestEtlGenIndex.push_back(closestEtlGenIndex);
    electron_closestEtlGenR.push_back(closestEtlGenR);
    electron_closestHgGenIndex.push_back(closestHgGenIndex);
    electron_closestHgGenR.push_back(closestHgGenR);
    

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM ECAL --------------------------------------"<<std::endl; 
    float weightedECALTimeCell = 0;
    float totalECALEnergyCell = 0;
    unsigned int ECALnCells = 0;
    if(fabs(electron_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromEcalCells(electron_iter, ecalRecHitsEB, weightedECALTimeCell, totalECALEnergyCell, ECALnCells);
    electron_ECALenergy.push_back(totalECALEnergyCell);
    electron_ECALnCells.push_back(ECALnCells);
    if(ECALnCells>0)
      electron_ECALtime.push_back(weightedECALTimeCell);
    else
      electron_ECALtime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL MAHI --------------------------------------"<<std::endl; 
    float weightedHCALTimeCell = 0;
    float totalHCALEnergyCell = 0;
    unsigned int HCALnCells = 0;
    if(fabs(electron_iter.eta())<1.4442)
    _jetTimingTools.jetTimeFromHcalCells(electron_iter, hcalRecHits, weightedHCALTimeCell, totalHCALEnergyCell, HCALnCells,false);
    electron_HCALenergy.push_back(totalHCALEnergyCell);
    electron_HCALnCells.push_back(HCALnCells);
    if(HCALnCells>0)
      electron_HCALtime.push_back(weightedHCALTimeCell);
    else
      electron_HCALtime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL TDC --------------------------------------"<<std::endl; 
    float weightedHCALTDCTimeCell = 0;
    float totalHCALTDCEnergyCell = 0;
    unsigned int HCALTDCnCells = 0;
    _jetTimingTools.jetTimeFromHcalCells(electron_iter, hcalRecHits, weightedHCALTDCTimeCell, totalHCALTDCEnergyCell, HCALTDCnCells,true);
    electron_HCALTDCenergy.push_back(totalHCALTDCEnergyCell);
    electron_HCALTDCnCells.push_back(HCALTDCnCells);
    if(HCALTDCnCells>0)
      electron_HCALTDCtime.push_back(weightedHCALTDCTimeCell);
    else
      electron_HCALTDCtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CELLS--------------------------------------"<<std::endl; 
    float weightedMTDTimeCell = 0;
    float totalMTDEnergyCell = 0;
    unsigned int MTDnCells = 0;
    if(fabs(electron_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDCells(electron_iter, mtdRecHitsBTL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,1);
    else if(fabs(electron_iter.eta())>1.4442 && fabs(electron_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDCells(electron_iter, mtdRecHitsETL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,0);
    }

    electron_MTDenergy.push_back(totalMTDEnergyCell);
    electron_MTDnCells.push_back(MTDnCells);
    if(MTDnCells>0)
      electron_MTDtime.push_back(weightedMTDTimeCell);
    else
      electron_MTDtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CLUSTERS--------------------------------------"<<std::endl; 
    float weightedMTDTimeClu = 0;
    float totalMTDEnergyClu = 0;
    unsigned int MTDnClus = 0;
    if(fabs(electron_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDClus(electron_iter, _btlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,1);
    else if(fabs(electron_iter.eta())>1.4442 && fabs(electron_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDClus(electron_iter, _etlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,0);
    }
    electron_MTDCluenergy.push_back(totalMTDEnergyClu);
    electron_MTDnClus.push_back(MTDnClus);
    if(MTDnClus>0)
      electron_MTDClutime.push_back(weightedMTDTimeClu);
    else
      electron_MTDClutime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HGCAL --------------------------------------"<<std::endl; 
    float weightedHGCALTimeTrackster = 0;
    float totalHGCALEnergyTrackster = 0;
    unsigned int HGCALnTracksters = 0;
    _jetTimingTools.jetTimeFromHgcalTracksters(electron_iter, tracksters, weightedHGCALTimeTrackster, totalHGCALEnergyTrackster, HGCALnTracksters);
    electron_HGCALenergy.push_back(totalHGCALEnergyTrackster);
    electron_HGCALnTracksters.push_back(HGCALnTracksters);
    if(HGCALnTracksters>0)
      electron_HGCALtime.push_back(weightedHGCALTimeTrackster);
    else
      electron_HGCALtime.push_back(-50);
  }
  for (const auto & muon_iter : *muonsH){

    if(muon_iter.pt()<10)continue;
    if(fabs(muon_iter.eta())>3)continue;
    
    nmuons++;
    muon_pt.push_back(muon_iter.pt());
    bool passLooseId = false;
    bool passTightId = false;
    bool passLooseIso = false;
    bool passTightIso = false;
    if (muon_iter.isGlobalMuon() && muon_iter.isPFMuon()){
	passLooseId = true;
	if (muon_iter.globalTrack()->normalizedChi2() < 10. &&  muon_iter.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && fabs(muon_iter.muonBestTrack()->dxy(primaryVertex.position())) < 0.2 &&  fabs(muon_iter.muonBestTrack()->dz(primaryVertex.position())) < 0.5 &&  muon_iter.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && muon_iter.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)
	{
	    passTightId = true;
	}
    }
    double pfIso = (muon_iter.pfIsolationR04().sumChargedHadronPt + std::max(0., muon_iter.pfIsolationR04().sumNeutralHadronEt + muon_iter.pfIsolationR04().sumPhotonEt - 0.5*muon_iter.pfIsolationR04().sumPUPt))/muon_iter.pt();
    double tkIso = muon_iter.isolationR03().sumPt/muon_iter.pt();
    if (pfIso < 0.25) passLooseIso = true;
    if (pfIso < 0.15) passTightIso = true;
    muon_tightId.push_back(passTightId);
    muon_looseId.push_back(passLooseId);
    muon_tightIso.push_back(passTightIso);
    muon_looseIso.push_back(passLooseIso);
    muon_eta.push_back(muon_iter.eta());
    muon_phi.push_back(muon_iter.phi());
    muon_e.push_back(muon_iter.energy());
    int closestGenIndex = -1;
    float closestGenR = 999;
    int closestEbGenIndex = -1;
    float closestEbGenR = 999;
    int closestEtlGenIndex = -1;
    float closestEtlGenR = 999;
    int closestHgGenIndex = -1;
    float closestHgGenR = 999;
    for (int ig = 0; ig < ngen; ig++){
	TLorentzVector jetVec;
	jetVec.SetPtEtaPhiM(muon_iter.pt(),muon_iter.eta(),muon_iter.phi(),0);
	TLorentzVector pos;
	pos.SetPtEtaPhiM(1,q_eta[ig],q_phi[ig],0);
	TLorentzVector ebPos;
	ebPos.SetPtEtaPhiM(1,q_ebeta[ig],q_ebphi[ig],0);
	TLorentzVector hgPos;
	hgPos.SetPtEtaPhiM(1,q_hgeta[ig],q_hgphi[ig],0);
	TLorentzVector etlPos;
	etlPos.SetPtEtaPhiM(1,q_etleta[ig],q_etlphi[ig],0);
	
	float ebDeltaR = jetVec.DeltaR(ebPos);
	if (ebDeltaR < 0.4 and closestEbGenR > ebDeltaR){
	    closestEbGenR = ebDeltaR;
	    closestEbGenIndex = ig;
	}
	float etlDeltaR = jetVec.DeltaR(etlPos);
	if (etlDeltaR < 0.4 and closestEtlGenR > etlDeltaR){
	    closestEtlGenR = etlDeltaR;
	    closestEtlGenIndex = ig;
	}
	float hgDeltaR = jetVec.DeltaR(hgPos);
	if (hgDeltaR < 0.4 and closestHgGenR > hgDeltaR){
	    closestHgGenR = hgDeltaR;
	    closestHgGenIndex = ig;
	}
	float genDeltaR = jetVec.DeltaR(pos);
	if (genDeltaR < 0.4 and closestGenR > genDeltaR){
	    closestGenR = genDeltaR;
	    closestGenIndex = ig;
	}
    }
    muon_closestGenIndex.push_back(closestGenIndex);
    muon_closestGenR.push_back(closestGenR);
    muon_closestEbGenIndex.push_back(closestEbGenIndex);
    muon_closestEbGenR.push_back(closestEbGenR);
    muon_closestEtlGenIndex.push_back(closestEtlGenIndex);
    muon_closestEtlGenR.push_back(closestEtlGenR);
    muon_closestHgGenIndex.push_back(closestHgGenIndex);
    muon_closestHgGenR.push_back(closestHgGenR);
    

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM ECAL --------------------------------------"<<std::endl; 
    float weightedECALTimeCell = 0;
    float totalECALEnergyCell = 0;
    unsigned int ECALnCells = 0;
    if(fabs(muon_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromEcalCells(muon_iter, ecalRecHitsEB, weightedECALTimeCell, totalECALEnergyCell, ECALnCells);
    muon_ECALenergy.push_back(totalECALEnergyCell);
    muon_ECALnCells.push_back(ECALnCells);
    if(ECALnCells>0)
      muon_ECALtime.push_back(weightedECALTimeCell);
    else
      muon_ECALtime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL MAHI --------------------------------------"<<std::endl; 
    float weightedHCALTimeCell = 0;
    float totalHCALEnergyCell = 0;
    unsigned int HCALnCells = 0;
    if(fabs(muon_iter.eta())<1.4442)
    _jetTimingTools.jetTimeFromHcalCells(muon_iter, hcalRecHits, weightedHCALTimeCell, totalHCALEnergyCell, HCALnCells,false);
    muon_HCALenergy.push_back(totalHCALEnergyCell);
    muon_HCALnCells.push_back(HCALnCells);
    if(HCALnCells>0)
      muon_HCALtime.push_back(weightedHCALTimeCell);
    else
      muon_HCALtime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL TDC --------------------------------------"<<std::endl; 
    float weightedHCALTDCTimeCell = 0;
    float totalHCALTDCEnergyCell = 0;
    unsigned int HCALTDCnCells = 0;
    _jetTimingTools.jetTimeFromHcalCells(muon_iter, hcalRecHits, weightedHCALTDCTimeCell, totalHCALTDCEnergyCell, HCALTDCnCells,true);
    muon_HCALTDCenergy.push_back(totalHCALTDCEnergyCell);
    muon_HCALTDCnCells.push_back(HCALTDCnCells);
    if(HCALTDCnCells>0)
      muon_HCALTDCtime.push_back(weightedHCALTDCTimeCell);
    else
      muon_HCALTDCtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CELLS--------------------------------------"<<std::endl; 
    float weightedMTDTimeCell = 0;
    float totalMTDEnergyCell = 0;
    unsigned int MTDnCells = 0;
    if(fabs(muon_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDCells(muon_iter, mtdRecHitsBTL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,1);
    else if(fabs(muon_iter.eta())>1.4442 && fabs(muon_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDCells(muon_iter, mtdRecHitsETL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,0);
    }

    muon_MTDenergy.push_back(totalMTDEnergyCell);
    muon_MTDnCells.push_back(MTDnCells);
    if(MTDnCells>0)
      muon_MTDtime.push_back(weightedMTDTimeCell);
    else
      muon_MTDtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CLUSTERS--------------------------------------"<<std::endl; 
    float weightedMTDTimeClu = 0;
    float totalMTDEnergyClu = 0;
    unsigned int MTDnClus = 0;
    if(fabs(muon_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDClus(muon_iter, _btlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,1);
    else if(fabs(muon_iter.eta())>1.4442 && fabs(muon_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDClus(muon_iter, _etlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,0);
    }
    muon_MTDCluenergy.push_back(totalMTDEnergyClu);
    muon_MTDnClus.push_back(MTDnClus);
    if(MTDnClus>0)
      muon_MTDClutime.push_back(weightedMTDTimeClu);
    else
      muon_MTDClutime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HGCAL --------------------------------------"<<std::endl; 
    float weightedHGCALTimeTrackster = 0;
    float totalHGCALEnergyTrackster = 0;
    unsigned int HGCALnTracksters = 0;
    _jetTimingTools.jetTimeFromHgcalTracksters(muon_iter, tracksters, weightedHGCALTimeTrackster, totalHGCALEnergyTrackster, HGCALnTracksters);
    muon_HGCALenergy.push_back(totalHGCALEnergyTrackster);
    muon_HGCALnTracksters.push_back(HGCALnTracksters);
    if(HGCALnTracksters>0)
      muon_HGCALtime.push_back(weightedHGCALTimeTrackster);
    else
      muon_HGCALtime.push_back(-50);
  }
  uint iJ = -1;
  for (const auto & recojet_iter : *_recoak4PFJetsH){

    if(recojet_iter.pt()<10)continue;
    if(fabs(recojet_iter.eta())>3)continue;
    iJ++;
    
    nrecojets++;
    recojet_pt.push_back(recojet_iter.pt());
    recojet_eta.push_back(recojet_iter.eta());
    recojet_phi.push_back(recojet_iter.phi());
    recojet_e.push_back(recojet_iter.energy());
    std::vector<reco::PFCandidatePtr> pfCands = recojet_iter.getPFConstituents();
      for(reco::PFCandidatePtr& pfCand : pfCands){
	  if (pfCand->pt() < 2) continue;
	  pfCand_pt.push_back(pfCand->pt());
	  pfCand_time.push_back(pfCand->time());
	  pfCand_timeError.push_back(pfCand->timeError());
	  pfCand_eta.push_back(pfCand->eta());
	  pfCand_phi.push_back(pfCand->phi());
	  pfCand_iJ.push_back(iJ);
	  npfCand++;
      }
    int closestGenIndex = -1;
    float closestGenR = 999;
    int closestEbGenIndex = -1;
    float closestEbGenR = 999;
    int closestEtlGenIndex = -1;
    float closestEtlGenR = 999;
    int closestHgGenIndex = -1;
    float closestHgGenR = 999;
    for (int ig = 0; ig < ngen; ig++){
	TLorentzVector jetVec;
	jetVec.SetPtEtaPhiM(recojet_iter.pt(),recojet_iter.eta(),recojet_iter.phi(),0);
	TLorentzVector pos;
	pos.SetPtEtaPhiM(1,q_eta[ig],q_phi[ig],0);
	TLorentzVector ebPos;
	ebPos.SetPtEtaPhiM(1,q_ebeta[ig],q_ebphi[ig],0);
	TLorentzVector hgPos;
	hgPos.SetPtEtaPhiM(1,q_hgeta[ig],q_hgphi[ig],0);
	TLorentzVector etlPos;
	etlPos.SetPtEtaPhiM(1,q_etleta[ig],q_etlphi[ig],0);
	
	float ebDeltaR = jetVec.DeltaR(ebPos);
	if (ebDeltaR < 0.4 and closestEbGenR > ebDeltaR){
	    closestEbGenR = ebDeltaR;
	    closestEbGenIndex = ig;
	}
	float etlDeltaR = jetVec.DeltaR(etlPos);
	if (etlDeltaR < 0.4 and closestEtlGenR > etlDeltaR){
	    closestEtlGenR = etlDeltaR;
	    closestEtlGenIndex = ig;
	}
	float hgDeltaR = jetVec.DeltaR(hgPos);
	if (hgDeltaR < 0.4 and closestHgGenR > hgDeltaR){
	    closestHgGenR = hgDeltaR;
	    closestHgGenIndex = ig;
	}
	float genDeltaR = jetVec.DeltaR(pos);
	if (genDeltaR < 0.4 and closestGenR > genDeltaR){
	    closestGenR = genDeltaR;
	    closestGenIndex = ig;
	}
    }
    recojet_closestGenIndex.push_back(closestGenIndex);
    recojet_closestGenR.push_back(closestGenR);
    recojet_closestEbGenIndex.push_back(closestEbGenIndex);
    recojet_closestEbGenR.push_back(closestEbGenR);
    recojet_closestEtlGenIndex.push_back(closestEtlGenIndex);
    recojet_closestEtlGenR.push_back(closestEtlGenR);
    recojet_closestHgGenIndex.push_back(closestHgGenIndex);
    recojet_closestHgGenR.push_back(closestHgGenR);
    

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

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL MAHI --------------------------------------"<<std::endl; 
    float weightedHCALTimeCell = 0;
    float totalHCALEnergyCell = 0;
    unsigned int HCALnCells = 0;
    if(fabs(recojet_iter.eta())<1.4442)
    _jetTimingTools.jetTimeFromHcalCells(recojet_iter, hcalRecHits, weightedHCALTimeCell, totalHCALEnergyCell, HCALnCells,false);
    recojet_HCALenergy.push_back(totalHCALEnergyCell);
    recojet_HCALnCells.push_back(HCALnCells);
    if(HCALnCells>0)
      recojet_HCALtime.push_back(weightedHCALTimeCell);
    else
      recojet_HCALtime.push_back(-50);

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HCAL TDC --------------------------------------"<<std::endl; 
    float weightedHCALTDCTimeCell = 0;
    float totalHCALTDCEnergyCell = 0;
    unsigned int HCALTDCnCells = 0;
    _jetTimingTools.jetTimeFromHcalCells(recojet_iter, hcalRecHits, weightedHCALTDCTimeCell, totalHCALTDCEnergyCell, HCALTDCnCells,true);
    recojet_HCALTDCenergy.push_back(totalHCALTDCEnergyCell);
    recojet_HCALTDCnCells.push_back(HCALTDCnCells);
    if(HCALTDCnCells>0)
      recojet_HCALTDCtime.push_back(weightedHCALTDCTimeCell);
    else
      recojet_HCALTDCtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CELLS--------------------------------------"<<std::endl; 
    float weightedMTDTimeCell = 0;
    float totalMTDEnergyCell = 0;
    unsigned int MTDnCells = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDCells(recojet_iter, mtdRecHitsBTL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,1);
    else if(fabs(recojet_iter.eta())>1.4442 && fabs(recojet_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDCells(recojet_iter, mtdRecHitsETL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,0);
    }

    recojet_MTDenergy.push_back(totalMTDEnergyCell);
    recojet_MTDnCells.push_back(MTDnCells);
    if(MTDnCells>0)
      recojet_MTDtime.push_back(weightedMTDTimeCell);
    else
      recojet_MTDtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CLUSTERS--------------------------------------"<<std::endl; 
    float weightedMTDTimeClu = 0;
    float totalMTDEnergyClu = 0;
    unsigned int MTDnClus = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDClus(recojet_iter, _btlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,1);
    else if(fabs(recojet_iter.eta())>1.4442 && fabs(recojet_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDClus(recojet_iter, _etlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,0);
    }
    recojet_MTDCluenergy.push_back(totalMTDEnergyClu);
    recojet_MTDnClus.push_back(MTDnClus);
    if(MTDnClus>0)
      recojet_MTDClutime.push_back(weightedMTDTimeClu);
    else
      recojet_MTDClutime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HGCAL --------------------------------------"<<std::endl; 
    float weightedHGCALTimeTrackster = 0;
    float totalHGCALEnergyTrackster = 0;
    unsigned int HGCALnTracksters = 0;
  // Handle<View<ticl::Trackster>> tracksterEMH;
  // Handle<View<ticl::Trackster>> tracksterMergeH;
  // Handle<View<ticl::Trackster>> tracksterHADH;
  // Handle<View<ticl::Trackster>> tracksterTrkEMH;
  // Handle<View<ticl::Trackster>> tracksterTrkH;
    // for (auto const & trackster : *tracksterHADH){
    //   tracksterHAD_time.push_back(trackster.time());
    //   tracksterHAD_timeError.push_back(trackster.timeError());
    //   tracksterHAD_e.push_back(trackster.regressed_energy());
    //   tracksterHAD_eta.push_back(trackster.barycenter().eta());
    //   tracksterHAD_phi.push_back(trackster.barycenter().phi());
    // }
    // for (auto const & trackster : *tracksterTrkEMH){
    //   tracksterTrkEM_time.push_back(trackster.time());
    //   tracksterTrkEM_timeError.push_back(trackster.timeError());
    //   tracksterTrkEM_e.push_back(trackster.regressed_energy());
    //   tracksterTrkEM_eta.push_back(trackster.barycenter().eta());
    //   tracksterTrkEM_phi.push_back(trackster.barycenter().phi());
    // }
    // for (auto const & trackster : *tracksterTrkH){
    //   tracksterTrk_time.push_back(trackster.time());
    //   tracksterTrk_timeError.push_back(trackster.timeError());
    //   tracksterTrk_e.push_back(trackster.regressed_energy());
    //   tracksterTrk_eta.push_back(trackster.barycenter().eta());
    //   tracksterTrk_phi.push_back(trackster.barycenter().phi());
    // }
    // tracksters.reserve(tracksters.size() + distance(tracksterEMH->begin(),tracksterEMH->end()));
    // tracksters.insert(tracksters.end(),tracksterEMH->begin(),tracksterEMH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterHADH->begin(),tracksterHADH->end()));
    // tracksters.insert(tracksters.end(),tracksterHADH->begin(),tracksterHADH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterTrkEMH->begin(),tracksterTrkEMH->end()));
    // tracksters.insert(tracksters.end(),tracksterTrkEMH->begin(),tracksterTrkEMH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterTrkH->begin(),tracksterTrkH->end()));
    // tracksters.insert(tracksters.end(),tracksterTrkH->begin(),tracksterTrkH->end());
    _jetTimingTools.jetTimeFromHgcalTracksters(recojet_iter, tracksters, weightedHGCALTimeTrackster, totalHGCALEnergyTrackster, HGCALnTracksters);
    recojet_HGCALenergy.push_back(totalHGCALEnergyTrackster);
    recojet_HGCALnTracksters.push_back(HGCALnTracksters);
    if(HGCALnTracksters>0)
      recojet_HGCALtime.push_back(weightedHGCALTimeTrackster);
    else
      recojet_HGCALtime.push_back(-50);


  }
  
  

  // --- setup tree values                                                                                                                                   
  initTreeStructure();
  clearVectors();
  tree_.ngen        = ngen;
  tree_.npvtrack = npvtrack;
  tree_.npfCand = npfCand;
  // for (uint it = 0; it < tracksterEM_time.size(); it++){
  //   tree_.tracksterEM_time.push_back(tracksterEM_time[it]);
  //   tree_.tracksterEM_timeError.push_back(tracksterEM_timeError[it]);
  //   tree_.tracksterEM_e.push_back(tracksterEM_e[it]);
  //   tree_.tracksterEM_eta.push_back(tracksterEM_eta[it]);
  //   tree_.tracksterEM_phi.push_back(tracksterEM_phi[it]);
  // }
  // for (uint it = 0; it < tracksterHAD_time.size(); it++){
  //   tree_.tracksterHAD_time.push_back(tracksterHAD_time[it]);
  //   tree_.tracksterHAD_timeError.push_back(tracksterHAD_timeError[it]);
  //   tree_.tracksterHAD_e.push_back(tracksterHAD_e[it]);
  //   tree_.tracksterHAD_eta.push_back(tracksterHAD_eta[it]);
  //   tree_.tracksterHAD_phi.push_back(tracksterHAD_phi[it]);
  // }
  for (uint it = 0; it < tracksterMerge_time.size(); it++){
    tree_.tracksterMerge_ticlIteration.push_back(tracksterMerge_ticlIteration[it]);
    tree_.tracksterMerge_time.push_back(tracksterMerge_time[it]);
    tree_.tracksterMerge_timeError.push_back(tracksterMerge_timeError[it]);
    tree_.tracksterMerge_e.push_back(tracksterMerge_e[it]);
    tree_.tracksterMerge_eta.push_back(tracksterMerge_eta[it]);
    tree_.tracksterMerge_phi.push_back(tracksterMerge_phi[it]);
    tree_.tracksterMerge_iJ.push_back(tracksterMerge_iJ[it]);
  }
  // for (uint it = 0; it < tracksterTrkEM_time.size(); it++){
  //   tree_.tracksterTrkEM_time.push_back(tracksterTrkEM_time[it]);
  //   tree_.tracksterTrkEM_timeError.push_back(tracksterTrkEM_timeError[it]);
  //   tree_.tracksterTrkEM_e.push_back(tracksterTrkEM_e[it]);
  //   tree_.tracksterTrkEM_eta.push_back(tracksterTrkEM_eta[it]);
  //   tree_.tracksterTrkEM_phi.push_back(tracksterTrkEM_phi[it]);
  // }
  for (int it = 0; it < npfCand; it++){
    tree_.pfCand_pt.push_back(pfCand_pt[it]);
    tree_.pfCand_time.push_back(pfCand_time[it]);
    tree_.pfCand_timeError.push_back(pfCand_timeError[it]);
    tree_.pfCand_eta.push_back(pfCand_eta[it]);
    tree_.pfCand_phi.push_back(pfCand_phi[it]); 
    tree_.pfCand_iJ.push_back(pfCand_iJ[it]); 
  }
  for (int it = 0; it < npvtrack; it++){
    tree_.pvtrack_pt.push_back(pvtrack_pt[it]);
    tree_.pvtrack_eta.push_back(pvtrack_eta[it]);
    tree_.pvtrack_phi.push_back(pvtrack_phi[it]); 
  }
  for (int ig = 0; ig < ngen; ig++){
    tree_.q_pt.push_back(q_pt[ig]);
    tree_.q_eta.push_back(q_eta[ig]);
    tree_.q_phi.push_back(q_phi[ig]); 
    tree_.q_pdgId.push_back(q_pdgId[ig]);
    tree_.q_ctau.push_back(q_ctau[ig]);
    tree_.q_vx.push_back(q_vx[ig]);
    tree_.q_vy.push_back(q_vy[ig]);
    tree_.q_vz.push_back(q_vz[ig]); 
    tree_.q_ebphi.push_back(q_ebphi[ig]); 
    tree_.q_ebeta.push_back(q_ebeta[ig]); 
    tree_.q_ebdelay.push_back(q_ebdelay[ig]); 
    tree_.q_hgphi.push_back(q_hgphi[ig]); 
    tree_.q_hgeta.push_back(q_hgeta[ig]); 
    tree_.q_pathdelay.push_back(q_pathdelay[ig]);
    tree_.q_decayInHGCAL.push_back(q_decayInHGCAL[ig]);
    tree_.q_hgdelay.push_back(q_hgdelay[ig]); 
  }
  tree_.nrecojets        = nrecojets;
  for (int ij = 0; ij < nrecojets; ij++){
    tree_.recojet_closestGenIndex.push_back(recojet_closestGenIndex[ij]);
    tree_.recojet_closestGenR.push_back(recojet_closestGenR[ij]);
    tree_.recojet_closestEbGenIndex.push_back(recojet_closestEbGenIndex[ij]);
    tree_.recojet_closestEbGenR.push_back(recojet_closestEbGenR[ij]);
    tree_.recojet_closestEtlGenIndex.push_back(recojet_closestEtlGenIndex[ij]);
    tree_.recojet_closestEtlGenR.push_back(recojet_closestEtlGenR[ij]);
    tree_.recojet_closestHgGenIndex.push_back(recojet_closestHgGenIndex[ij]);
    tree_.recojet_closestHgGenR.push_back(recojet_closestHgGenR[ij]);
    tree_.recojet_pt.push_back(recojet_pt[ij]);
    tree_.recojet_eta.push_back(recojet_eta[ij]);
    tree_.recojet_phi.push_back(recojet_phi[ij]); 
    tree_.recojet_e.push_back(recojet_e[ij]);
    tree_.recojet_ECALtime.push_back(recojet_ECALtime[ij]);
    tree_.recojet_ECALenergy.push_back(recojet_ECALenergy[ij]);
    tree_.recojet_ECALnCells.push_back(recojet_ECALnCells[ij]);
    tree_.recojet_HCALtime.push_back(recojet_HCALtime[ij]);
    tree_.recojet_HCALenergy.push_back(recojet_HCALenergy[ij]);
    tree_.recojet_HCALnCells.push_back(recojet_HCALnCells[ij]);
    tree_.recojet_HCALTDCtime.push_back(recojet_HCALTDCtime[ij]);
    tree_.recojet_HCALTDCenergy.push_back(recojet_HCALTDCenergy[ij]);
    tree_.recojet_HCALTDCnCells.push_back(recojet_HCALTDCnCells[ij]);
    tree_.recojet_MTDtime.push_back(recojet_MTDtime[ij]);
    tree_.recojet_MTDenergy.push_back(recojet_MTDenergy[ij]);
    tree_.recojet_MTDnCells.push_back(recojet_MTDnCells[ij]);
    tree_.recojet_MTDClutime.push_back(recojet_MTDClutime[ij]);
    tree_.recojet_MTDCluenergy.push_back(recojet_MTDCluenergy[ij]);
    tree_.recojet_MTDnClus.push_back(recojet_MTDnClus[ij]);
    tree_.recojet_HGCALtime.push_back(recojet_HGCALtime[ij]);
    tree_.recojet_HGCALenergy.push_back(recojet_HGCALenergy[ij]);
    tree_.recojet_HGCALnTracksters.push_back(recojet_HGCALnTracksters[ij]);
  }
  tree_.nmuons        = nmuons;
  for (int ij = 0; ij < nmuons; ij++){
    tree_.muon_closestGenIndex.push_back(muon_closestGenIndex[ij]);
    tree_.muon_closestGenR.push_back(muon_closestGenR[ij]);
    tree_.muon_closestEbGenIndex.push_back(muon_closestEbGenIndex[ij]);
    tree_.muon_closestEbGenR.push_back(muon_closestEbGenR[ij]);
    tree_.muon_closestEtlGenIndex.push_back(muon_closestEtlGenIndex[ij]);
    tree_.muon_closestEtlGenR.push_back(muon_closestEtlGenR[ij]);
    tree_.muon_closestHgGenIndex.push_back(muon_closestHgGenIndex[ij]);
    tree_.muon_closestHgGenR.push_back(muon_closestHgGenR[ij]);
    tree_.muon_pt.push_back(muon_pt[ij]);
    tree_.muon_tightId.push_back(muon_tightId[ij]);
    tree_.muon_looseId.push_back(muon_looseId[ij]);
    tree_.muon_tightIso.push_back(muon_tightIso[ij]);
    tree_.muon_looseIso.push_back(muon_looseIso[ij]);
    tree_.muon_eta.push_back(muon_eta[ij]);
    tree_.muon_phi.push_back(muon_phi[ij]); 
    tree_.muon_e.push_back(muon_e[ij]);
    tree_.muon_ECALtime.push_back(muon_ECALtime[ij]);
    tree_.muon_ECALenergy.push_back(muon_ECALenergy[ij]);
    tree_.muon_ECALnCells.push_back(muon_ECALnCells[ij]);
    tree_.muon_HCALtime.push_back(muon_HCALtime[ij]);
    tree_.muon_HCALenergy.push_back(muon_HCALenergy[ij]);
    tree_.muon_HCALnCells.push_back(muon_HCALnCells[ij]);
    tree_.muon_HCALTDCtime.push_back(muon_HCALTDCtime[ij]);
    tree_.muon_HCALTDCenergy.push_back(muon_HCALTDCenergy[ij]);
    tree_.muon_HCALTDCnCells.push_back(muon_HCALTDCnCells[ij]);
    tree_.muon_MTDtime.push_back(muon_MTDtime[ij]);
    tree_.muon_MTDenergy.push_back(muon_MTDenergy[ij]);
    tree_.muon_MTDnCells.push_back(muon_MTDnCells[ij]);
    tree_.muon_MTDClutime.push_back(muon_MTDClutime[ij]);
    tree_.muon_MTDCluenergy.push_back(muon_MTDCluenergy[ij]);
    tree_.muon_MTDnClus.push_back(muon_MTDnClus[ij]);
    tree_.muon_HGCALtime.push_back(muon_HGCALtime[ij]);
    tree_.muon_HGCALenergy.push_back(muon_HGCALenergy[ij]);
    tree_.muon_HGCALnTracksters.push_back(muon_HGCALnTracksters[ij]);
  }
  tree_.nelectrons        = nelectrons;
  for (int ij = 0; ij < nelectrons; ij++){
    tree_.electron_closestGenIndex.push_back(electron_closestGenIndex[ij]);
    tree_.electron_closestGenR.push_back(electron_closestGenR[ij]);
    tree_.electron_closestEbGenIndex.push_back(electron_closestEbGenIndex[ij]);
    tree_.electron_closestEbGenR.push_back(electron_closestEbGenR[ij]);
    tree_.electron_closestEtlGenIndex.push_back(electron_closestEtlGenIndex[ij]);
    tree_.electron_closestEtlGenR.push_back(electron_closestEtlGenR[ij]);
    tree_.electron_closestHgGenIndex.push_back(electron_closestHgGenIndex[ij]);
    tree_.electron_closestHgGenR.push_back(electron_closestHgGenR[ij]);
    tree_.electron_pt.push_back(electron_pt[ij]);
    tree_.electron_eta.push_back(electron_eta[ij]);
    tree_.electron_phi.push_back(electron_phi[ij]); 
    tree_.electron_e.push_back(electron_e[ij]);
    tree_.electron_ECALtime.push_back(electron_ECALtime[ij]);
    tree_.electron_ECALenergy.push_back(electron_ECALenergy[ij]);
    tree_.electron_ECALnCells.push_back(electron_ECALnCells[ij]);
    tree_.electron_HCALtime.push_back(electron_HCALtime[ij]);
    tree_.electron_HCALenergy.push_back(electron_HCALenergy[ij]);
    tree_.electron_HCALnCells.push_back(electron_HCALnCells[ij]);
    tree_.electron_HCALTDCtime.push_back(electron_HCALTDCtime[ij]);
    tree_.electron_HCALTDCenergy.push_back(electron_HCALTDCenergy[ij]);
    tree_.electron_HCALTDCnCells.push_back(electron_HCALTDCnCells[ij]);
    tree_.electron_MTDtime.push_back(electron_MTDtime[ij]);
    tree_.electron_MTDenergy.push_back(electron_MTDenergy[ij]);
    tree_.electron_MTDnCells.push_back(electron_MTDnCells[ij]);
    tree_.electron_MTDClutime.push_back(electron_MTDClutime[ij]);
    tree_.electron_MTDCluenergy.push_back(electron_MTDCluenergy[ij]);
    tree_.electron_MTDnClus.push_back(electron_MTDnClus[ij]);
    tree_.electron_HGCALtime.push_back(electron_HGCALtime[ij]);
    tree_.electron_HGCALenergy.push_back(electron_HGCALenergy[ij]);
    tree_.electron_HGCALnTracksters.push_back(electron_HGCALnTracksters[ij]);
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
  // tree->Branch("tracksterEM_time",&tree_.tracksterEM_time);
  // tree->Branch("tracksterEM_timeError",&tree_.tracksterEM_timeError);
  // tree->Branch("tracksterEM_e",&tree_.tracksterEM_e);
  // tree->Branch("tracksterEM_eta",&tree_.tracksterEM_eta);
  // tree->Branch("tracksterEM_phi",&tree_.tracksterEM_phi);
  // tree->Branch("tracksterHAD_time",&tree_.tracksterHAD_time);
  // tree->Branch("tracksterHAD_timeError",&tree_.tracksterHAD_timeError);
  // tree->Branch("tracksterHAD_e",&tree_.tracksterHAD_e);
  // tree->Branch("tracksterHAD_eta",&tree_.tracksterHAD_eta);
  // tree->Branch("tracksterHAD_phi",&tree_.tracksterHAD_phi);
  tree->Branch("tracksterMerge_time",&tree_.tracksterMerge_time);
  tree->Branch("tracksterMerge_ticlIteration",&tree_.tracksterMerge_ticlIteration);
  tree->Branch("tracksterMerge_timeError",&tree_.tracksterMerge_timeError);
  tree->Branch("tracksterMerge_e",&tree_.tracksterMerge_e);
  tree->Branch("tracksterMerge_eta",&tree_.tracksterMerge_eta);
  tree->Branch("tracksterMerge_phi",&tree_.tracksterMerge_phi);
  tree->Branch("tracksterMerge_iJ",&tree_.tracksterMerge_iJ);
  // tree->Branch("tracksterTrkEM_time",&tree_.tracksterTrkEM_time);
  // tree->Branch("tracksterTrkEM_timeError",&tree_.tracksterTrkEM_timeError);
  // tree->Branch("tracksterTrkEM_e",&tree_.tracksterTrkEM_e);
  // tree->Branch("tracksterTrkEM_eta",&tree_.tracksterTrkEM_eta);
  // tree->Branch("tracksterTrkEM_phi",&tree_.tracksterTrkEM_phi);
  tree->Branch("npvtrack",              &tree_.npvtrack,                "npvtrack/I");
  tree->Branch("npfCand",              &tree_.npfCand,                "npfCand/I");
  tree->Branch("ngen",              &tree_.ngen,                "ngen/I");
  tree->Branch("pvtrack_eta", &tree_.pvtrack_eta);
  tree->Branch("pvtrack_phi", &tree_.pvtrack_phi);
  tree->Branch("pvtrack_pt", &tree_.pvtrack_pt);
  tree->Branch("pfCand_eta", &tree_.pfCand_eta);
  tree->Branch("pfCand_phi", &tree_.pfCand_phi);
  tree->Branch("pfCand_iJ", &tree_.pfCand_iJ);
  tree->Branch("pfCand_pt", &tree_.pfCand_pt);
  tree->Branch("pfCand_time", &tree_.pfCand_time);
  tree->Branch("pfCand_timeError", &tree_.pfCand_timeError);
  tree->Branch("q_eta", &tree_.q_eta);
  tree->Branch("q_phi", &tree_.q_phi);
  tree->Branch("q_pt", &tree_.q_pt);
  tree->Branch("q_ebeta", &tree_.q_ebeta);
  tree->Branch("q_ebphi", &tree_.q_ebphi);
  tree->Branch("q_ebdelay", &tree_.q_ebdelay);
  tree->Branch("q_hgeta", &tree_.q_hgeta);
  tree->Branch("q_decayInHGCAL", &tree_.q_decayInHGCAL);
  tree->Branch("q_pathdelay", &tree_.q_pathdelay);
  tree->Branch("q_hgphi", &tree_.q_hgphi);
  tree->Branch("q_hgdelay", &tree_.q_hgdelay);
  tree->Branch("q_pdgId", &tree_.q_pdgId);
  tree->Branch("q_ctau", &tree_.q_ctau);
  tree->Branch("q_vx", &tree_.q_vx);
  tree->Branch("q_vy", &tree_.q_vy);
  tree->Branch("q_vz", &tree_.q_vz);

  tree->Branch("nrecojets",              &tree_.nrecojets,                "nrecojets/I");
  tree->Branch("recoJet_pt",             &tree_.recojet_pt);
  tree->Branch("recoJet_eta",             &tree_.recojet_eta);
  tree->Branch("recoJet_phi",             &tree_.recojet_phi);
  tree->Branch("recoJet_e",             &tree_.recojet_e);
  tree->Branch("recoJet_ECALtime",             &tree_.recojet_ECALtime);
  tree->Branch("recoJet_ECALenergy",             &tree_.recojet_ECALenergy);
  tree->Branch("recoJet_ECALnCells",             &tree_.recojet_ECALnCells);
  tree->Branch("recoJet_HCALtime",             &tree_.recojet_HCALtime);
  tree->Branch("recoJet_HCALenergy",             &tree_.recojet_HCALenergy);
  tree->Branch("recoJet_HCALnCells",             &tree_.recojet_HCALnCells);
  tree->Branch("recoJet_HCALTDCtime",             &tree_.recojet_HCALTDCtime);
  tree->Branch("recoJet_HCALTDCenergy",             &tree_.recojet_HCALTDCenergy);
  tree->Branch("recoJet_HCALTDCnCells",             &tree_.recojet_HCALTDCnCells);
  tree->Branch("recoJet_MTDtime",             &tree_.recojet_MTDtime);
  tree->Branch("recoJet_MTDenergy",             &tree_.recojet_MTDenergy);
  tree->Branch("recoJet_MTDnCells",             &tree_.recojet_MTDnCells);
  tree->Branch("recoJet_MTDClutime",             &tree_.recojet_MTDClutime);
  tree->Branch("recoJet_MTDCluenergy",             &tree_.recojet_MTDCluenergy);
  tree->Branch("recoJet_MTDnClus",             &tree_.recojet_MTDnClus);
  tree->Branch("recoJet_HGCALtime",             &tree_.recojet_HGCALtime);
  tree->Branch("recoJet_HGCALenergy",             &tree_.recojet_HGCALenergy);
  tree->Branch("recoJet_HGCALnTracksters",             &tree_.recojet_HGCALnTracksters);
  tree->Branch("recoJet_closestGenIndex",&tree_.recojet_closestGenIndex);
  tree->Branch("recoJet_closestGenR",&tree_.recojet_closestGenR);
  tree->Branch("recoJet_closestEbGenIndex",&tree_.recojet_closestEbGenIndex);
  tree->Branch("recoJet_closestEbGenR",&tree_.recojet_closestEbGenR);
  tree->Branch("recoJet_closestEtlGenIndex",&tree_.recojet_closestEtlGenIndex);
  tree->Branch("recoJet_closestEtlGenR",&tree_.recojet_closestEtlGenR);
  tree->Branch("recoJet_closestHgGenIndex",&tree_.recojet_closestHgGenIndex);
  tree->Branch("recoJet_closestHgGenR",&tree_.recojet_closestHgGenR);
 
  tree->Branch("nmuons",              &tree_.nmuons,                "nmuons/I");
  tree->Branch("muon_pt",             &tree_.muon_pt);
  tree->Branch("muon_tightId",             &tree_.muon_tightId);
  tree->Branch("muon_tightIso",             &tree_.muon_tightIso);
  tree->Branch("muon_looseId",             &tree_.muon_looseId);
  tree->Branch("muon_looseIso",             &tree_.muon_looseIso);
  tree->Branch("muon_eta",             &tree_.muon_eta);
  tree->Branch("muon_phi",             &tree_.muon_phi);
  tree->Branch("muon_e",             &tree_.muon_e);
  tree->Branch("muon_ECALtime",             &tree_.muon_ECALtime);
  tree->Branch("muon_ECALenergy",             &tree_.muon_ECALenergy);
  tree->Branch("muon_ECALnCells",             &tree_.muon_ECALnCells);
  tree->Branch("muon_HCALtime",             &tree_.muon_HCALtime);
  tree->Branch("muon_HCALenergy",             &tree_.muon_HCALenergy);
  tree->Branch("muon_HCALnCells",             &tree_.muon_HCALnCells);
  tree->Branch("muon_HCALTDCtime",             &tree_.muon_HCALTDCtime);
  tree->Branch("muon_HCALTDCenergy",             &tree_.muon_HCALTDCenergy);
  tree->Branch("muon_HCALTDCnCells",             &tree_.muon_HCALTDCnCells);
  tree->Branch("muon_MTDtime",             &tree_.muon_MTDtime);
  tree->Branch("muon_MTDenergy",             &tree_.muon_MTDenergy);
  tree->Branch("muon_MTDnCells",             &tree_.muon_MTDnCells);
  tree->Branch("muon_MTDClutime",             &tree_.muon_MTDClutime);
  tree->Branch("muon_MTDCluenergy",             &tree_.muon_MTDCluenergy);
  tree->Branch("muon_MTDnClus",             &tree_.muon_MTDnClus);
  tree->Branch("muon_HGCALtime",             &tree_.muon_HGCALtime);
  tree->Branch("muon_HGCALenergy",             &tree_.muon_HGCALenergy);
  tree->Branch("muon_HGCALnTracksters",             &tree_.muon_HGCALnTracksters);
  tree->Branch("muon_closestGenIndex",&tree_.muon_closestGenIndex);
  tree->Branch("muon_closestGenR",&tree_.muon_closestGenR);
  tree->Branch("muon_closestEbGenIndex",&tree_.muon_closestEbGenIndex);
  tree->Branch("muon_closestEbGenR",&tree_.muon_closestEbGenR);
  tree->Branch("muon_closestEtlGenIndex",&tree_.muon_closestEtlGenIndex);
  tree->Branch("muon_closestEtlGenR",&tree_.muon_closestEtlGenR);
  tree->Branch("muon_closestHgGenIndex",&tree_.muon_closestHgGenIndex);
  tree->Branch("muon_closestHgGenR",&tree_.muon_closestHgGenR);
  tree->Branch("nelectrons",              &tree_.nelectrons,                "nelectrons/I");
  tree->Branch("electron_pt",             &tree_.electron_pt);
  tree->Branch("electron_eta",             &tree_.electron_eta);
  tree->Branch("electron_phi",             &tree_.electron_phi);
  tree->Branch("electron_e",             &tree_.electron_e);
  tree->Branch("electron_ECALtime",             &tree_.electron_ECALtime);
  tree->Branch("electron_ECALenergy",             &tree_.electron_ECALenergy);
  tree->Branch("electron_ECALnCells",             &tree_.electron_ECALnCells);
  tree->Branch("electron_HCALtime",             &tree_.electron_HCALtime);
  tree->Branch("electron_HCALenergy",             &tree_.electron_HCALenergy);
  tree->Branch("electron_HCALnCells",             &tree_.electron_HCALnCells);
  tree->Branch("electron_HCALTDCtime",             &tree_.electron_HCALTDCtime);
  tree->Branch("electron_HCALTDCenergy",             &tree_.electron_HCALTDCenergy);
  tree->Branch("electron_HCALTDCnCells",             &tree_.electron_HCALTDCnCells);
  tree->Branch("electron_MTDtime",             &tree_.electron_MTDtime);
  tree->Branch("electron_MTDenergy",             &tree_.electron_MTDenergy);
  tree->Branch("electron_MTDnCells",             &tree_.electron_MTDnCells);
  tree->Branch("electron_MTDClutime",             &tree_.electron_MTDClutime);
  tree->Branch("electron_MTDCluenergy",             &tree_.electron_MTDCluenergy);
  tree->Branch("electron_MTDnClus",             &tree_.electron_MTDnClus);
  tree->Branch("electron_HGCALtime",             &tree_.electron_HGCALtime);
  tree->Branch("electron_HGCALenergy",             &tree_.electron_HGCALenergy);
  tree->Branch("electron_HGCALnTracksters",             &tree_.electron_HGCALnTracksters);
  tree->Branch("electron_closestGenIndex",&tree_.electron_closestGenIndex);
  tree->Branch("electron_closestGenR",&tree_.electron_closestGenR);
  tree->Branch("electron_closestEbGenIndex",&tree_.electron_closestEbGenIndex);
  tree->Branch("electron_closestEbGenR",&tree_.electron_closestEbGenR);
  tree->Branch("electron_closestEtlGenIndex",&tree_.electron_closestEtlGenIndex);
  tree->Branch("electron_closestEtlGenR",&tree_.electron_closestEtlGenR);
  tree->Branch("electron_closestHgGenIndex",&tree_.electron_closestHgGenIndex);
  tree->Branch("electron_closestHgGenR",&tree_.electron_closestHgGenR);
}


// ------------ initialize trees ------------                                                                                                                 
void Phase2TimingAnalyzer::initTreeStructure()
{
}


void Phase2TimingAnalyzer::clearVectors()
{
  tree_.electron_pt.clear();
  tree_.electron_eta.clear();
  tree_.electron_phi.clear();
  tree_.electron_e.clear();
  tree_.electron_ECALtime.clear();
  tree_.electron_ECALenergy.clear();
  tree_.electron_ECALnCells.clear();
  tree_.electron_HCALtime.clear();
  tree_.electron_HCALenergy.clear();
  tree_.electron_HCALnCells.clear();
  tree_.electron_HCALTDCtime.clear();
  tree_.electron_HCALTDCenergy.clear();
  tree_.electron_HCALTDCnCells.clear();
  tree_.electron_MTDtime.clear();
  tree_.electron_MTDenergy.clear();
  tree_.electron_MTDnCells.clear();
  tree_.electron_MTDClutime.clear();
  tree_.electron_MTDCluenergy.clear();
  tree_.electron_MTDnClus.clear();
  tree_.electron_HGCALtime.clear();
  tree_.electron_HGCALenergy.clear();
  tree_.electron_HGCALnTracksters.clear();
  tree_.electron_closestGenIndex.clear();
  tree_.electron_closestGenR.clear();
  tree_.electron_closestEbGenIndex.clear();
  tree_.electron_closestEbGenR.clear();
  tree_.electron_closestEtlGenIndex.clear();
  tree_.electron_closestEtlGenR.clear();
  tree_.electron_closestHgGenIndex.clear();
  tree_.electron_closestHgGenR.clear();

  tree_.muon_pt.clear();
  tree_.muon_tightId.clear();
  tree_.muon_looseId.clear();
  tree_.muon_tightIso.clear();
  tree_.muon_looseIso.clear();
  tree_.muon_eta.clear();
  tree_.muon_phi.clear();
  tree_.muon_e.clear();
  tree_.muon_ECALtime.clear();
  tree_.muon_ECALenergy.clear();
  tree_.muon_ECALnCells.clear();
  tree_.muon_HCALtime.clear();
  tree_.muon_HCALenergy.clear();
  tree_.muon_HCALnCells.clear();
  tree_.muon_HCALTDCtime.clear();
  tree_.muon_HCALTDCenergy.clear();
  tree_.muon_HCALTDCnCells.clear();
  tree_.muon_MTDtime.clear();
  tree_.muon_MTDenergy.clear();
  tree_.muon_MTDnCells.clear();
  tree_.muon_MTDClutime.clear();
  tree_.muon_MTDCluenergy.clear();
  tree_.muon_MTDnClus.clear();
  tree_.muon_HGCALtime.clear();
  tree_.muon_HGCALenergy.clear();
  tree_.muon_HGCALnTracksters.clear();
  tree_.muon_closestGenIndex.clear();
  tree_.muon_closestGenR.clear();
  tree_.muon_closestEbGenIndex.clear();
  tree_.muon_closestEbGenR.clear();
  tree_.muon_closestEtlGenIndex.clear();
  tree_.muon_closestEtlGenR.clear();
  tree_.muon_closestHgGenIndex.clear();
  tree_.muon_closestHgGenR.clear();

  tree_.recojet_pt.clear();
  tree_.recojet_eta.clear();
  tree_.recojet_phi.clear();
  tree_.recojet_e.clear();
  tree_.recojet_ECALtime.clear();
  tree_.recojet_ECALenergy.clear();
  tree_.recojet_ECALnCells.clear();
  tree_.recojet_HCALtime.clear();
  tree_.recojet_HCALenergy.clear();
  tree_.recojet_HCALnCells.clear();
  tree_.recojet_HCALTDCtime.clear();
  tree_.recojet_HCALTDCenergy.clear();
  tree_.recojet_HCALTDCnCells.clear();
  tree_.recojet_MTDtime.clear();
  tree_.recojet_MTDenergy.clear();
  tree_.recojet_MTDnCells.clear();
  tree_.recojet_MTDClutime.clear();
  tree_.recojet_MTDCluenergy.clear();
  tree_.recojet_MTDnClus.clear();
  tree_.recojet_HGCALtime.clear();
  tree_.recojet_HGCALenergy.clear();
  tree_.recojet_HGCALnTracksters.clear();
  tree_.recojet_closestGenIndex.clear();
  tree_.recojet_closestGenR.clear();
  tree_.recojet_closestEbGenIndex.clear();
  tree_.recojet_closestEbGenR.clear();
  tree_.recojet_closestEtlGenIndex.clear();
  tree_.recojet_closestEtlGenR.clear();
  tree_.recojet_closestHgGenIndex.clear();
  tree_.recojet_closestHgGenR.clear();

  // tree_.tracksterEM_time.clear();
  // tree_.tracksterHAD_time.clear();
  // tree_.tracksterTrkEM_time.clear();
  tree_.tracksterMerge_time.clear();
  tree_.tracksterMerge_ticlIteration.clear();
  // tree_.tracksterEM_timeError.clear();
  // tree_.tracksterHAD_timeError.clear();
  // tree_.tracksterTrkEM_timeError.clear();
  tree_.tracksterMerge_timeError.clear();
  // tree_.tracksterEM_e.clear();
  // tree_.tracksterHAD_e.clear();
  // tree_.tracksterTrkEM_e.clear();
  tree_.tracksterMerge_e.clear();
  // tree_.tracksterEM_eta.clear();
  // tree_.tracksterHAD_eta.clear();
  // tree_.tracksterTrkEM_eta.clear();
  tree_.tracksterMerge_eta.clear();
  // tree_.tracksterEM_phi.clear();
  // tree_.tracksterHAD_phi.clear();
  // tree_.tracksterTrkEM_phi.clear();
  tree_.tracksterMerge_phi.clear();
  tree_.tracksterMerge_iJ.clear();
  tree_.pvtrack_pt.clear();
  tree_.pvtrack_eta.clear();
  tree_.pvtrack_phi.clear();
  tree_.pfCand_pt.clear();
  tree_.pfCand_time.clear();
  tree_.pfCand_timeError.clear();
  tree_.pfCand_eta.clear();
  tree_.pfCand_phi.clear();
  tree_.pfCand_iJ.clear();
  tree_.q_pt.clear();
  tree_.q_eta.clear();
  tree_.q_phi.clear();
  tree_.q_pdgId.clear();
  tree_.q_ctau.clear();
  tree_.q_vx.clear();
  tree_.q_vy.clear();
  tree_.q_vz.clear();
  tree_.q_ebphi.clear();
  tree_.q_ebeta.clear();
  tree_.q_ebdelay.clear();
  tree_.q_hgphi.clear();
  tree_.q_hgeta.clear();
  tree_.q_pathdelay.clear();
  tree_.q_decayInHGCAL.clear();
  tree_.q_hgdelay.clear();
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
