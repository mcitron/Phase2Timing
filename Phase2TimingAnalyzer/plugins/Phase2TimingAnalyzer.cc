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

//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TMath.h"
#include "TTree.h"


//
// class declaration
//
struct tree_struc_{
    std::vector<float> q_eta;
    std::vector<float> q_phi;
    std::vector<float> q_pt;
    std::vector<float> q_vx;
    std::vector<float> q_vy;
    std::vector<float> q_vz;
    std::vector<float> q_ebeta;
    std::vector<float> q_ebphi;
    std::vector<float> q_ebdelay;
    std::vector<float> q_hgeta;
    std::vector<float> q_hgphi;
    std::vector<float> q_hgdelay;
  int                           nrecojets;
  int                           ngen;
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
  std::vector<float>  recojet_closestHgGenIndex;
  std::vector<float>  recojet_closestHgGenR;
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
  _genParticles(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  _genParticlesH(),
  _recoak4PFJets(consumes< edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("recoak4PFJets"))),
  _recoak4PFJetsH(),
  ecalRecHitsEBToken_{consumes<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>>(												       iConfig.getParameter<edm::InputTag>("ebRecHitsColl"))},
  _ecalRecHitsEBH(),
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
  Handle<View<ticl::Trackster>> tracksterEMH;
  Handle<View<ticl::Trackster>> tracksterMergeH;
  Handle<View<ticl::Trackster>> tracksterHADH;
  Handle<View<ticl::Trackster>> tracksterTrkEMH;
  Handle<View<ticl::Trackster>> tracksterTrkH;

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
  int ngen = 0;
  std::vector<float> q_eta;
  std::vector<float> q_phi;
  std::vector<float> q_pt;
  std::vector<float> q_vx;
  std::vector<float> q_vy;
  std::vector<float> q_vz;
  std::vector<float> q_ebeta;
  std::vector<float> q_ebphi;
  std::vector<float> q_ebdelay;
  std::vector<float> q_hgeta;
  std::vector<float> q_hgphi;
  std::vector<float> q_hgdelay;
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
  std::vector<float>  recojet_closestHgGenIndex;
  std::vector<float>  recojet_closestHgGenR;

  
  bool debug=0;


  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON GENPARTICLES --------------------------------------"<<std::endl; 
  for (const auto & genpar_iter : *_genParticlesH){

    if (genpar_iter.mother(0) == NULL)continue;
    if(abs(genpar_iter.pdgId()) !=5  || genpar_iter.status() != 23)continue;
    float vx = genpar_iter.vertex().x();
    float vy = genpar_iter.vertex().y();
    float vz = genpar_iter.vertex().z();
    interestingquark++;
    reco::GenParticle * genParticleMother = (reco::GenParticle *) genpar_iter.mother();
    std::vector<double> ecalIntersection = _jetTimingTools.surfaceIntersection(genpar_iter,*genParticleMother,130);
    std::vector<double> hgcalIntersection = _jetTimingTools.endCapIntersection(genpar_iter,*genParticleMother,300,520);
    ngen++;
    q_pt.push_back(genpar_iter.pt());
    q_eta.push_back(genpar_iter.eta());
    q_phi.push_back(genpar_iter.phi());
    q_vx.push_back(vx);
    q_vy.push_back(vy);
    q_vz.push_back(vz);
    q_ebphi.push_back(ecalIntersection[1]);
    q_ebeta.push_back(ecalIntersection[0]);
    q_ebdelay.push_back(ecalIntersection[3]);
    q_hgphi.push_back(hgcalIntersection[1]);
    q_hgeta.push_back(hgcalIntersection[0]);
    q_hgdelay.push_back(hgcalIntersection[3]);
  }

  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);
  auto const& mtdRecHitsBTL = iEvent.get(mtdRecHitsBTLToken_);
  auto const& mtdRecHitsETL = iEvent.get(mtdRecHitsETLToken_);
  //  auto const& mtdClusBTL = iEvent.get(btlRecCluToken_);
  //  auto const& mtdClusETL = iEvent.get(btlRecCluToken_);

  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON RECO JETS --------------------------------------"<<std::endl; 
  for (const auto & recojet_iter : *_recoak4PFJetsH){

    if(recojet_iter.pt()<20)continue;
    if(fabs(recojet_iter.eta())>3)continue;
    
    nrecojets++;
    recojet_pt.push_back(recojet_iter.pt());
    recojet_eta.push_back(recojet_iter.eta());
    recojet_phi.push_back(recojet_iter.phi());
    recojet_e.push_back(recojet_iter.energy());
    int closestGenIndex = -1;
    float closestGenR = 999;
    int closestEbGenIndex = -1;
    float closestEbGenR = 999;
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
	
	float ebDeltaR = jetVec.DeltaR(ebPos);
	if (ebDeltaR < 0.4 and closestEbGenR > ebDeltaR){
	    closestEbGenR = ebDeltaR;
	    closestEbGenIndex = ig;
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
    std::vector<ticl::Trackster> tracksters;
    // tracksters.reserve(tracksters.size() + distance(tracksterEMH->begin(),tracksterEMH->end()));
    // tracksters.insert(tracksters.end(),tracksterEMH->begin(),tracksterEMH->end());
    tracksters.reserve(tracksters.size() + distance(tracksterMergeH->begin(),tracksterMergeH->end()));
    tracksters.insert(tracksters.end(),tracksterMergeH->begin(),tracksterMergeH->end());
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
  for (int ig = 0; ig < ngen; ig++){
    tree_.q_pt.push_back(q_pt[ig]);
    tree_.q_eta.push_back(q_eta[ig]);
    tree_.q_phi.push_back(q_phi[ig]); 
    tree_.q_vx.push_back(q_vx[ig]);
    tree_.q_vy.push_back(q_vy[ig]);
    tree_.q_vz.push_back(q_vz[ig]); 
    tree_.q_ebphi.push_back(q_ebphi[ig]); 
    tree_.q_ebeta.push_back(q_ebeta[ig]); 
    tree_.q_ebdelay.push_back(q_ebdelay[ig]); 
    tree_.q_hgphi.push_back(q_hgphi[ig]); 
    tree_.q_hgeta.push_back(q_hgeta[ig]); 
    tree_.q_hgdelay.push_back(q_hgdelay[ig]); 
  }
  tree_.nrecojets        = nrecojets;
  for (int ij = 0; ij < nrecojets; ij++){
    tree_.recojet_closestGenIndex.push_back(recojet_closestGenIndex[ij]);
    tree_.recojet_closestGenR.push_back(recojet_closestGenR[ij]);
    tree_.recojet_closestEbGenIndex.push_back(recojet_closestEbGenIndex[ij]);
    tree_.recojet_closestEbGenR.push_back(recojet_closestEbGenR[ij]);
    tree_.recojet_closestHgGenIndex.push_back(recojet_closestHgGenIndex[ij]);
    tree_.recojet_closestHgGenR.push_back(recojet_closestHgGenR[ij]);
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
    tree_.recojet_MTDClutime.push_back(recojet_MTDClutime[ij]);
    tree_.recojet_MTDCluenergy.push_back(recojet_MTDCluenergy[ij]);
    tree_.recojet_MTDnClus.push_back(recojet_MTDnClus[ij]);
    tree_.recojet_HGCALtime.push_back(recojet_HGCALtime[ij]);
    tree_.recojet_HGCALenergy.push_back(recojet_HGCALenergy[ij]);
    tree_.recojet_HGCALnTracksters.push_back(recojet_HGCALnTracksters[ij]);
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
  tree->Branch("ngen",              &tree_.ngen,                "ngen/I");
  tree->Branch("q_eta", &tree_.q_eta);
  tree->Branch("q_phi", &tree_.q_phi);
  tree->Branch("q_pt", &tree_.q_pt);
  tree->Branch("q_vx", &tree_.q_vx);
  tree->Branch("q_ebeta", &tree_.q_ebeta);
  tree->Branch("q_ebphi", &tree_.q_ebphi);
  tree->Branch("q_ebdelay", &tree_.q_ebdelay);
  tree->Branch("q_hgeta", &tree_.q_hgeta);
  tree->Branch("q_hgphi", &tree_.q_hgphi);
  tree->Branch("q_hgdelay", &tree_.q_hgdelay);
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
  tree->Branch("recoJet_closestHgGenIndex",&tree_.recojet_closestHgGenIndex);
  tree->Branch("recoJet_closestHgGenR",&tree_.recojet_closestHgGenR);
}


// ------------ initialize trees ------------                                                                                                                 
void Phase2TimingAnalyzer::initTreeStructure()
{
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
  tree_.recojet_closestHgGenIndex.clear();
  tree_.recojet_closestHgGenR.clear();
  tree_.q_pt.clear();
  tree_.q_eta.clear();
  tree_.q_phi.clear();
  tree_.q_vx.clear();
  tree_.q_vy.clear();
  tree_.q_vz.clear();
  tree_.q_ebphi.clear();
  tree_.q_ebeta.clear();
  tree_.q_ebdelay.clear();
  tree_.q_hgphi.clear();
  tree_.q_hgeta.clear();
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
