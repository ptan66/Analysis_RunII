#include "WZEdmAnalyzer.h"
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include "math.h"


// user include files
#include "DataFormats/Math/interface/Point3D.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


// Math
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"


// MC truth matching
//#include "RecoBTag/MCTools/interface/JetFlavour.h"
//#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include "HepMC/GenEvent.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
//#include "DataFormats/BTauReco/interface/JetTagFwd.h"


#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"


// luminosity
//#include "DataFormats/Luminosity/interface/LumiSummary.h"



// UserCode
//#include "QGLikelihoodCalculator.h"




// 
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TRandom.h"

//#include "kinematics.h"
#include "AnalysisTools.h"
#include "EdmAnalysisTools.h"
#include "lhapdfcc.h"


using namespace std;

using namespace edm;
using namespace reco;
using namespace l1extra ;
using namespace lhef;

WZEdmAnalyzer::WZEdmAnalyzer(const edm::ParameterSet& iConfig) :
  _is_debug(                    iConfig.getParameter<bool>("DEBUG")),
  _is_data(                     iConfig.getParameter<bool>("DATA")),
  _gen_only(                    iConfig.getParameter<bool>("GEN_ONLY")),
  _check_jecref(                iConfig.getParameter<bool>("CHECK_JECREF")),
  _save_allevents(              iConfig.getParameter<bool>("SAVE_ALLEVENTS")),
  _vertexing(                   iConfig.getParameter<bool>("VERTEXING")),
  _smoothing(                   iConfig.getParameter<bool>("SMOOTHING")),
  _kvfPSet(                     iConfig.getParameter<edm::ParameterSet>("KVFParameters")),
  _reco(                        iConfig.getParameter<bool>("RECO")),
  _reco_selection(              iConfig.getParameter<std::string>("RECOSELECTION")),
  BeamSpotToken_(               consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))), 
  VerticesToken_(               consumes<reco::VertexCollection> (iConfig.getParameter<std::string>("Vertices"))),
  MuonCollectionToken_(         consumes<reco::MuonCollection> (iConfig.getParameter<std::string>("Muons"))),
  ElectronCollectionToken_(     consumes<reco::GsfElectronCollection> ( iConfig.getParameter<std::string>("Electrons"))),
  _effectiveAreas(             (iConfig.getParameter<edm::FileInPath>("EffAreasConfigFile")).fullPath()),
  EleLooseIdMapToken_(          consumes<edm::ValueMap<bool> >(          iConfig.getParameter<edm::InputTag>("EleLooseIdMap"))),
  EleMediumIdMapToken_(         consumes<edm::ValueMap<bool> >(         iConfig.getParameter<edm::InputTag>("EleMediumIdMap"))),
  EleTightIdMapToken_(          consumes<edm::ValueMap<bool> >(          iConfig.getParameter<edm::InputTag>("EleTightIdMap"))), 
  ElectronEcalPFClusterIsolationProducerToken_(      consumes<edm::ValueMap<float> >(      iConfig.getParameter<edm::InputTag>("ElectronEcalPFClusterIsolationProducer"))),
  ElectronHcalPFClusterIsolationProducerToken_(      consumes<edm::ValueMap<float> >(      iConfig.getParameter<edm::InputTag>("ElectronHcalPFClusterIsolationProducer"))),
  TrigMvaValuesMapToken_(       consumes<edm::ValueMap<float> >(      iConfig.getParameter<edm::InputTag>("TrigMvaValuesMap"))),
  TrigMvaCategoriesMapToken_(   consumes<edm::ValueMap<int> >(    iConfig.getParameter<edm::InputTag>("TrigMvaCategoriesMap"))),
  TrigMvaMediumIdMapsToken_(    consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("TrigMvaMediumIdMaps"))),
  TrigMvaTightIdMapsToken_(     consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("TrigMvaTightIdMaps"))),
  NonTrigMvaValuesMapToken_(    consumes<edm::ValueMap<float> >(   iConfig.getParameter<edm::InputTag>("NonTrigMvaValuesMap"))),
  NonTrigMvaCategoriesMapToken_(consumes<edm::ValueMap<int> >( iConfig.getParameter<edm::InputTag>("NonTrigMvaCategoriesMap"))),
  NonTrigMvaMediumIdMapsToken_( consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("NonTrigMvaMediumIdMaps"))),
  NonTrigMvaTightIdMapsToken_(  consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("NonTrigMvaTightIdMaps"))),
  PFCHSJetToken_(               consumes<reco::PFJetCollection>( iConfig.getParameter<std::string>("PFCHSJets"))),
  pfchsJetFlavourInfosToken_(   consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("PFCHSJetFlavourInfos") ) ), 
  PFCHSJetTagInfos_(            iConfig.getParameter<std::vector< std::string> >("PFCHSJetTagInfos")), 
  // temporarily set to ak5 jet, 
  ak5PFJetToken_(               consumes<reco::PFJetCollection>( iConfig.getParameter<std::string>("AK5PFJets"))),
  ak5PFJetFlavourInfosToken_(   consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("AK5PFJetFlavourInfos") ) ), 
  ak5PFJetTagInfos_(            iConfig.getParameter<std::vector< std::string> >("AK5PFJetTagInfos")), 
  CaloJetToken_(                consumes<reco::CaloJetCollection>( iConfig.getParameter<std::string>("CaloJets"))),
  recoCaloToken_(               consumes<reco::JetFlavourMatchingCollection>(iConfig.getParameter<edm::InputTag>("CaloJetFlavourInfos"))), 
  CaloJetTagInfos_(             iConfig.getParameter<std::vector< std::string> >("CaloJetTagInfos")), 
  JPTJetToken_(                 consumes<reco::JPTJetCollection>( iConfig.getParameter<std::string>("JPTJets"))),
  recoJPTToken_(                consumes<reco::JetFlavourMatchingCollection>(iConfig.getParameter<edm::InputTag>("JPTJetFlavourInfos"))), 
  JPTJetTagInfos_(              iConfig.getParameter<std::vector< std::string> >("JPTJetTagInfos")), 
  PFJetToken_(                  consumes<reco::PFJetCollection>( iConfig.getParameter<std::string>("PFJets"))),
  pfJetFlavourInfosToken_(      consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("PFJetFlavourInfos") ) ), 
  PFJetTagInfos_(               iConfig.getParameter<std::vector< std::string> >("PFJetTagInfos")), 
  JetTagCollectionTags_(        iConfig.getParameter<std::vector< std::string> >("JetTagCollections")),
  jetMinPt(                     iConfig.getParameter<double>("JetMinPt")),
  leptonThreshold(              iConfig.getParameter<double>("LeptonThreshold")),
  inputJetIDValueMap(           iConfig.getParameter<edm::InputTag>("InputJetIDValueMap")), 
  pfchsJetCorrToken_(           consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("PFCHSJetCorrectionToken"))), 
  caloJetCorrToken_(            consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("CaloJetCorrectionToken"))), 
  jptJetCorrToken_(             consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("JPTJetCorrectionToken"))), 
  pfJetCorrToken_(              consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("PFJetCorrectionToken"))), 
  FixGridRhoToken_(             consumes<double> (iConfig.getParameter<edm::InputTag>("FixGridRho"))), 
  RhoSrcToken_(                 consumes<double> (iConfig.getParameter<edm::InputTag>("RhoSrc"))), 
  SigmaSrcToken_(               consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaSrc"))), 
  RhoSrcTokenCHS_(              consumes<double> (iConfig.getParameter<edm::InputTag>("RhoSrcCHS"))), 
  SigmaSrcTokenCHS_(            consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaSrcCHS"))), 
  RhoSrcTokenCalo_(             consumes<double> (iConfig.getParameter<edm::InputTag>("RhoSrcCalo"))), 
  SigmaSrcTokenCalo_(           consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaSrcCalo"))), 
  RhoIsoSrcToken_(              consumes<double> (iConfig.getParameter<edm::InputTag>("RhoIsoSrc"))), 
  SigmaIsoSrcToken_(            consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaIsoSrc"))), 
  RhoChSrcToken_(               consumes<double> (iConfig.getParameter<edm::InputTag>("RhoChSrc"))), 
  SigmaChSrcToken_(             consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaChSrc"))), 
  RhoCh2p4SrcToken_(            consumes<double> (iConfig.getParameter<edm::InputTag>("RhoCh2p4Src"))), 
  SigmaCh2p4SrcToken_(          consumes<double> (iConfig.getParameter<edm::InputTag>("SigmaCh2p4Src"))), 
  TrackCollectionToken_(        consumes<reco::TrackCollection>(iConfig.getParameter<std::string>("Tracks"))),
  trackMinPtWithMCTruth(        iConfig.getParameter<double>("TrackMinPtWithMCTruth")),
  leptonMinPtForComposition(    iConfig.getParameter<double>("LeptonMinPtForComposition")),
  PhotonCollectionToken_(       consumes<reco::PhotonCollection> (iConfig.getParameter<std::string>("PhotonCollection"))),
  L1GTReadoutRecordToken_(      consumes<L1GlobalTriggerReadoutRecord> (iConfig.getParameter<edm::InputTag>("L1GTReadoutRecordLabel"))),
  TriggerResultsToken_(         consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("TriggerResultsLabel"))),
  TriggerSummaryToken_(         consumes<trigger::TriggerEvent> (iConfig.getParameter<edm::InputTag>("TriggerSummaryLabel"))),
  TriggerRefPath_(              iConfig.getParameter<std::string>("TriggerRefPath")),
  tmpHLTTriggerMuons_(          iConfig.getParameter<std::vector< std::string> >("HLTTriggerMuons")), 
  tmpHLTTriggerMuonElectrons_(  iConfig.getParameter<std::vector< std::string> >("HLTTriggerMuonElectrons")), 
  tmpHLTTriggerElectrons_(      iConfig.getParameter<std::vector< std::string> >("HLTTriggerElectrons")), 
  GeneratorLevelToken_(         consumes<GenEventInfoProduct> (iConfig.getParameter<std::string>("GeneratorLevelTag"))),
  LHEEventProductToken_(        consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("LHEEventProductTag"))),
  GenJetAlgorithmToken_(        consumes<reco::GenJetCollection> (iConfig.getParameter<std::string>("GenJets"))),
  akGenJetAlgorithmToken_(      consumes<reco::GenJetCollection> (iConfig.getParameter<std::string>("akGenJets"))),
  genJetFlavourInfosToken_(     consumes<reco::JetFlavourInfoMatchingCollection> (iConfig.getParameter<edm::InputTag>("akGenJetFlavourInfos") ) ), 

  genJetMinPt(                  iConfig.getParameter<double>("GenJetMinPt")),
  //  SimTrackTags_(                iConfig.getParameter<std::string>("SimTracks")),
  diLeptonMinMass(              iConfig.getParameter<double>("DiLeptonMinMass")) {


  // for accessing trigger information
  //  TriggerProcess_ = TriggerResultsLabel_.process();
  TriggerProcess_ = (iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")).process();

  std::cout << "trigger process label: " << TriggerProcess_ << std::endl;

  if (!_gen_only) {
    fitter = new KalmanVertexFitter(_kvfPSet , _smoothing);
    // fitter->setMaximumDistance(vertexingMaxDistance);
    //fitter->setMaximumNumberOfIterations(vertexingMaxNumOfIterations);
  }


  lumiSummaryToken_                = consumes<LumiSummary,edm::InLumi>(edm::InputTag("lumiProducer") );



  jetID_ValueMapToken_             = consumes< edm::ValueMap<reco::JetID> >(inputJetIDValueMap);
  SuperClusterCollectionToken_     = consumes<reco::SuperClusterCollection>(edm::InputTag("superClusters"));
  ConversionCollectionToken_       = consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));
  GenParticleCollectionToken_      = consumes<reco::GenParticleCollection>(edm::InputTag( "genParticles" ));


  // btagging tokens
  jetTagsToken_                    = consumes<reco::JetTagCollection>( edm::InputTag( JetTagCollectionTags_[0] ) );
  jetTagsCSVToken_                 = consumes<reco::JetTagCollection>( edm::InputTag( JetTagCollectionTags_[1] ) );

  myPFCHSJetTagsTCHPToken_         = consumes<reco::JetTagCollection>( edm::InputTag( PFCHSJetTagInfos_[0] ) );
  myPFCHSJetTagsJPToken_           = consumes<reco::JetTagCollection>( edm::InputTag( PFCHSJetTagInfos_[1] ) );
  myPFCHSJetTagsCSVToken_          = consumes<reco::JetTagCollection>( edm::InputTag( PFCHSJetTagInfos_[2] ) );


  myAK5PFJetTagsTCHPToken_         = consumes<reco::JetTagCollection>( edm::InputTag( ak5PFJetTagInfos_[0] ) );
  myAK5PFJetTagsJPToken_           = consumes<reco::JetTagCollection>( edm::InputTag( ak5PFJetTagInfos_[1] ) );
  myAK5PFJetTagsCSVToken_          = consumes<reco::JetTagCollection>( edm::InputTag( ak5PFJetTagInfos_[2] ) );


  myCaloJetTagsTCHPToken_          = consumes<reco::JetTagCollection>( edm::InputTag( CaloJetTagInfos_[0] ) );
  myCaloJetTagsJPToken_            = consumes<reco::JetTagCollection>( edm::InputTag( CaloJetTagInfos_[1] ) );
  myCaloJetTagsCSVToken_           = consumes<reco::JetTagCollection>( edm::InputTag( CaloJetTagInfos_[2] ) );


  myJPTJetTagsTCHPToken_           = consumes<reco::JetTagCollection>( edm::InputTag( JPTJetTagInfos_[0] ) );
  myJPTJetTagsJPToken_             = consumes<reco::JetTagCollection>( edm::InputTag( JPTJetTagInfos_[1] ) );
  myJPTJetTagsCSVToken_            = consumes<reco::JetTagCollection>( edm::InputTag( JPTJetTagInfos_[2] ) );

  myPFJetTagsTCHPToken_            = consumes<reco::JetTagCollection>( edm::InputTag( PFJetTagInfos_[0] ) );
  myPFJetTagsJPToken_              = consumes<reco::JetTagCollection>( edm::InputTag( PFJetTagInfos_[1] ) );
  myPFJetTagsCSVToken_             = consumes<reco::JetTagCollection>( edm::InputTag( PFJetTagInfos_[2] ) );


  QGTagsHandleMLPToken_            = consumes<edm::ValueMap<float> > (  edm::InputTag("QGTagger","qgMLP") );
  QGTagsHandleLikelihoodToken_     = consumes<edm::ValueMap<float> > (  edm::InputTag("QGTagger","qgLikelihood") );


  // pujetID
  puJetIdFlagToken_                = consumes<edm::ValueMap<int> > (edm::InputTag("pfPileupJetId","fullId") );
  puJetIdMvaToken_                 = consumes<edm::ValueMap<float> > (edm::InputTag("pfPileupJetId","fullDiscriminant") );

  // ak5 pujetID
  puJetIdFlagAK5Token_             = consumes<edm::ValueMap<int> > (edm::InputTag("ak5PFPileupJetId","fullId") );
  puJetIdMvaAK5Token_              = consumes<edm::ValueMap<float> > (edm::InputTag("ak5PFPileupJetId","fullDiscriminant") );


  puJetIdFlagCHSToken_             = consumes<edm::ValueMap<int> > (edm::InputTag("pfchsPileupJetId","fullId") );
  puJetIdMvaCHSToken_              = consumes<edm::ValueMap<float> > (edm::InputTag("pfchsPileupJetId","fullDiscriminant") );



  ecalEBRecHitToken_               = consumes<EBRecHitCollection> (edm::InputTag("reducedEcalRecHitsEB") );
  ecalEERecHitToken_               = consumes<EERecHitCollection> (edm::InputTag("reducedEcalRecHitsEE") );
 

  tevMapH1Token_                   = consumes<reco::TrackToTrackMap>( edm::InputTag("tevMuons", "default") ); 
  tevMapH2Token_                   = consumes<reco::TrackToTrackMap>( edm::InputTag("tevMuons", "firstHit") ); 
  tevMapH3Token_                   = consumes<reco::TrackToTrackMap>( edm::InputTag("tevMuons", "picky") ); 

  PupInfoToken_                    = consumes<std::vector< PileupSummaryInfo > > ( edm::InputTag("addPileupInfo") );



  genEventInfoToken_               = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  theGenToken_                     = consumes<reco::JetFlavourMatchingCollection>(edm::InputTag("flavourByValGenJet"));
  genMETToken_                     = consumes<GenMETCollection>(edm::InputTag("genMetTrue"));


  tcMETToken_                      = consumes<METCollection>(edm::InputTag("tcMet") );
  caloMETToken_                    = consumes<CaloMETCollection>(edm::InputTag("caloMet") );
  pfMETToken_                      = consumes<PFMETCollection>(edm::InputTag("pfMet") );

  calibratedElectronsToken_        = consumes<reco::GsfElectronCollection> ( edm::InputTag("calibratedElectrons" ));






  this->displayConfig();
  totalProcessedEvts = 0;
 
 
  // initialize PDF packages for some calculation
  /*
  lhaPDFPath = getenv("LHAPATH");
  std::string pdfSet(lhaPDFPath);    pdfSet.append("/");     pdfSet.append(pdf); 
  initpdfset_((char *)pdfSet.data(), pdfSet.size());        initpdf_(subset);
  getxmin_(subset, &xmin);           getxmax_(subset, &xmax);
  getq2min_(subset, &qmin);          getq2max_(subset, &qmax);

  
  if (_is_debug) {
    std::cout  << "x range [" << xmin << ", " << xmax << "]" << std::endl;
    std::cout << "Q range [" << qmin << ", " << qmax << "]" << std::endl;
    std::cout  << "        PDF initialized" << std::endl;
    std::cout << "**************************************" << std::endl;
  }
  */

  //  if (_is_debug) std::cout  << "check point ... create output file - " << out.c_str() << std::endl;


  TClass::GetClass("_event_filterBit_")->SetCanSplit(true);
  TClass::GetClass("_gen_eventInfo_")->SetCanSplit(true);
  TClass::GetClass("_gen_ttbar_")->SetCanSplit(true);
  TClass::GetClass("_gen_DrellYan_")->SetCanSplit(true);
  TClass::GetClass("_mc_process_")->SetCanSplit(true);
  TClass::GetClass("_genwz_")->SetCanSplit(true);
  TClass::GetClass("_vec4_")->SetCanSplit(true);
  TClass::GetClass("_trg_bits_")->SetCanSplit(true);
  TClass::GetClass("_hlt_info_")->SetCanSplit(true);
  TClass::GetClass("_met_")->SetCanSplit(true);
  TClass::GetClass("_mets_")->SetCanSplit(true);
  TClass::GetClass("_dileadingjets_")->SetCanSplit(true);
  TClass::GetClass("_run_info_")->SetCanSplit(true);
  TClass::GetClass("_vertex_")->SetCanSplit(true);
  TClass::GetClass("_l1_obj_")->SetCanSplit(true);
  TClass::GetClass("_supercluster_")->SetCanSplit(true);
  TClass::GetClass("_photon_")->SetCanSplit(true);
  TClass::GetClass("_electron_")->SetCanSplit(true);
  TClass::GetClass("_beam_spot_")->SetCanSplit(true);
  TClass::GetClass("_track_")->SetCanSplit(true);
  TClass::GetClass("_muon_")->SetCanSplit(true);
  TClass::GetClass("_jet_")->SetCanSplit(true);
  TClass::GetClass("_di_jet_")->SetCanSplit(true);
  TClass::GetClass("_gen_jet_")->SetCanSplit(true);
  TClass::GetClass("_W_")->SetCanSplit(true);
  TClass::GetClass("_di_lepton_")->SetCanSplit(true);
  TClass::GetClass("_tri_lepton_")->SetCanSplit(true);
  TClass::GetClass("_quar_lepton_")->SetCanSplit(true);
  TClass::GetClass("_lepton_photon_")->SetCanSplit(true);
  TClass::GetClass("_dilepton_photon_")->SetCanSplit(true);
  TClass::GetClass("_event_")->SetCanSplit(true);

  // use TFile service?
  //  edm::Service< TFileService > fs_;
  ntuple = fs->make<TTree>("ntuple","ntuple");



  // Define output file and ntuple without TFileService
  //  hFile  = new TFile( out.c_str(), open.c_str() );
  //  ntuple = new TTree("ntuple", "");
  
  myEvent = new _event_();
  ntuple->Branch("ewk_asym", "_event_", &myEvent, 32000, 99); 
  if (_is_debug) std::cout << "check point ... finishing BeginJob() " << std::endl; 
}


WZEdmAnalyzer::~WZEdmAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


/*************************************************************************                                                             
 *                                                                                                                                         
 * main event loop                                                        
 * few levels:
 * 1) gen candidates only
 * 2)                                                                  
 *                                                                                                                                         
 **************************************************************************/                                                               

void 
WZEdmAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  if (totalProcessedEvts<=10) std::cout << "is real data ? " << iEvent.isRealData()<< std::endl;
  if (_is_debug)  std::cout << "check point ... process event " << std::endl;
  totalProcessedEvts ++;


  myEvent->Clear();  
  _is_save = false;  // default will not save the event
  

  //    if (iEvent.id().run() != 257613) return;



  // std::cout << "check point ... process event " << std::endl;

  myEvent->setEventNum(             iEvent.id().event());
  myEvent->setRunNum(               iEvent.id().run());
  myEvent->setLumiBlock(            iEvent.luminosityBlock() );
  myEvent->setBunchCrossing(        iEvent.bunchCrossing());
  myEvent->setOrbitNum(             iEvent.orbitNumber());


  //  std::cout << " average luminosity in ["
  //    << setw( 8) 
  //	    <<  iEvent.luminosityBlock() << "] = " 
  //    << avgInstLumi << std::endl;

  
  const edm::Timestamp timeStamp =  iEvent.time();
  unsigned int high = timeStamp.value() >> 32;       // seconds
  unsigned int low  = 0xffffffff & timeStamp.value(); // microseconds  
  myEvent->setTimeHigh( high );     myEvent->setTimeHigh( low );
  if (_is_debug) {
    std::cout << "Run number: "   << iEvent.id().event() << std::endl;
    std::cout << "Event number: " << iEvent.id().run() << std::endl;
  }
  


  // (++). Process only generator information
  if (_gen_only) {
 

    myMCTruth =  myEvent->getMCInfo();    
    myGenWZ   =  myEvent->getGenWZ(); 

    if (_reco) {

      iEvent.getByToken(GeneratorLevelToken_,     mcTruth);
      this->fillMCInfo(mcTruth,                 myMCTruth);
      this->fillGenWZ(mcTruth,                  myGenWZ);

    } else {



      iEvent.getByToken( GenParticleCollectionToken_,      genParticles );
      bool hasLHE = iEvent.getByToken( LHEEventProductToken_, lheEventInfo );

      if (hasLHE)  copyLHEweights( myEvent, lheEventInfo.product() );


      this->fillMCInfo(genParticles,            myMCTruth);
      this->fillGenWZ(genParticles,             myGenWZ);
      this->fillGenTTbar(genParticles,          myEvent->getGenTTbar() );
      //      this->fillGenDrellYan(genParticles,       myEvent->getGenDrellYan() );
      if (hasLHE) {
	this->fillGenDrellYan(genParticles,       lheEventInfo.product() , myEvent->getGenDrellYan() );
      } else {

	this->fillGenDrellYan(genParticles,       0 , myEvent->getGenDrellYan() );

      }
    }

    _is_save = true; // if generator level, all events will be save.
    ntuple->Fill();
    return;
  }




  // -- event filters;
  //edm::Handle<bool> trackingFailure;
  // iEvent.getByToken("trackingFailureFilter", trackingFailure);
  //  iEvent.getByLabel("trackingFailureFilter", trackingFailure);
  //myEvent->getEventFilterBit()->trackingFailureFilter =  *(trackingFailure);



  /*
  Handle<bool> EcalDeadCellTriggerPrimitiveFilterHandle;
  iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter", EcalDeadCellTriggerPrimitiveFilterHandle);
  if (EcalDeadCellTriggerPrimitiveFilterHandle.isValid()) myEvent->getEventFilterBit()->EcalDeadCellTriggerPrimitiveFilter  = (Bool_t)(*EcalDeadCellTriggerPrimitiveFilterHandle);
  

  Handle<bool> ecalLaserCorrFilterHandle;
  iEvent.getByLabel("ecalLaserCorrFilter", ecalLaserCorrFilterHandle);
  if (ecalLaserCorrFilterHandle.isValid()) myEvent->getEventFilterBit()->ecalLaserCorrFilter  = (Bool_t)(*ecalLaserCorrFilterHandle);
 

  Handle<bool> eeBadScFilterHandle;
  iEvent.getByLabel("eeBadScFilter", eeBadScFilterHandle);
  if (eeBadScFilterHandle.isValid()) myEvent->getEventFilterBit()->eeBadScFilter  = (Bool_t)(*eeBadScFilterHandle);
  
 Handle<bool> hcalLaserEventFilterHandle;
  iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventFilterHandle);
  if (hcalLaserEventFilterHandle.isValid()) myEvent->getEventFilterBit()->hcalLaserEventFilter  = (Bool_t)(*hcalLaserEventFilterHandle);
 

Handle<bool> CSCTightHaloFilterHandle;
  iEvent.getByLabel("CSCTightHaloFilter", CSCTightHaloFilterHandle);
  if (CSCTightHaloFilterHandle.isValid()) myEvent->getEventFilterBit()->CSCTightHaloFilter  = (Bool_t)(*CSCTightHaloFilterHandle);
  
*/








  // -- Magnetic field
  ESHandle<MagneticField> MF;
  iSetup.get<IdealMagneticFieldRecord>().get(MF);
  const MagneticField* theMagneticField = MF.product();
  myEvent->setBField( fabs(theMagneticField->inTesla(GlobalPoint(0,0,0)).z()));
  bField =  fabs(theMagneticField->inTesla(GlobalPoint(0,0,0)).z());


  iEvent.getByToken(BeamSpotToken_, recoBeamSpotHandle);              
  //  iEvent.getByLabel(BeamSpotTags_, recoBeamSpotHandle);  
            
  if (recoBeamSpotHandle.isValid() ) {
    vertexBeamSpot = recoBeamSpotHandle.product();
  } else {
    vertexBeamSpot = 0;
  } 


  /*********************************************************************
   *
   *
   *
   *
   * (++). Acess the jet ID value map
   *
   *********************************************************************/
  iEvent.getByToken(jetID_ValueMapToken_,jetID_ValueMap_Handle);


  // (++). Access the jet correction; using JetCorrector since 76x release
  iEvent.getByToken(pfchsJetCorrToken_, pfchsJetCorr);
  //  iEvent.getByToken(caloJetCorrToken_,  caloJetCorr);
  // iEvent.getByToken(jptJetCorrToken_,   jptJetCorr);
  iEvent.getByToken(pfJetCorrToken_,    pfJetCorr);



  // (++) Access the uncertainty of the jet energy correction
  // pf jet 
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",pfchsJetCorParColl);
  if ( pfchsJetCorParColl.isValid() ) {
    pfchsJetUnc = new JetCorrectionUncertainty( (*pfchsJetCorParColl)["Uncertainty"]  );
 
  } else {

    pfchsJetUnc = 0;
  }

  /*
  iSetup.get<JetCorrectionsRecord>().get("AK5Calo",caloJetCorParColl);
  if (caloJetCorParColl.isValid()) {
    caloJetUnc = new JetCorrectionUncertainty( (*caloJetCorParColl)["Uncertainty"]  );

  } else {
    caloJetUnc = 0;
  }
  jptJetUnc  = 0;
  */

  // pf jet 
  iSetup.get<JetCorrectionsRecord>().get("AK4PF",pfJetCorParColl);
  if (pfJetCorParColl.isValid() ) {
    pfJetUnc   = new JetCorrectionUncertainty( (*pfJetCorParColl)["Uncertainty"]  );
  } else {

    pfJetUnc = 0;
  }



  // get the rho values from fast jet 
  //edm::Handle<double> rhoHandle;
  // default rho value for ak4PFjet
  // NOTE74

  iEvent.getByToken(FixGridRhoToken_, fixGridRhoHandle);  

  iEvent.getByToken(RhoSrcToken_,          rhoHandle);         
  iEvent.getByToken(SigmaSrcToken_,        sigmaHandle);         
  if (rhoHandle.isValid() && sigmaHandle.isValid() ) {
    myEvent->rho      = *rhoHandle;
    myEvent->sigma    = *sigmaHandle;
  }

  // default rho value for ak4PFjetCHS
  iEvent.getByToken(RhoSrcTokenCHS_,          rhoHandleCHS);         
  iEvent.getByToken(SigmaSrcTokenCHS_,        sigmaHandleCHS);  
  if (rhoHandleCHS.isValid() && sigmaHandleCHS.isValid() ) {

    myEvent->rhoCHS      = *rhoHandleCHS;
    myEvent->sigmaCHS    = *sigmaHandleCHS;
  }

  
  // default rho value for ak4CaloJet
  iEvent.getByToken(RhoSrcTokenCalo_,          rhoHandleCalo);         
  iEvent.getByToken(SigmaSrcTokenCalo_,        sigmaHandleCalo);         
  if (rhoHandleCalo.isValid() && sigmaHandleCalo.isValid() ) {

    myEvent->rhoCalo      = *rhoHandleCalo;
    myEvent->sigmaCalo    = *sigmaHandleCalo;
  }

  // default rho value for track jet
  //  iEvent.getByToken(RhoSrcTokenTrack_,          rhoHandleTrack);         
  // iEvent.getByToken(SigmaSrcTokenTrack_,        sigmaHandleTrack);         
  // if (rhoHandleTrack.isValid() && sigmaHandleTrack.isValid() ) {
  //  myEvent->rhoTrack      = *rhoHandleTrack;
  //  myEvent->sigmaTrack    = *sigmaHandleTrack;
  // }

  // others
  iEvent.getByToken(RhoIsoSrcToken_,       rhoIsoHandle);   
  iEvent.getByToken(SigmaIsoSrcToken_,     sigmaIsoHandle);   
  myEvent->rhoIso   = *rhoIsoHandle;
  myEvent->sigmaIso = *sigmaIsoHandle;



  iEvent.getByToken(RhoChSrcToken_,       rhoChHandle);   
  iEvent.getByToken(SigmaChSrcToken_,     sigmaChHandle);   
  myEvent->rhoCh   = *rhoChHandle;
  myEvent->sigmaCh = *sigmaChHandle;



  // eta up to 2.4
  iEvent.getByToken(RhoCh2p4SrcToken_,       rhoCh2p4Handle);   
  iEvent.getByToken(SigmaCh2p4SrcToken_,     sigmaCh2p4Handle);   
  myEvent->rhoCh2p4   = *rhoCh2p4Handle;
  myEvent->sigmaCh2p4 = *sigmaCh2p4Handle;



  // (++). Setup for vertex
  // Be careful, maybe do not work well on CMSSW_3XX
  //  if (_vertexing) 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transientTrackBuilder);




  // (++). Extract trigger information
  iEvent.getByToken(L1GTReadoutRecordToken_,          l1TriggerReadoutRec);
  iSetup.get<L1GtTriggerMenuRcd>().get(               l1GTTriggerMenu);


  iEvent.getByToken(TriggerResultsToken_,             triggerResults);
  iEvent.getByToken(TriggerSummaryToken_,             triggerObj);

  // (++). Define and extract container item
  iEvent.getByToken(VerticesToken_,                   recVtxs);
  iEvent.getByToken(TrackCollectionToken_,            tracks);
  iEvent.getByToken(MuonCollectionToken_,              muons);
  iEvent.getByToken(ElectronCollectionToken_,          electrons);
  iEvent.getByToken(PFCHSJetToken_,                    pfchsJets);
  iEvent.getByToken(ak5PFJetToken_,                    ak5pfJets);

  // iEvent.getByToken(CaloJetToken_,                     caloJets);
  // iEvent.getByToken(JPTJetToken_,                      jptJets);
  iEvent.getByToken(PFJetToken_,                       pfJets);
  iEvent.getByToken(PhotonCollectionToken_,            photons);
  iEvent.getByToken(SuperClusterCollectionToken_,      superclusters);
  iEvent.getByToken(ConversionCollectionToken_,        convCol);


  iEvent.getByToken(calibratedElectronsToken_,         calibratedElectrons);

  /*********************************************************************
   *
   *
   *
   *
   * (++). Acess electron ID values maps
   *       considering cut-based and MVA electron ID
   *
   *********************************************************************/
  iEvent.getByToken(EleLooseIdMapToken_ ,             loose_id_decisions);
  iEvent.getByToken(EleMediumIdMapToken_,             medium_id_decisions);
  iEvent.getByToken(EleTightIdMapToken_,              tight_id_decisions);

  iEvent.getByToken(ElectronEcalPFClusterIsolationProducerToken_, electronEcalPFClusterIsolation);
  iEvent.getByToken(ElectronHcalPFClusterIsolationProducerToken_, electronHcalPFClusterIsolation);


  iEvent.getByToken(TrigMvaValuesMapToken_,           trigMvaValues);
  iEvent.getByToken(TrigMvaCategoriesMapToken_,       trigMvaCategories);
  iEvent.getByToken(TrigMvaMediumIdMapsToken_,        trigMvaMedium_id_decisions);
  iEvent.getByToken(TrigMvaTightIdMapsToken_,         trigMvaTight_id_decisions);


  iEvent.getByToken(NonTrigMvaValuesMapToken_,        nonTrigMvaValues);
  iEvent.getByToken(NonTrigMvaCategoriesMapToken_,    nonTrigMvaCategories);
  iEvent.getByToken(NonTrigMvaMediumIdMapsToken_,     nonTrigMvaMedium_id_decisions);
  iEvent.getByToken(NonTrigMvaTightIdMapsToken_,      nonTrigMvaTight_id_decisions);


  // regression electron calibration


  // btagging 
  iEvent.getByToken(jetTagsToken_,             jetTags);
  iEvent.getByToken(jetTagsCSVToken_,          jetTagsCSV);


  // ak4pfchs
  iEvent.getByToken(myPFCHSJetTagsTCHPToken_,  myPFCHSJetTagsTCHP);
  iEvent.getByToken(myPFCHSJetTagsJPToken_,    myPFCHSJetTagsJP);
  iEvent.getByToken(myPFCHSJetTagsCSVToken_,   myPFCHSJetTagsCSV);



  iEvent.getByToken(myAK5PFJetTagsTCHPToken_,     myAK5PFJetTagsTCHP);
  iEvent.getByToken(myAK5PFJetTagsJPToken_,       myAK5PFJetTagsJP);
  iEvent.getByToken(myAK5PFJetTagsCSVToken_,      myAK5PFJetTagsCSV);



  /*
  iEvent.getByToken(myCaloJetTagsTCHPToken_,   myCaloJetTagsTCHP);
  iEvent.getByToken(myCaloJetTagsJPToken_,     myCaloJetTagsJP);
  iEvent.getByToken(myCaloJetTagsCSVToken_,    myCaloJetTagsCSV);


  iEvent.getByToken(myJPTJetTagsTCHPToken_,    myJPTJetTagsTCHP);
  iEvent.getByToken(myJPTJetTagsJPToken_,      myJPTJetTagsJP);
  iEvent.getByToken(myJPTJetTagsCSVToken_,     myJPTJetTagsCSV);
  */


  iEvent.getByToken(myPFJetTagsTCHPToken_,     myPFJetTagsTCHP);
  iEvent.getByToken(myPFJetTagsJPToken_,       myPFJetTagsJP);
  iEvent.getByToken(myPFJetTagsCSVToken_,      myPFJetTagsCSV);


  // quark gluon likelihood separator
  iEvent.getByToken(QGTagsHandleMLPToken_,        QGTagsHandleMLP);
  iEvent.getByToken(QGTagsHandleLikelihoodToken_, QGTagsHandleLikelihood);



  // access corrected METs
  //  iEvent.getByLabel("met",                            rawMEThandle);
  //  iEvent.getByLabel("corMetGlobalMuons",              muCorrMEThandle);
  //  iEvent.getByLabel("tcMet",                          tcMEThandle);
  //  iEvent.getByLabel("pfMet",                          pfMEThandle);
  //  iEvent.getByLabel("metMuonJESCorAK5",               muJESCorrMEThandle);

  iEvent.getByToken(ecalEBRecHitToken_,           ecalEBRecHitHandle);
  iEvent.getByToken(ecalEERecHitToken_,           ecalEERecHitHandle);
 

  recoPhotons = 0;
  recoTracks  = 0;
  recoMuons   = 0;
  recoJets    = 0;
  caloTowers  = 0;

  recoTracks       =   tracks.product();
  recoMuons        =   muons.product();
  recoPhotons      =   photons.product();
  recoElectrons    =   electrons.product(); 
  recoJetTags      =   jetTags.product();


  // TeV Muon Refit
  iEvent.getByToken(tevMapH1Token_, tevMapH1);   tevMap1 = *(tevMapH1.product());
  iEvent.getByToken(tevMapH2Token_, tevMapH2);   tevMap2 = *(tevMapH2.product());
  iEvent.getByToken(tevMapH3Token_, tevMapH3);   tevMap3 = *(tevMapH3.product());
  

  // filling the leading primary vertex
  if (recVtxs->size() > 0) {

    reco::VertexRef vtx(recVtxs, 0);    

    myEvent->vtxPosX      = vtx->x();
    myEvent->vtxPosY      = vtx->y();
    myEvent->vtxPosZ      = vtx->z();
    myEvent->vtxPosXError = vtx->xError();
    myEvent->vtxPosYError = vtx->yError();
    myEvent->vtxPosZError = vtx->zError();
  }


  //  return;

  /************************************************************************
   * (++). Filter some events. Add the filter requirement Here
   * filter some events.
   * for signal MC, every event will be processed; otherwise,
   * only events with at least one muon candidate will be saved
   * to save space  
   * !!!!!!!!!!!!now we are saving every MC event processed !!!!!!!!!!!!!!!!!
   ************************************************************************/       
  _is_save = true;

   
  // (++). Process Monte Carlo events  
  if (!_is_data)   { 

    // pile up reweighting

    if ( iEvent.getByToken(PupInfoToken_, PupInfo) ) {
      

      std::vector<PileupSummaryInfo>::const_iterator PVInfo;
      for (PVInfo = PupInfo->begin(); PVInfo != PupInfo->end(); ++PVInfo) {
	
	int BX = PVInfo->getBunchCrossing();
	
	if(BX == -1) {
	  myEvent->mcVertexNumm1    = PVInfo->getPU_NumInteractions();
	}
	if(BX == 0) {
	  // this is for 2011 mc 3D reweighting, etc. 
	  myEvent->mcVertexNum      = PVInfo->getPU_NumInteractions();
	  myEvent->mcVertexNumTruth = PVInfo->getTrueNumInteractions();
	}
	if(BX == 1) {
	  myEvent->mcVertexNump1    = PVInfo->getPU_NumInteractions();
	}
	
      }
    }
    
    myMCTruth        =  myEvent->getMCInfo();             
    myGenWZ          =  myEvent->getGenWZ();  
    
    hasGenJets = iEvent.getByToken( GenJetAlgorithmToken_,    genJets );
    iEvent.getByToken( akGenJetAlgorithmToken_,               akGenJets );    
    // iEvent.getByLabel( "genParticleCandidates",           genParticles );
    // iEvent.getByLabel(GeneratorLevelTag_,                 mcTruth);
    // iEvent.getByLabel(SimTrackTags_,                      simTracks);

    iEvent.getByToken(theGenToken_ ,                         theGenTag);
    iEvent.getByToken(genEventInfoToken_,                    genEventInfo);


    /*************************************************************************
     *
     * jet flavor, for 76x
     *
     *************************************************************************/
    //    iEvent.getByToken( recoCaloToken_,              theRecoCaloTag);
    //    iEvent.getByToken( recoJPTToken_,               theRecoJPTTag);
    iEvent.getByToken( pfchsJetFlavourInfosToken_,  thePFCHSJetFlavourInfos );
    iEvent.getByToken( ak5PFJetFlavourInfosToken_,  theAK5PFJetFlavourInfos);
    iEvent.getByToken( pfJetFlavourInfosToken_,     thePFJetFlavourInfos);
    iEvent.getByToken( genJetFlavourInfosToken_,    theGenJetFlavourInfos);



    myEvent->setEventWeight(genEventInfo->weight());
    iEvent.getByToken(genMETToken_,                 genMEThandle);

    myEvent->getMETs()->genMET.pt                  = (genMEThandle->front() ).et();
    myEvent->getMETs()->genMET.phi                 = (genMEThandle->front() ).phi();
    myEvent->getMETs()->genMET.sumEt               = (genMEThandle->front() ).sumEt();
    myEvent->getMETs()->genMET.e_longitudinal      = (genMEThandle->front() ).e_longitudinal();
    myEvent->getMETs()->genMET.MuonEtFraction      = (genMEThandle->front() ).MuonEtFraction();
    myEvent->getMETs()->genMET.InvisibleEtFraction = (genMEThandle->front() ).InvisibleEtFraction();


    // Feb 07, 2013: only the aod version is being kept, in general _reco is false. 
    if (_reco) {

      iEvent.getByToken(GeneratorLevelToken_,           mcTruth);
      this->fillMCInfo(mcTruth,                       myMCTruth);
      this->fillGenWZ(mcTruth,                        myGenWZ);
 

    } else {

      
      iEvent.getByToken( GenParticleCollectionToken_,              genParticles );
      bool hasLHE = iEvent.getByToken( LHEEventProductToken_, lheEventInfo );

       
      if (hasLHE)  copyLHEweights( myEvent, lheEventInfo.product() );

      this->fillMCInfo(genParticles,                  myMCTruth);
      this->fillGenWZ(genParticles,                   myGenWZ);
      this->fillGenEventInfo(genEventInfo,            myEvent->getGenEventInfo() );
      this->fillGenTTbar(genParticles,                myEvent->getGenTTbar() );
      

      if (hasLHE) {
	this->fillGenDrellYan(genParticles,       lheEventInfo.product() , myEvent->getGenDrellYan() );
      } else {
	
	this->fillGenDrellYan(genParticles,       0 , myEvent->getGenDrellYan() );

      }


    }

    this->fillGenJets();
    this->fillAKGenJets();
  }


  // (++). Out put cross section information
  this->fillRunInfo();  
  if (_is_debug) std::cout << "checking point .. finish filling run information" << std::endl;

  

  if (_is_debug) std::cout << "checking point .. finish candidate composition" << std::endl;


  // (++). Copy some of the event information for both real data and MC
  this->fillEventInfo(iEvent, iSetup);



  if (_is_save)  ntuple->Fill();


  if (pfchsJetUnc) delete pfchsJetUnc;
  // if (jptJetUnc)   delete jptJetUnc;
  // if (caloJetUnc)  delete caloJetUnc;
  if (pfJetUnc)    delete pfJetUnc;
  return;
}


void 
WZEdmAnalyzer::beginJob() 
{ 
}

void
WZEdmAnalyzer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{

  bool changed(true);

  // initialize the HLT trigger config. 

  std::cout <<"Trigger Results: " << TriggerProcess_ <<std::endl;

  if (hltConfig.init(iRun,iSetup,  TriggerProcess_ , changed)) {
    if (changed) {



      //      refIndex = hltConfig.triggerIndex(TriggerRefPath_);
      // if (refIndex >= hltConfig.size()) {


      //	refIndex = 0;

      //	edm::LogWarning("HLTrigReport") << "requested reference path '"+TriggerRefPath_+"' not in HLT menu. " << "using HLTriggerFinalPath instead." ;
 
      //	TriggerRefPath_ = "HLTriggerFinalPath";

      //	refIndex = hltConfig.triggerIndex(TriggerRefPath_);
      //	if (refIndex>= hltConfig.size()) {


      //	  refIndex = 0;
      //	  edm::LogWarning("HLTrigReport") << "requested reference path '"+TriggerRefPath_+"' not in HLT menu. " << "using first path instead." ;

      //	}
      //  }
    }
    //    hltConfig.dump("Triggers");

    //    hltConfig.dump("Streams");
    //  hltConfig.dump("Datasets");
    // hltConfig.dump("PrescaleTable");
    // hltConfig.dump("ProcessPSet");

    triggerNames  = hltConfig.triggerNames();


    checkHLTVersion( triggerNames, tmpHLTTriggerMuons_,         HLTTriggerMuons_); 
    checkHLTVersion( triggerNames, tmpHLTTriggerElectrons_,     HLTTriggerElectrons_); 
    checkHLTVersion( triggerNames, tmpHLTTriggerMuonElectrons_, HLTTriggerMuonElectrons_); 

    //    for (std::vector< std::string >::const_iterator it = triggerNames.begin(); it != triggerNames.end(); it ++) {
    //
    //      std::cout << *it << std::endl;
    // }

  } else {
    cout << "HLTEventAnalyzerAOD::analyze:"
	 << " config extraction failure with process name "
	 << "HLT" << endl;
  }
}

void 
WZEdmAnalyzer::endJob() 
{

  //  write results into disk
  std::cout << "MY INFORMATION: total processed events " << totalProcessedEvts << std::endl;
  
  return;
}



// ------------ method called when starting to processes a luminosity block  ------------
void 
WZEdmAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&iLumiBlock, edm::EventSetup const&iSetup)
{


  // access lumi-block averaged 
  avgInstLumi = 0;
  //  myavginstlumi=0;
  if (_is_data) {

    //    edm::EDGetTokenT<LumiSummary> lumiSummaryToken_;
    //    lumiSummaryToken_ = iC.consumes<LumiSummary,edm::InLumi>(lumiInputTag_);

    edm::Handle<LumiSummary>                 lumiSummaryHandle;
    iLumiBlock.getByToken(lumiSummaryToken_, lumiSummaryHandle);


    //  edm::EDGetTokenT<GenEventInfoProduct>     GeneratorLevelToken_;
    //    iLumiBlock.getByLabel("lumiProducer", lumiSummaryHandle); 
    //    iLumiBlock.getByLabel(LumiSummaryToken_, lumiSummaryHandle); 
    // Check that there is something
    if (lumiSummaryHandle.isValid())  avgInstLumi=lumiSummaryHandle->avgInsDelLumi();
  }

}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
WZEdmAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&iLumiBlock, edm::EventSetup const&iSetup)
{
}

/**************************************************************************                                                                
 *                                                                                                                                         
 *                                                                                                                                         
 *                                                                                                                                         
 **************************************************************************/                                                               
void 
WZEdmAnalyzer::displayConfig(void)
{

  std::cout << std::endl;
  std::cout <<setw(12) <<""<< "   configuration settings    "         << std::endl;
  std::cout <<setw(12) <<"*"<< "******************************"         << "*" << std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Debug flag ------ "<< setw(5)<< _is_debug    <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Data  flag ------ "<< setw(5)<< _is_data     <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Save all events - "  << setw(5)<< _save_allevents   <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Gen event flag -- "  << setw(5)<< _gen_only    <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Vertexing flag -- "  << setw(5)<< _vertexing   <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<<setw(15)<< left<< "  Reco flag- " << right<< setw(15)<< _reco_selection.data() <<"*"<< std::endl;
  std::cout <<right<< setw(12) <<"*"<< "******************************"         <<"*"<< std::endl;

  return;
}


void 
WZEdmAnalyzer::checkHLTVersion(std::vector< std::string > & triggerNames, std::vector<std::string> &atmpset,  std::vector<std::string> &aset) {

  aset.resize(0);

  std::string ahlt, aref;
  for (std::vector< std::string>::const_iterator it = atmpset.begin(); it != atmpset.end(); ++ it) {

    ahlt.resize(0);
    ahlt += *it;

    //    std::cout << ahlt<< std::endl;
    if (ahlt.size() ==0) {
      aset.push_back( ahlt ); continue;
    }

    if (ahlt.find("VERSION") != string::npos )  ahlt.replace(ahlt.find("VERSION"), std::string("VERSION").size(), ""); // we only use * to represent the version number

    //    std::cout << ahlt<< std::endl;


    for (std::vector< std::string >::const_iterator it2 = triggerNames.begin(); it2!= triggerNames.end(); ++ it2) {
      aref.resize(0);
      aref += (*it2);

      if (aref.find( ahlt) == string::npos) continue;
      //   if (aref.compare( ahlt) >=0) {
      ahlt.resize(0); ahlt += aref;
      break;
	// }
    }

    //    std::cout << ahlt<< std::endl;
    aset.push_back( ahlt );

  }

  for (std::vector< std::string>::const_iterator it = aset.begin(); it != aset.end(); it ++) {

    std::cout << *it << std::endl;
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(WZEdmAnalyzer);
