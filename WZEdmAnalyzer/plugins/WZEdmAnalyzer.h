#ifndef _WZEdmAnalyzer_
#define _WZEdmAnalyzer_


#include <memory>
#include <string>
#include <iostream>
#include "math.h"
#include <map>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include "HepMC/GenEvent.h"
#include "HepMC/SimpleVector.h"




// L1 triggers
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"


// jets
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
 

// tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

// muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/Vector3D.h"

// common between electrons and photons
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


// photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"


// electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"


// btagging
#include "DataFormats/BTauReco/interface/JetTag.h"


#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"



#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"


#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"


// vertex fitter
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/CachingVertex.h"


// vertex stuff
#include <DataFormats/VertexReco/interface/VertexFwd.h>



// MC jet flavor truth matching
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"


#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>



#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// for jet met corrections
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"


// jet correction
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"



// TeV Muon Reconstruction
//#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"




// luminosity
#include "DataFormats/Luminosity/interface/LumiSummary.h"




// UserCode
//#include "QGLikelihoodCalculator.h"


// ROOT libraries and user defined
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TRandom.h"

#include "kinematics.h"
#include "AnalysisTools.h"
#include "lhapdfcc.h"

using namespace edm;
using namespace reco;
using namespace std;


  


class WZEdmAnalyzer : public edm::EDAnalyzer {



 protected:

  void   setDebugPoint(  void) {  _is_debug = true;}
  void   clearDebugPoint(void){ _is_debug = false;}

 public:



  //  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  // typedef std::vector< edm::Handle< edm::ValueMap<double> > >           IsoDepositVals;

  
  explicit WZEdmAnalyzer(const edm::ParameterSet&);
  ~WZEdmAnalyzer();


 private:


  virtual void beginJob() ;
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginLuminosityBlock(edm::LuminosityBlock const&iLumiBlock, edm::EventSetup const&iSetup);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&iLumiBlock, edm::EventSetup const&iSetup);

  

  void   displayConfig(void);

  // keep some Monte Carlo information if necessary
  void   fillGenJets(  void);
  void   fillAKGenJets(void);
  void   fillSimTracks(void);


  void   fillGenEventInfo(edm::Handle<GenEventInfoProduct> &genEvtInfo, _gen_eventInfo_ *myGenEvtInfo);

  void   copyLHEweights(_event_ *myevt,  const LHEEventProduct * LHEevt);


  void   fillMCInfo(Handle<edm::HepMCProduct>  &mcTruth,   _mc_process_ *mc);
  void   fillMCInfo(Handle<reco::GenParticleCollection> &genParticles,  _mc_process_ *mc);
  void   fillGenTTbar(Handle<reco::GenParticleCollection> &genParticles,  _gen_ttbar_ *genttbar);
  void   fillGenDrellYan(Handle<reco::GenParticleCollection> &genParticles, const LHEEventProduct * evt,  _gen_DrellYan_ *gendrellyan);


  const Candidate *bornLevelParticle( const Candidate *init_p, bool address_down, bool match_initId=true, bool isdebug=false);
  const Candidate *genLevelLeptons(   const Candidate *born_level, math::PtEtaPhiMLorentzVector &dressed);


  void   fillGenWZ(Handle<edm::HepMCProduct> &mcTruth, _genwz_*);
  void   fillGenWZ(Handle<reco::GenParticleCollection> &genParticles, _genwz_*);


  void   fillRunInfo(  void);


  // keep selected event information
  void   fillEventInfo(      const edm::Event& iEvent,  const edm::EventSetup& iSetup);
  void   copyBeamSpotInfo(   const reco::BeamSpot *aBeamSpot, 
			     _beam_spot_ * myBeamSpot);
  void   copyTrackInfo(      const reco::Track *fwTrack, 
			     _track_ *myTrack);
  void   copyMuonInfo(       reco::MuonCollection::const_iterator fwMuon,  
			     _muon_  *myMuon);

  void   copyCaloJetInfo(    const edm::Event& iEvent,
			     const edm::EventSetup& iSetup, 
			     reco::CaloJetCollection::const_iterator jet,
			     edm::RefToBase<reco::Jet> &jetRef, 
			     JetCorrectionUncertainty *jetCorUnc, 
			     double scale, 
			     edm::Handle<reco::JetFlavourMatchingCollection> & theRecoTag, 
			     //   edm::Handle<reco::JetTagCollection> & jetTags, 
			     _jet_ *myJet);

  void   copyJPTJetInfo(     const edm::Event& iEvent,
			     const edm::EventSetup& iSetup, 
			     reco::JPTJetCollection::const_iterator jet, 
			     edm::RefToBase<reco::Jet> &jetRef, 
			     JetCorrectionUncertainty *jetCorUnc, 
			     double scale, 
			     edm::Handle<reco::JetFlavourMatchingCollection> & theRecoTag, 
			     //			     edm::Handle<reco::JetTagCollection> & jetTags, 
			     _jet_ *myJet);


  void   copyPFCHSJetInfo(      const edm::Event& iEvent,
				const edm::EventSetup& iSetup, 
				reco::PFJetCollection::const_iterator jet, 
				edm::RefToBase<reco::Jet> &jetRef, 
				JetCorrectionUncertainty *jetCorUnc, 
				double scale, 
				edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
				_jet_ *myJet);

  void   copyPFJetInfo(      const edm::Event& iEvent,
			     const edm::EventSetup& iSetup, 
			     reco::PFJetCollection::const_iterator jet, 
			     edm::RefToBase<reco::Jet> &jetRef, 
			     JetCorrectionUncertainty *jetCorUnc, 
			     double scale, 
			     edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
			     _jet_ *myJet);
  void   copyPFJetInfoCommon(      const edm::Event& iEvent,
				   const edm::EventSetup& iSetup, 
				   reco::PFJetCollection::const_iterator jet, 
				   edm::RefToBase<reco::Jet> &jetRef, 
				   JetCorrectionUncertainty *jetCorUnc, 
				   double scale, 
				   edm::Handle<reco::JetFlavourInfoMatchingCollection> & theRecoTag, 
				   _jet_ *myJet);
  

  void copySuperclusterInfo( reco::SuperClusterCollection::const_iterator supercluster, _supercluster_ *mySupercluster); 


  void   copyPhotonInfo(     reco::PhotonCollection::const_iterator photon, 
			     _photon_ *myPhoton);
  void   copyElectronInfo(   reco::GsfElectronCollection::const_iterator electron,
			     _electron_ *myElectron,
			     reco::GsfElectronRef electronRef, 
			     int ele_index = -1);
  void   copyVertexInfo(     const reco::VertexCollection::const_iterator vertex,   const reco::MuonCollection *muons,  int leadingMuonIndex,  
			     _vertex_ *myVertex);
  void   copyTrgBits(        const edm::Event& iEvent, 
			     _trg_bits_ * trgBits);
  void   copyHLTInfo(        const edm::Event& iEvent, 
			     _hlt_info_ * hltInfo);
  //  void   copyMETs(           _mets_ *mymet);
  void   copyMET(            const edm::Event& iEvent, const char *mettype, _met_ &amet);


 protected:

  // control flags from the configuration files
  bool              _is_debug;  // output debug information
  bool              _is_data;   // data or MC
  bool              _gen_only;  // dealing with gen events only
  bool              _save_allevents; // determine if it is signal MC
  bool              _vertexing; // determine if performing vertex fitting
  bool              _smoothing; // determine if performing vertex fitting
  //double vertexingMaxDistance;
  //int vertexingMaxNumOfIterations;
  edm::ParameterSet _kvfPSet ;

  bool              _reco;      // determine if there is any need for reco/refit
  bool              _is_save;
  std::string       _reco_selection; // determine if performing the dilepton reco

  float              avgInstLumi;

  double             bField;

  edm::InputTag      BeamSpotTags_;
  std::string        Vertices_;
  edm::InputTag      MuonCollectionTags_;
  edm::InputTag      ElectronCollectionTags_;

  EffectiveAreas _effectiveAreas;


  edm::EDGetTokenT<edm::ValueMap<bool> > EleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > EleTightIdMapToken_;
  
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 


  edm::Handle<reco::ConversionCollection> conversions_h;


  edm::EDGetTokenT<edm::ValueMap<float> > ElectronEcalPFClusterIsolationProducerToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > ElectronHcalPFClusterIsolationProducerToken_;

  edm::Handle<edm::ValueMap<float> >      electronEcalPFClusterIsolation;
  edm::Handle<edm::ValueMap<float> >      electronHcalPFClusterIsolation;

  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float> > TrigMvaValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> >   TrigMvaCategoriesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  TrigMvaMediumIdMapsToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  TrigMvaTightIdMapsToken_;


  
  edm::EDGetTokenT<edm::ValueMap<float> > NonTrigMvaValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> >   NonTrigMvaCategoriesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  NonTrigMvaMediumIdMapsToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  NonTrigMvaTightIdMapsToken_;



  edm::Handle<edm::ValueMap<float> >      trigMvaValues;
  edm::Handle<edm::ValueMap<int> >        trigMvaCategories;
  edm::Handle<edm::ValueMap<bool> >       trigMvaMedium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >       trigMvaTight_id_decisions;


  edm::Handle<edm::ValueMap<float> >      nonTrigMvaValues;
  edm::Handle<edm::ValueMap<int> >        nonTrigMvaCategories;
  edm::Handle<edm::ValueMap<bool> >       nonTrigMvaMedium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >       nonTrigMvaTight_id_decisions;
   

  /*************************************************************************
   *
   * jet flavor, btagging, etc. for 74x
   *
   *************************************************************************/
  edm::InputTag                                            PFCHSJetTags_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> pfchsJetFlavourInfosToken_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection>      thePFCHSJetFlavourInfos;
  std::vector< std::string  >                              PFCHSJetTagInfos_;
  edm::Handle<reco::JetTagCollection>                      myPFCHSJetTagsJP;
  edm::Handle<reco::JetTagCollection>                      myPFCHSJetTagsTCHP;
  edm::Handle<reco::JetTagCollection>                      myPFCHSJetTagsCSV;

  
  edm::InputTag                                            CaloJetTags_;
  edm::Handle<reco::JetFlavourMatchingCollection>          theRecoCaloTag;
  edm::InputTag                                            recoCaloTag_;
  std::vector< std::string  >                              CaloJetTagInfos_;
  edm::Handle<reco::JetTagCollection>                      myCaloJetTagsJP;
  edm::Handle<reco::JetTagCollection>                      myCaloJetTagsTCHP;
  edm::Handle<reco::JetTagCollection>                      myCaloJetTagsCSV;

   
  
  edm::InputTag                                            JPTJetTags_;
  edm::Handle<reco::JetFlavourMatchingCollection>          theRecoJPTTag;
  edm::InputTag                                            recoJPTTag_;
  std::vector< std::string  >                              JPTJetTagInfos_;
  edm::Handle<reco::JetTagCollection>                      myJPTJetTagsJP;
  edm::Handle<reco::JetTagCollection>                      myJPTJetTagsTCHP;
  edm::Handle<reco::JetTagCollection>                      myJPTJetTagsCSV;

  
  edm::InputTag                                            PFJetTags_; 
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> pfJetFlavourInfosToken_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection>      thePFJetFlavourInfos;
  std::vector< std::string  >                              PFJetTagInfos_;
  edm::Handle<reco::JetTagCollection>                      myPFJetTagsJP;
  edm::Handle<reco::JetTagCollection>                      myPFJetTagsTCHP;
  edm::Handle<reco::JetTagCollection>                      myPFJetTagsCSV;


  std::vector< std::string  >                              JetTagCollectionTags_;
  edm::Handle<reco::JetTagCollection>                      jetTags;
  edm::Handle<reco::JetTagCollection>                      jetTagsCSV;

   
  double                                                   jetMinPt;
  double                                                   leptonThreshold;


  // jet id
  edm::InputTag                                            inputJetIDValueMap;
  edm::EDGetTokenT<edm::ValueMap <reco::JetID> >           jetID_ValueMapToken_;
  edm::Handle< edm::ValueMap<reco::JetID> >                jetID_ValueMap_Handle;


  std::string                                              pfchsJetCorrectionService;
  std::string                                              caloJetCorrectionService;
  std::string                                              jptJetCorrectionService;
  std::string                                              pfJetCorrectionService;

  edm::ESHandle<JetCorrectorParametersCollection>          pfchsJetCorParColl;
  edm::ESHandle<JetCorrectorParametersCollection>          caloJetCorParColl;
  edm::ESHandle<JetCorrectorParametersCollection>          jptJetCorParColl;
  edm::ESHandle<JetCorrectorParametersCollection>          pfJetCorParColl;


  const JetCorrector*                                      pfchsJetCorr;
  const JetCorrector*                                      caloJetCorr;
  const JetCorrector*                                      jptJetCorr;
  const JetCorrector*                                      pfJetCorr;


  JetCorrectionUncertainty                                *pfchsJetUnc;
  JetCorrectionUncertainty                                *caloJetUnc;
  JetCorrectionUncertainty                                *jptJetUnc;
  JetCorrectionUncertainty                                *pfJetUnc;


  // default rho and sigma
  edm::EDGetTokenT<double>   FixGridRhoToken_;
  edm::Handle<double>        fixGridRhoHandle;


  edm::InputTag              RhoSrcLabel_;
  edm::Handle<double>        rhoHandle;

  edm::InputTag              SigmaSrcLabel_;
  edm::Handle<double>        sigmaHandle;

  // CHS
  edm::InputTag              RhoSrcLabelCHS_;
  edm::Handle<double>        rhoHandleCHS;

  edm::InputTag              SigmaSrcLabelCHS_;
  edm::Handle<double>        sigmaHandleCHS;

  // Calo
  edm::InputTag              RhoSrcLabelCalo_;
  edm::Handle<double>        rhoHandleCalo;

  edm::InputTag              SigmaSrcLabelCalo_;
  edm::Handle<double>        sigmaHandleCalo;


  // Track
  edm::InputTag              RhoSrcLabelTrack_;
  edm::Handle<double>        rhoHandleTrack;

  edm::InputTag              SigmaSrcLabelTrack_;
  edm::Handle<double>        sigmaHandleTrack;


  // isolation rho and sigma
  edm::InputTag              RhoIsoSrcLabel_;
  edm::Handle<double>        rhoIsoHandle;

  edm::InputTag              SigmaIsoSrcLabel_;
  edm::Handle<double>        sigmaIsoHandle;


  //charge hadron rho and sigma
  // eta is within 2.0
  edm::InputTag              RhoChSrcLabel_;
  edm::Handle<double>        rhoChHandle;


  edm::InputTag              SigmaChSrcLabel_;
  edm::Handle<double>        sigmaChHandle;

  // eta up to 2.4
  edm::InputTag              RhoCh2p4SrcLabel_;
  edm::Handle<double>        rhoCh2p4Handle;


  edm::InputTag              SigmaCh2p4SrcLabel_;
  edm::Handle<double>        sigmaCh2p4Handle;





  edm::InputTag      TrackCollectionTags_;
  double             trackMinPtWithMCTruth;
  double             leptonMinPtForComposition;
  edm::InputTag      PhotonCollectionTags_;
  edm::InputTag      L1ParticleMapCollectionTags_;
  edm::InputTag      L1GTReadoutRecordLabel_;
  edm::InputTag      HLTL1GTObjectMapLabel_;
  edm::InputTag      TriggerResultsLabel_;
  edm::InputTag      TriggerSummaryLabel_;


  std::vector< std::string>  triggerNames;


  // for HLT triggers
  unsigned int       refIndex ;
  std::string        TriggerRefPath_;
  std::string        TriggerProcess_;


  std::vector< std::string  >       HLTTriggerMuons_;
  std::vector< std::string  >       HLTTriggerMuonElectrons_;
  std::vector< std::string  >       HLTTriggerElectrons_;

  std::vector< std::string  >       tmpHLTTriggerMuons_;
  std::vector< std::string  >       tmpHLTTriggerMuonElectrons_;
  std::vector< std::string  >       tmpHLTTriggerElectrons_;


  std::string        HLTTriggerMuon_;
  std::string        HLTTriggerMuonLongName_;
  std::string        HLTTriggerMuonL_;
  std::string        HLTTriggerMuonLLongName_;
  std::string        HLTTriggerElectron_;
  std::string        HLTTriggerElectronLongName_;
  std::string        HLTTriggerElectronL_;
  std::string        HLTTriggerElectronLLongName_;


  edm::InputTag      GeneratorLevelTag_;
  edm::InputTag      LHEEventProductTag_;

  edm::InputTag              GenJetAlgorithmTags_;
  edm::InputTag              akGenJetAlgorithmTags_;

  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genJetFlavourInfosToken_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection>      theGenJetFlavourInfos;


  double                     genJetMinPt;


  edm::InputTag              SimTrackTags_;
  double                     diLeptonMinMass;

  //  std::string                out, open;
  // std::string                pdf;
  // int                        subset;
  double                     xmin, xmax, qmin, qmax;



  bool                                            hasGenJets;

  Handle<edm::HepMCProduct>                       mcTruth;
  // Handle<reco::GenParticleCollection>          genParticlesAOD;
  Handle<reco::GenParticleCollection>             genParticles;
  edm::Handle<GenEventInfoProduct>                genEventInfo;

  Handle<LHEEventProduct>                         lheEventInfo;


  //  handles to access collection product. 
  //  Handle< L1ParticleMapCollection >           l1MapColl ;
  Handle<reco::PhotonCollection>                  photons; 
  Handle<reco::SuperClusterCollection>            superclusters; 
  Handle<reco::MuonCollection>                    muons;
  Handle<reco::GsfElectronCollection>             electrons;
  edm::Handle<reco::ConversionCollection>         convCol;
  Handle<reco::GsfElectronCollection > calibratedElectrons;
  edm::Handle<edm::ValueMap<double>> regEne_handle;
  edm::Handle<edm::ValueMap<double>> regErr_handle;
 
  edm::Handle<edm::ValueMap<float>> mvaTrigV0_handle;
  edm::Handle<edm::ValueMap<float>> mvaNonTrigV0_handle;
   





  Handle<GenJetCollection>                        genJets;
  Handle<GenJetCollection>                        akGenJets;
  Handle<reco::PFJetCollection>                   pfchsJets;
  Handle<reco::CaloJetCollection>                 caloJets;
  Handle<reco::JPTJetCollection>                  jptJets;
  Handle<reco::PFJetCollection>                   pfJets;

  Handle<edm::SimTrackContainer>                  simTracks;
  Handle<reco::TrackCollection>                   tracks;



  Handle<CaloTowerCollection>                     towers;

  edm::Handle<L1GlobalTriggerReadoutRecord>       l1TriggerReadoutRec;
  edm::Handle<L1GlobalTriggerObjectMapRecord>     hltL1GTObjMap; 
  edm::ESHandle<L1GtTriggerMenu>                  l1GTTriggerMenu;
  edm::Handle<TriggerResults>                     triggerResults;
  edm::Handle<trigger::TriggerEvent>              triggerObj;  


  // pile up
  Handle<std::vector< PileupSummaryInfo > >       PupInfo;

  //
  // mets
  edm::Handle< GenMETCollection >      genMEThandle;

  //  edm::Handle< CaloMETCollection >     rawMEThandle;
  //  edm::Handle< METCollection >         tcMEThandle;
  //  edm::Handle< CaloMETCollection >     muCorrMEThandle;
  //  edm::Handle< PFMETCollection >       pfMEThandle;
  //  edm::Handle< CaloMETCollection >     muJESCorrMEThandle;



  // tracker geometry
  edm::ESHandle<TrackerGeometry>       tracker;

  // beam spot
  edm::Handle<reco::BeamSpot>          recoBeamSpotHandle;

  edm::Handle<EcalRecHitCollection>    ecalEBRecHitHandle;
  edm::Handle<EcalRecHitCollection>    ecalEERecHitHandle;
  edm::Handle<EcalRecHitCollection>    ecalRecHitHandle;

  Handle< int >                        genProcessID;
  Handle< double >                     genEventScale;
  Handle< GenRunInfoProduct >          runInfo;
  Handle< double >                     evtWeight;
  Handle<reco::VertexCollection>       recVtxs;

  const reco::Vertex                   *recVtx;
  const reco::BeamSpot                 *vertexBeamSpot;
  const HepMC::GenEvent                *genEvt;
  const CaloTowerCollection            *caloTowers;
  const reco::TrackCollection          *recoTracks, *recoIsoTracks;
  const reco::MuonCollection           *recoMuons;
  const reco::GsfElectronCollection    *recoElectrons;
  const reco::MuonCollection           *recoTrackerMuons;
  const edm::SimTrackContainer         *recoSimTracks;
  const reco::PFJetCollection          *recoJets;
  const reco::JetTagCollection         *recoJetTags;
  const reco::PhotonCollection         *recoPhotons;
  //  const reco::modules::JetFlavourIdentifier     *jetFlavorIdentifier;


  edm::Handle<reco::JetFlavourMatchingCollection> theGenTag;
  edm::Handle<reco::JetFlavourMatchingCollection> theRecoPFTag;





  // TeV Muon Refit
  Handle <reco::TrackToTrackMap>       tevMapH1;   
  reco::TrackToTrackMap                tevMap1;
  Handle <reco::TrackToTrackMap>       tevMapH2;
  reco::TrackToTrackMap                tevMap2;
  Handle <reco::TrackToTrackMap>       tevMapH3;
  reco::TrackToTrackMap                tevMap3;



  // user code for quark-gluon likelihood separator
  // QGLikelihoodCalculator              *qgLikelihoodCal;
  edm::Handle<edm::ValueMap<float> >  QGTagsHandleMLP;
  edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;



  
  
  // output data defintion/access
  // necessary global variables
  edm::Service< TFileService > fs;
  //  TFile                    *hFile;
  TTree                    *ntuple;
  _event_                  *myEvent;        // event setup
  

  // these are not really needed in global
  _photon_                 *myPhoton;       // photon list
  _electron_               *myElectron;
  _track_                  *myTrack;        // track 
  _track_                  *mySimTrack;     // sim tracks-structure same as track
  _muon_                   *myMuon;         // muon
  _W_                      *myW;            // W candidate, one to one to muon
  _jet_                    *myJet;          // jet
  _di_lepton_              *myDilepton;     // lepton pairs (mumu only)
  _l1_obj_                 *myL1Obj;        // l1 objects
  _mc_process_             *myMCTruth;      // MC production
  _genwz_                  *myGenWZ;        // MC WZ Information
  _run_info_               *myRunInfo;      // run information


  //  some other variables 
  const char               *lhaPDFPath;
  Int_t                     evtNum, runNum;
  Int_t                     totalProcessedEvts;


  // vertex fitter
  KalmanVertexFitter                        *fitter;
  vector<TransientTrack>                     vertexTracks;
  edm::ESHandle<TransientTrackBuilder>       transientTrackBuilder;

  // hlt information
  HLTConfigProvider                          hltConfig;

  // additional functions
  // void copyHLTInfo(const edm::Event& iEvent, 
  //	   edm::Handle<TriggerResults>         &triggerResults, 
  //	   std::string &triggerProcessName, 
  //	   edm::Handle<trigger::TriggerEvent>  &triggerObj, 
  //	   HLTConfigProvider &hltConfig,
  //	   bool &trgBit, const char *hltname,  
  //	   _event_ *aEvent,  
  //	   _vec4_ *(*ptr2AddHlt)(_event_ *));


  // additional functions
  void copyHLTInfo(const edm::Event& iEvent, 
		   edm::Handle<TriggerResults>         &triggerResults, 
		   std::string &triggerProcessName, 
		   edm::Handle<trigger::TriggerEvent>  &triggerObj, 
		   HLTConfigProvider &hltConfig,
		   bool &trgBit, const char *hltname,  
		   _event_ *aEvent,  
		   _vec4_ *(*ptr2AddHlt)(_event_ *), 
		   _vec4_ *(*ptr2AddHltLeg1)(_event_ *) = NULL, 
		   _vec4_ *(*ptr2AddHltLeg2)(_event_ *) = NULL) ;
  

  void checkHLTVersion(std::vector< std::string > & triggerNames, std::vector<std::string> &atmpset, std::vector<std::string> &aset);

//  void copyHLTInfo(const edm::Event& iEvent, bool &trgBit, const char *hltname, HLTConfigProvider &hltConfig , _event_ *aEvent, _vec4_ *(*ptr2AddHlt)(_event_ *) );


};

#endif
