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
#include "EdmAnalysisTools.h"
#include "lhapdfcc.h"


using namespace std;

using namespace edm;
using namespace reco;
using namespace l1extra ;
using namespace lhef;

WZEdmAnalyzer::WZEdmAnalyzer(const edm::ParameterSet& iConfig) :
  _is_debug(                   iConfig.getParameter<bool>("DEBUG")),
  _is_data(                    iConfig.getParameter<bool>("DATA")),
  _gen_only(                   iConfig.getParameter<bool>("GEN_ONLY")),
  _mc_signal(                  iConfig.getParameter<bool>("MC_SIGNAL")),
  _vertexing(                  iConfig.getParameter<bool>("VERTEXING")),
  _smoothing(                  iConfig.getParameter<bool>("SMOOTHING")),
  _kvfPSet(                    iConfig.getParameter<edm::ParameterSet>("KVFParameters")),
  _reco(                       iConfig.getParameter<bool>("RECO")),
  _reco_selection(             iConfig.getParameter<std::string>("RECOSELECTION")),
  BeamSpotTags_(               iConfig.getParameter<edm::InputTag>("BeamSpot")), 
  Vertices_(                   iConfig.getParameter<std::string>("Vertices")),
  MuonCollectionTags_(         iConfig.getParameter<std::string>("Muons")),
  ElectronCollectionTags_(     iConfig.getParameter<std::string>("Electrons")),
  _effectiveAreas(            (iConfig.getParameter<edm::FileInPath>("EffAreasConfigFile")).fullPath()),
  EleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(  iConfig.getParameter<edm::InputTag>("EleLooseIdMap"))),
  EleMediumIdMapToken_(consumes<edm::ValueMap<bool> >( iConfig.getParameter<edm::InputTag>("EleMediumIdMap"))),
  EleTightIdMapToken_(consumes<edm::ValueMap<bool> >(  iConfig.getParameter<edm::InputTag>("EleTightIdMap"))),
  JetTags_(                    iConfig.getParameter<std::string>("Jets")),
  CaloJetTags_(                iConfig.getParameter<std::string>("CaloJets")),
  JPTJetTags_(                 iConfig.getParameter<std::string>("JPTJets")),
  PFJetTags_(                  iConfig.getParameter<std::string>("PFJets")),
  jetMinPt(                    iConfig.getParameter<double>("JetMinPt")),
  leptonThreshold(             iConfig.getParameter<double>("LeptonThreshold")),
  inputJetIDValueMap(          iConfig.getParameter<edm::InputTag>("InputJetIDValueMap")), 
  // calojetIDHelperConfig(       iConfig.getParameter<edm::ParameterSet>( "JetIDParams" ) ), 
  jetCorrectionService(        iConfig.getParameter<std::string>("JetCorrectionService") ), 
  caloJetCorrectionService(    iConfig.getParameter<std::string>("CaloJetCorrectionService") ), 
  jptJetCorrectionService(     iConfig.getParameter<std::string>("JPTJetCorrectionService") ), 
  pfJetCorrectionService(      iConfig.getParameter<std::string>("PFJetCorrectionService") ), 
  FixGridRhoToken_(consumes<double> ( iConfig.getParameter<edm::InputTag>("FixGridRho"))), 
  RhoSrcLabel_(                iConfig.getParameter<edm::InputTag>("RhoSrc")), 
  SigmaSrcLabel_(              iConfig.getParameter<edm::InputTag>("SigmaSrc")), 
  RhoIsoSrcLabel_(             iConfig.getParameter<edm::InputTag>("RhoIsoSrc")), 
  SigmaIsoSrcLabel_(           iConfig.getParameter<edm::InputTag>("SigmaIsoSrc")), 
  RhoChSrcLabel_(              iConfig.getParameter<edm::InputTag>("RhoChSrc")), 
  SigmaChSrcLabel_(            iConfig.getParameter<edm::InputTag>("SigmaChSrc")), 
  RhoCh2p4SrcLabel_(           iConfig.getParameter<edm::InputTag>("RhoCh2p4Src")), 
  SigmaCh2p4SrcLabel_(         iConfig.getParameter<edm::InputTag>("SigmaCh2p4Src")), 
  TrackCollectionTags_(        iConfig.getParameter<std::string>("Tracks")),
  trackMinPtWithMCTruth(       iConfig.getParameter<double>("TrackMinPtWithMCTruth")),
  leptonMinPtForComposition(   iConfig.getParameter<double>("LeptonMinPtForComposition")),
  JetTagCollectionTags_(       iConfig.getParameter<std::string>("JetTagCollection")),
  PhotonCollectionTags_(       iConfig.getParameter<std::string>("PhotonCollection")),
  L1ParticleMapCollectionTags_(iConfig.getParameter<std::string>("L1ParticleMapCollection")),
  L1GTReadoutRecordLabel_(     iConfig.getParameter<edm::InputTag>("L1GTReadoutRecordLabel")),
  HLTL1GTObjectMapLabel_(      iConfig.getParameter<edm::InputTag>("HLTL1GTObjectMapLabel")),
  TriggerResultsLabel_(        iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  TriggerSummaryLabel_(        iConfig.getParameter<edm::InputTag>("TriggerSummaryLabel")),
  TriggerRefPath_(             iConfig.getParameter<std::string>("TriggerRefPath")),
  tmpHLTTriggerMuons_(         iConfig.getParameter<std::vector< std::string> >("HLTTriggerMuons")), 
  tmpHLTTriggerMuonElectrons_( iConfig.getParameter<std::vector< std::string> >("HLTTriggerMuonElectrons")), 
  tmpHLTTriggerElectrons_(     iConfig.getParameter<std::vector< std::string> >("HLTTriggerElectrons")), 
  //  HLTTriggerMuon_(             iConfig.getParameter<std::string>("HLTTriggerMuon")),
  //  HLTTriggerMuonL_(            iConfig.getParameter<std::string>("HLTTriggerMuonL")),
  //  HLTTriggerElectron_(         iConfig.getParameter<std::string>("HLTTriggerElectron")),
  //  HLTTriggerElectronL_(        iConfig.getParameter<std::string>("HLTTriggerElectronL")),
  GeneratorLevelTag_(          iConfig.getParameter<std::string>("GeneratorLevelTag")),
  LHEEventProductTag_(         iConfig.getParameter<edm::InputTag>("LHEEventProductTag")),
  GenJetAlgorithmTags_(        iConfig.getParameter<std::string>("GenJets")),
  akGenJetAlgorithmTags_(      iConfig.getParameter<std::string>("akGenJets")),
  genJetMinPt(                 iConfig.getParameter<double>("GenJetMinPt")),
  SimTrackTags_(               iConfig.getParameter<std::string>("SimTracks")),
  diLeptonMinMass(             iConfig.getParameter<double>("DiLeptonMinMass")),
  out(                         iConfig.getParameter<std::string>("out")),
  open(                        iConfig.getParameter<std::string>("open")),
  pdf(                         iConfig.getParameter<std::string>("pdf")),
  subset(                      iConfig.getParameter<int>("subset")) {



  // for accessing trigger information
  TriggerProcess_ = TriggerResultsLabel_.process();


  if (!_gen_only) {
    fitter = new KalmanVertexFitter(_kvfPSet , _smoothing);
    // fitter->setMaximumDistance(vertexingMaxDistance);
    //fitter->setMaximumNumberOfIterations(vertexingMaxNumOfIterations);
  }



  jetID_ValueMapToken_= consumes< edm::ValueMap<reco::JetID> >(inputJetIDValueMap);






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

  if (_is_debug) std::cout  << "check point ... create output file - " << out.c_str() << std::endl;


  // Define output file and ntuple
  hFile  = new TFile( out.c_str(), open.c_str() );
  ntuple = new TTree("ntuple", "");
  
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

      iEvent.getByLabel(GeneratorLevelTag_,     mcTruth);
      this->fillMCInfo(mcTruth,                 myMCTruth);
      this->fillGenWZ(mcTruth,                  myGenWZ);

    } else {


      iEvent.getByLabel( "genParticles",        genParticles );
      bool hasLHE = //iEvent.getByType( lheEventInfo );
	iEvent.getByLabel( LHEEventProductTag_, lheEventInfo );

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
  edm::Handle<bool> trackingFailure;
  iEvent.getByLabel("trackingFailureFilter", trackingFailure);
  myEvent->getEventFilterBit()->trackingFailureFilter =  *(trackingFailure);



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


  iEvent.getByLabel(BeamSpotTags_, recoBeamSpotHandle);              
  if (recoBeamSpotHandle.isValid() ) {
    vertexBeamSpot = recoBeamSpotHandle.product();
  } else {
    vertexBeamSpot = 0;
  } 

  // (++). Acess the jet ID value map
  iEvent.getByToken(jetID_ValueMapToken_,jetID_ValueMap_Handle);


  // (++). Access the jet correction
  jetCorr     = JetCorrector::getJetCorrector(jetCorrectionService,    iSetup);
  caloJetCorr = JetCorrector::getJetCorrector(caloJetCorrectionService,iSetup);
  jptJetCorr  = JetCorrector::getJetCorrector(jptJetCorrectionService, iSetup);
  pfJetCorr   = JetCorrector::getJetCorrector(pfJetCorrectionService,  iSetup);


  // (++) Access the uncertainty of the jet energy correction
  iSetup.get<JetCorrectionsRecord>().get("AK5Calo",jetCorParColl);
  jetUnc = new JetCorrectionUncertainty( (*jetCorParColl)["Uncertainty"]  );


  // pf jet 
  iSetup.get<JetCorrectionsRecord>().get("AK5PF",pfJetCorParColl);
  pfJetUnc = new JetCorrectionUncertainty( (*pfJetCorParColl)["Uncertainty"]  );
  

  // get the rho values from fast jet 
  //edm::Handle<double> rhoHandle;
  // default rho value for ak4PFjet
  // NOTE74

  iEvent.getByToken(FixGridRhoToken_, fixGridRhoHandle);  

  iEvent.getByLabel(RhoSrcLabel_,          rhoHandle);         
  iEvent.getByLabel(SigmaSrcLabel_,        sigmaHandle);         
  if (rhoHandle.isValid() && sigmaHandle.isValid() ) {
    myEvent->rho      = *rhoHandle;
    myEvent->sigma    = *sigmaHandle;
  }

  // default rho value for ak4PFjetCHS
  iEvent.getByLabel(RhoSrcLabelCHS_,          rhoHandleCHS);         
  iEvent.getByLabel(SigmaSrcLabelCHS_,        sigmaHandleCHS);  
  if (rhoHandleCHS.isValid() && sigmaHandleCHS.isValid() ) {

    myEvent->rhoCHS      = *rhoHandleCHS;
    myEvent->sigmaCHS    = *sigmaHandleCHS;
  }

  
  // default rho value for ak4CaloJet
  iEvent.getByLabel(RhoSrcLabelCalo_,          rhoHandleCalo);         
  iEvent.getByLabel(SigmaSrcLabelCalo_,        sigmaHandleCalo);         
  if (rhoHandleCalo.isValid() && sigmaHandleCalo.isValid() ) {

    myEvent->rhoCalo      = *rhoHandleCalo;
    myEvent->sigmaCalo    = *sigmaHandleCalo;
  }

  // default rho value for track jet
  iEvent.getByLabel(RhoSrcLabelTrack_,          rhoHandleTrack);         
  iEvent.getByLabel(SigmaSrcLabelTrack_,        sigmaHandleTrack);         
  if (rhoHandleTrack.isValid() && sigmaHandleTrack.isValid() ) {
    myEvent->rhoTrack      = *rhoHandleTrack;
    myEvent->sigmaTrack    = *sigmaHandleTrack;
  }

  // others
  iEvent.getByLabel(RhoIsoSrcLabel_,       rhoIsoHandle);   
  iEvent.getByLabel(SigmaIsoSrcLabel_,     sigmaIsoHandle);   
  myEvent->rhoIso   = *rhoIsoHandle;
  myEvent->sigmaIso = *sigmaIsoHandle;



  iEvent.getByLabel(RhoChSrcLabel_,       rhoChHandle);   
  iEvent.getByLabel(SigmaChSrcLabel_,     sigmaChHandle);   
  myEvent->rhoCh   = *rhoChHandle;
  myEvent->sigmaCh = *sigmaChHandle;



  // eta up to 2.4
  iEvent.getByLabel(RhoCh2p4SrcLabel_,       rhoCh2p4Handle);   
  iEvent.getByLabel(SigmaCh2p4SrcLabel_,     sigmaCh2p4Handle);   
  myEvent->rhoCh2p4   = *rhoCh2p4Handle;
  myEvent->sigmaCh2p4 = *sigmaCh2p4Handle;



  // (++). Setup for vertex
  // Be careful, maybe do not work well on CMSSW_3XX
  //  if (_vertexing) 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transientTrackBuilder);




  // (++). Extract trigger information
  iEvent.getByLabel(L1GTReadoutRecordLabel_,          l1TriggerReadoutRec);
  iEvent.getByLabel("hltL1GtObjectMap",               hltL1GTObjMap);
  iSetup.get<L1GtTriggerMenuRcd>().get(               l1GTTriggerMenu);


  iEvent.getByLabel(TriggerResultsLabel_,             triggerResults);
  iEvent.getByLabel(TriggerSummaryLabel_,             triggerObj);

  // (++). Define and extract container item
  iEvent.getByLabel(Vertices_,                        recVtxs);
  iEvent.getByLabel(TrackCollectionTags_,             tracks);
  iEvent.getByLabel(MuonCollectionTags_,              muons);
  iEvent.getByLabel(ElectronCollectionTags_,          electrons);
  iEvent.getByLabel(JetTags_,                         jets);
  iEvent.getByLabel(CaloJetTags_,                     caloJets);
  iEvent.getByLabel(JPTJetTags_,                      jptJets);
  iEvent.getByLabel(PFJetTags_,                       pfJets);
  iEvent.getByLabel(PhotonCollectionTags_,            photons);
  iEvent.getByLabel("superClusters",                  superclusters);
  iEvent.getByLabel("allConversions",                 convCol);



  // access cut-based electron ID value map
  iEvent.getByToken(EleLooseIdMapToken_ ,             loose_id_decisions);
  iEvent.getByToken(EleMediumIdMapToken_,             medium_id_decisions);
  iEvent.getByToken(EleTightIdMapToken_,              tight_id_decisions);


  // regression electron calibration
  /*  NOTE74: comment out for now Oct. 19, 2015
  iEvent.getByLabel(edm::InputTag("calibratedElectrons","calibratedGsfElectrons"), calibratedElectrons);
  iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneRegForGsfEle"), regEne_handle);
  iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);


  iEvent.getByLabel("mvaTrigV0", mvaTrigV0_handle);
  iEvent.getByLabel("mvaNonTrigV0", mvaNonTrigV0_handle);

  if (_is_debug) {
    for (reco::GsfElectronCollection::const_iterator electron = (*calibratedElectrons.product()).begin(); electron != (*calibratedElectrons.product()).end(); electron ++) {
      
      std::cout << electron->pt() << ", " << electron->eta() << ", " << electron->phi() << ", " << electron->superCluster()->rawEnergy() << std::endl;
      
    }
  }

  */


  // btagging 
  iEvent.getByLabel(JetTagCollectionTags_,                   jetTags);
  iEvent.getByLabel("combinedSecondaryVertexBJetTags",       jetTagsCSV);

  iEvent.getByLabel("MyJetBProbabilityBJetTags",             myJetTagsJP);
  iEvent.getByLabel("MyTrackCountingHighPurBJetTags",        myJetTagsTCHP);
  iEvent.getByLabel("MyCombinedSecondaryVertexBJetTags",     myJetTagsCSV);


  iEvent.getByLabel("MyJetBProbabilityCaloBJetTags",         myCaloJetTagsJP);
  iEvent.getByLabel("MyTrackCountingHighPurCaloBJetTags",    myCaloJetTagsTCHP);
  iEvent.getByLabel("MyCombinedSecondaryVertexCaloBJetTags", myCaloJetTagsCSV);


  // quark gluon likelihood separator
  iEvent.getByLabel("QGTagger","qgMLP", QGTagsHandleMLP);
  iEvent.getByLabel("QGTagger","qgLikelihood", QGTagsHandleLikelihood);



  // access corrected METs
  //  iEvent.getByLabel("met",                            rawMEThandle);
  //  iEvent.getByLabel("corMetGlobalMuons",              muCorrMEThandle);
  //  iEvent.getByLabel("tcMet",                          tcMEThandle);
  //  iEvent.getByLabel("pfMet",                          pfMEThandle);
  //  iEvent.getByLabel("metMuonJESCorAK5",               muJESCorrMEThandle);

  // ecal hits handle
  // if (_reco) {
  iEvent.getByLabel("reducedEcalRecHitsEB",         ecalEBRecHitHandle);
  iEvent.getByLabel("reducedEcalRecHitsEE",         ecalEERecHitHandle);
    // }

  //  EcalClusterLazyTools *lazyTools(iEvent, iSetup, ecalEBRecHitHandle.product(),ecalEERecHitHandle.product() );



  
  // electron pf isolations
  //  iEvent.getByLabel(ElectronIsoValsTags_[0] ,         electronIsoValsCh);
  //  iEvent.getByLabel(ElectronIsoValsTags_[1] ,         electronIsoValsPhoton);
  //  iEvent.getByLabel(ElectronIsoValsTags_[2] ,         electronIsoValsNeutral); 

  //  edm::Handle<reco::ConversionCollection> conversions_h;
  // iEvent.getByLabel("allConversions", conversions_h);

 



  recoPhotons = 0;
  recoTracks  = 0;
  recoMuons   = 0;
  recoJets    = 0;
  caloTowers  = 0;

  recoTracks       =   tracks.product();
  recoMuons        =   muons.product();
  recoJets         =   jets.product();
  //  caloTowers       =   towers.product();
  recoPhotons      =   photons.product();
  recoPhotons      =   photons.product(); 
  recoElectrons    =   electrons.product(); 
  recoJetTags      =   jetTags.product();


  // TeV Muon Refit
  iEvent.getByLabel("tevMuons", "default", tevMapH1); tevMap1 = *(tevMapH1.product());
  iEvent.getByLabel("tevMuons", "firstHit", tevMapH2);tevMap2 = *(tevMapH2.product());
  iEvent.getByLabel("tevMuons", "picky", tevMapH3);   tevMap3 = *(tevMapH3.product());
  

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

  /*
  if (_is_data) {

    _is_save = true; // no any filter this time

  } else {  // MC background filter events with at least one electron or muon 

    _is_save = false;
    for (reco::GsfElectronCollection::const_iterator electron = (*recoElectrons).begin(); electron != (*recoElectrons).end(); electron ++) {
      
      if ( electron->energy() > leptonThreshold) _is_save = true;
    }


    for (reco::MuonCollection::const_iterator muon = (*recoMuons).begin(); 
	 muon != (*recoMuons).end(); 
	 muon ++) {
      
      if (muon->pt() > leptonThreshold) _is_save = true;

    }

    _is_save = true;
    if (!_is_save && !_mc_signal) return;
  }

  */

  
  
  // (++). Process Monte Carlo events  
  if (!_is_data)   { 

    // pile up reweighting
    if ( iEvent.getByLabel("addPileupInfo", PupInfo) ) {
      

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
    

    hasGenJets = iEvent.getByLabel( GenJetAlgorithmTags_,    genJets );
    iEvent.getByLabel( akGenJetAlgorithmTags_,               akGenJets );    
    // iEvent.getByLabel( "genParticleCandidates",           genParticles );
    // iEvent.getByLabel(GeneratorLevelTag_,                 mcTruth);
    // iEvent.getByLabel(SimTrackTags_,                      simTracks);

    iEvent.getByLabel ("flavourByValGenJet" ,                theGenTag);
    iEvent.getByLabel ("flavourByValPF",                     theRecoPFTag);
    iEvent.getByLabel ("flavourByValCalo",                   theRecoCaloTag);
    //    iEvent.getRun().getByLabel("source",                     runInfo); 
    iEvent.getByLabel("generator",                           genEventInfo);

    myEvent->setEventWeight(genEventInfo->weight());



    iEvent.getByLabel("genMetTrue",                   genMEThandle);

    myEvent->getMETs()->genMET.pt                  = (genMEThandle->front() ).et();
    myEvent->getMETs()->genMET.phi                 = (genMEThandle->front() ).phi();
    myEvent->getMETs()->genMET.sumEt               = (genMEThandle->front() ).sumEt();
    myEvent->getMETs()->genMET.e_longitudinal      = (genMEThandle->front() ).e_longitudinal();
    myEvent->getMETs()->genMET.MuonEtFraction      = (genMEThandle->front() ).MuonEtFraction();
    myEvent->getMETs()->genMET.InvisibleEtFraction = (genMEThandle->front() ).InvisibleEtFraction();


    // Feb 07, 2013: only the aod version is being kept, in general _reco is false. 
    if (_reco) {

      iEvent.getByLabel(GeneratorLevelTag_,           mcTruth);
      iEvent.getByLabel(SimTrackTags_,                simTracks);
      recoSimTracks                                =  simTracks.product();

      this->fillMCInfo(mcTruth,                       myMCTruth);
      this->fillGenWZ(mcTruth,                        myGenWZ);
      this->fillSimTracks();

    } else {
      
      iEvent.getByLabel( "genParticles",              genParticles );
      bool hasLHE = iEvent.getByLabel("source", lheEventInfo );
      //    iEvent.getByLabel( "source", lheEventInfo );
       
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

    //  this->fillGenJets();
    this->fillAKGenJets();
  }


  // (++). Out put cross section information
  this->fillRunInfo();  
  if (_is_debug) std::cout << "checking point .. finish filling run information" << std::endl;

  

  if (_is_debug) std::cout << "checking point .. finish candidate composition" << std::endl;


  // (++). Copy some of the event information for both real data and MC
  this->fillEventInfo(iEvent, iSetup);



  if (_is_save)  ntuple->Fill();


  delete jetUnc;
  delete pfJetUnc;
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

  std::cout <<"Trigger Results: " << TriggerResultsLabel_.process() <<std::endl;

  if (hltConfig.init(iRun,iSetup, TriggerResultsLabel_.process(), changed)) {
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
  
  hFile ->cd();
  ntuple->Write();
  hFile ->Write();
  hFile ->Close();

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
    edm::Handle<LumiSummary> l;
    iLumiBlock.getByLabel("lumiProducer", l); 
    // Check that there is something
    if (l.isValid())  avgInstLumi=l->avgInsDelLumi();
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
  std::cout <<right<< setw(12) <<"*"<<setw(25)<< left<< "  Signal MC flag -- "  << setw(5)<< _mc_signal   <<"*"<< std::endl;
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
