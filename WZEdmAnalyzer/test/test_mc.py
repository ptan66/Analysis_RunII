import FWCore.ParameterSet.Config as cms
import string



process = cms.Process("asym")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("./mymessage.log")
process.MessageLogger.suppressInfo = cms.untracked.vstring("ERROR")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
             reportEvery = cms.untracked.int32(1000),
             limit       = cms.untracked.int32(10000000) )

#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')


process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")



## global tag for MC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = 'MCRUN2_74_V9'
#process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'
## global tag for 2015C
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'
## global tag for 2015D
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'



process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#######################################################################
#
#  global control flag
#
#######################################################################
useAOD = True



process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource", 
       	fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root'


    )
)


# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    fileName = cms.untracked.string('step2_RAW2DIGI_L1Reco_RECO_PU.root'),
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('AOD'),
#        filterName = cms.untracked.string('')
#    )
#)
#process.out_step = cms.EndPath(process.output)





###########################################################################
#
# common filters suggested to be included
# Feb08, 2013
#
###########################################################################
process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')



## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)


## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

## The tracking POG filters __________________________________________________||
process.load('RecoMET.METFilters.trackingPOGFilters_cff')
## NOTE: to make tagging mode of the tracking POG filters (three of them), please do:
process.manystripclus53X.taggedMode = cms.untracked.bool(True)
process.manystripclus53X.forcedValue = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
## Also the stored boolean for the three filters is opposite to what we usually
## have for other filters, i.e., true means rejected bad events while false means 
## good events.


###############################################################
#
#
#  jet correction
#
###############################################################
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#from RecoJets.Configuration.RecoPFJets_cff import *
#from JetMETCorrections.Configuration.DefaultJEC_cff import *

from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectors_cff import *



#process.myjecs = cms.Sequence(
#    process.ak4PFJetsL1FastL2L3
#)

#for isolation purpose???
process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)


# calculate rho within tracker coverager. 
# from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedParticles_cfi import pfAllChargedParticles
process.pfPileUpAllChargedParticlesClone =  cms.EDFilter("PdgIdPFCandidateSelector",
     src = cms.InputTag("particleFlow"),
     pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212,11,-11,13,-13)
  #  makeClones = cms.bool(True)
   )


process.kt6PFJetsForCh =  process.kt6PFJets.clone(
    src = cms.InputTag("pfPileUpAllChargedParticlesClone"),
    Ghost_EtaMax = cms.double(2.6),
    Rho_EtaMax = cms.double(2.0)
    )

process.kt6PFJetsForCh2p4 =  process.kt6PFJets.clone(
    src = cms.InputTag("pfPileUpAllChargedParticlesClone"),
    Ghost_EtaMax = cms.double(3.0),
    Rho_EtaMax = cms.double(2.4)
    )



#NOTE74
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
#process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
#process.corrPfMetType1.jetCorrLabel = cms.InputTag('ak4PFCHSL1FastL2L3')
#process.corrPfMetType1.src='ak4PFJetsCHS'
#process.corrPfMetType1.skipMuons=False
#process.mymets=cms.Sequence(    
#    process.corrPfMetType1 +
#    process.correctionTermsPfMetType0PFCandidate + 
#    process.correctionTermsPfMetShiftXY +
#    process.pfMetT0pc +
#    process.pfMetT0pcT1 +
#    process.pfMetT0pcTxy +
#    process.pfMetT0pcT1Txy +
#    process.pfMetT1 +
#    process.pfMetT1Txy
#)



#btagging
process.MyAk4PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
   process.j2tParametersVX,
   jets = cms.InputTag("ak4PFJets")
)
process.MyImpactParameterPFTagInfos = process.impactParameterTagInfos.clone(
  jetTracks = "MyAk4PFJetTracksAssociatorAtVertex"
)
# SV tag info, use IP tag info as input
process.MySecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone(
  trackIPTagInfos = cms.InputTag("MyImpactParameterPFTagInfos"),
)


process.MyTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterPFTagInfos"))
)

process.MyCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterPFTagInfos"), 
			cms.InputTag("MySecondaryVertexTagInfos")			
	)
)
process.MyJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterPFTagInfos"))
)

process.myBTaggers = cms.Sequence(
  process.MyAk4PFJetTracksAssociatorAtVertex *
  process.MyImpactParameterPFTagInfos *
  (
     process.MyTrackCountingHighPurBJetTags +
     process.MyJetBProbabilityBJetTags + 
     (  process.MySecondaryVertexTagInfos *
        process.MyCombinedSecondaryVertexBJetTags
        )
     )
)





#btagging with calojets
process.MyAk4CaloJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
   process.j2tParametersVX,
   jets = cms.InputTag("ak4CaloJets")
)
process.MyImpactParameterCaloTagInfos = process.impactParameterTagInfos.clone(
  jetTracks = "MyAk4CaloJetTracksAssociatorAtVertex"
)
# SV tag info, use IP tag info as input
process.MySecondaryVertexCaloTagInfos = process.secondaryVertexTagInfos.clone(
  trackIPTagInfos = cms.InputTag("MyImpactParameterCaloTagInfos"),
)


process.MyTrackCountingHighPurCaloBJetTags = process.trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterCaloTagInfos"))
)

process.MyCombinedSecondaryVertexCaloBJetTags = process.combinedSecondaryVertexBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterCaloTagInfos"), 
			cms.InputTag("MySecondaryVertexCaloTagInfos")			
	)
)
process.MyJetBProbabilityCaloBJetTags = process.jetBProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterCaloTagInfos"))
)

process.myCaloBTaggers = cms.Sequence(
  process.MyAk4CaloJetTracksAssociatorAtVertex *
  process.MyImpactParameterCaloTagInfos *
  (
     process.MyTrackCountingHighPurCaloBJetTags +
     process.MyJetBProbabilityCaloBJetTags + 
     (  process.MySecondaryVertexCaloTagInfos *
        process.MyCombinedSecondaryVertexCaloBJetTags
        )
     )
)





# MC flavor identification
process.myPartons = cms.EDProducer("PartonSelector",
     withLeptons = cms.bool(False),
     src         = cms.InputTag("genParticles")
)


#for reco pf jet
process.flavourByRefPF = cms.EDProducer("JetPartonMatcher",
     jets = cms.InputTag("ak4PFJets"),
     coneSizeToAssociate = cms.double(0.3),
     partons = cms.InputTag("myPartons")
 )

process.flavourByValPF = cms.EDProducer("JetFlavourIdentifier",
     srcByReference = cms.InputTag("flavourByRefPF"),
     physicsDefinition = cms.bool(False)
)




#for reco Calo jet
process.flavourByRefCalo = cms.EDProducer("JetPartonMatcher",
     jets = cms.InputTag("ak4CaloJets"),
     coneSizeToAssociate = cms.double(0.3),
     partons = cms.InputTag("myPartons")
 )

process.flavourByValCalo = cms.EDProducer("JetFlavourIdentifier",
     srcByReference = cms.InputTag("flavourByRefCalo"),
     physicsDefinition = cms.bool(False)
)




#gen jet flavor identification
process.flavourByRefGenJet = cms.EDProducer("JetPartonMatcher",
     jets = cms.InputTag("ak4GenJets"),
     coneSizeToAssociate = cms.double(0.3),
     partons = cms.InputTag("myPartons")
 )
 
process.flavourByValGenJet = cms.EDProducer("JetFlavourIdentifier",
     srcByReference = cms.InputTag("flavourByRefGenJet"),
     physicsDefinition = cms.bool(False)
)



# data 
#SUPERCLUSTER_COLL_EB = "hybridSuperClusters"
#SUPERCLUSTER_COLL_EE = "multi5x5SuperClustersWithPreshower"

# mc
SUPERCLUSTER_COLL_EB = "correctedHybridSuperClusters"
SUPERCLUSTER_COLL_EE = "correctedMulti5x5SuperClustersWithPreshower"

#  SuperClusters  ################
process.superClusters = cms.EDProducer("SuperClusterMerger",
      src = cms.VInputTag(cms.InputTag( SUPERCLUSTER_COLL_EB ,"", "RECO"),
      cms.InputTag( SUPERCLUSTER_COLL_EE ,"", "RECO") )  
)



# quark-gluon likelihood separator
#NOTE74
# n/a
#process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')  
#process.QGTagger.srcJets = cms.InputTag('ak4PFJets')
#process.QGTagger.jec     = cms.untracked.string('ak4PFL1FastL2L3')




#NOTE74
#from RecoJets.JetProducers.pujetidsequence_cff import *
#process.recoPuJetId = puJetId.clone(
#   jets = cms.InputTag("ak4PFJets"),
#   applyJec = cms.bool(True),
#   inputIsCorrected = cms.bool(False),                
#)

#process.recoPuJetMva = puJetMva.clone(
#   jets = cms.InputTag("ak4PFJets"),
#   jetids = cms.InputTag("recoPuJetId"),
#   applyJec = cms.bool(True),
#   inputIsCorrected = cms.bool(False),                
#)




###################################################################
#
#   set up cut-based electron ID, using code from Ilya Kravchenko
#
#
###################################################################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)






#NOTE74
#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    calibratedElectrons = cms.PSet(
#        initialSeed = cms.untracked.uint32(1),
#        engineName = cms.untracked.string('TRandom3')
#    ),
#)

#process.load("EgammaAnalysis.ElectronTools.calibratedElectrons_cfi")

# dataset to correct
#process.calibratedElectrons.isMC = cms.bool(True)
#process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
#process.calibratedElectrons.updateEnergyError = cms.bool(True)
#process.calibratedElectrons.correctionsType = cms.int32(2)
#process.calibratedElectrons.combinationType = cms.int32(3)
#process.calibratedElectrons.lumiRatio = cms.double(0.607)
#process.calibratedElectrons.verbose = cms.bool(False)
#process.calibratedElectrons.synchronization = cms.bool(False)


# this was commnented out
#if (isRealData):
#  process.calibratedElectrons.isMC = cms.bool(False)
#  process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")
#else:
#  process.calibratedElectrons.isMC = cms.bool(True)
#  process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")



#process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
#process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
#process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
#process.eleRegressionEnergy.produceValueMaps = cms.bool(True)
#process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root")
#process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)




#process.load('EgammaAnalysis/ElectronTools/electronIdMVAProducer_cfi')





process.analyzer = cms.EDAnalyzer(
    "WZEdmAnalyzer",
    #parameters
    DEBUG                     = cms.bool(True),
    DATA                      = cms.bool(False),
    GEN_ONLY                  = cms.bool(False),
    MC_SIGNAL                 = cms.bool(True),
    VERTEXING                 = cms.bool(True),
    SMOOTHING                 = cms.bool(True),
    KVFParameters = cms.PSet(
    maxDistance = cms.double(0.01),
    maxNbrOfIterations = cms.int32(10)
    ),
    RECO                      = cms.bool(False),
    RECOSELECTION             = cms.string("DILEPTON"),
#choice among "BTAG", "DILEPTON", etc.
    BeamSpot                  = cms.InputTag("offlineBeamSpot"), 
    Vertices                  = cms.string(  "offlinePrimaryVertices"),
    Muons                     = cms.string(  "muons"),
    Electrons                 = cms.string(  "gedGsfElectrons"),
    EleLooseIdMap             = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    EleMediumIdMap            = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    EleTightIdMap             = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

    ElectronIsoVals	      = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
    Jets                      = cms.string("ak4CaloJets"),
    CaloJets                  = cms.string("ak4CaloJets"),
    JPTJets                   = cms.string("JetPlusTrackZSPCorJetAntiKt4"),
    PFJets                    = cms.string("ak4PFJets"),
    JetMinPt                  = cms.double(10),    
    LeptonThreshold           = cms.double(10),    
    InputJetIDValueMap         = cms.InputTag("ak4JetID"), 
    JetCorrectionService      = cms.string('ak4CaloL1FastL2L3'),
    CaloJetCorrectionService  = cms.string('ak4CaloL1FastL2L3'),
    JPTJetCorrectionService   = cms.string('ak4JPTL1FastL2L3'),
    PFJetCorrectionService    = cms.string('ak4PFL1FastL2L3'),
    RhoSrc                    = cms.InputTag('ak4PFJets', 'rho'),
    SigmaSrc                  = cms.InputTag('ak4PFJets', 'sigma'),
    RhoSrcCHS                 = cms.InputTag('ak4PFJetsCHS', 'rho'),
    SigmaSrcCHS               = cms.InputTag('ak4PFJetsCHS', 'sigma'),
    RhoSrcCalo                = cms.InputTag('ak4CaloJets', 'rho'),
    SigmaSrcCalo              = cms.InputTag('ak4CaloJets', 'sigma'),
    RhoSrcTrack               = cms.InputTag('ak4TrackJets', 'rho'),
    SigmaSrcTrack             = cms.InputTag('ak4TrackJets', 'sigma'),

    RhoIsoSrc                 = cms.InputTag('kt6PFJetsForIso', 'rho'),
    SigmaIsoSrc               = cms.InputTag('kt6PFJetsForIso', 'sigma'),
    RhoChSrc                  = cms.InputTag('kt6PFJetsForCh',  'rho'),
    SigmaChSrc                = cms.InputTag('kt6PFJetsForCh',  'sigma'),
    RhoCh2p4Src               = cms.InputTag('kt6PFJetsForCh2p4',  'rho'),
    SigmaCh2p4Src             = cms.InputTag('kt6PFJetsForCh2p4',  'sigma'), 
    Tracks                    = cms.string("generalTracks"),
    TrackMinPtWithMCTruth     = cms.double(10),
    LeptonMinPtForComposition = cms.double(15),
    JetTagCollection          = cms.string("trackCountingHighEffBJetTags"),
    PhotonCollection          = cms.string("photons"),
    L1ParticleMapCollection   = cms.string("l1extraParticleMap"),
    L1GTReadoutRecordLabel    = cms.InputTag('gtDigis'),
    HLTL1GTObjectMapLabel     = cms.InputTag('hltL1GtObjectMap'), 
    TriggerResultsLabel       = cms.InputTag("TriggerResults","","HLT"),
    TriggerSummaryLabel       = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    TriggerRefPath            = cms.string("HLTriggerFinalPath"), 
    HLTTriggerMuons           = cms.vstring("HLT_Mu40_vVERSION", "HLT_IsoMu24_v15", "HLT_Mu17_Mu8_v17"),
    HLTTriggerMuonElectrons   = cms.vstring("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", ""),
    HLTTriggerElectrons       = cms.vstring("", "HLT_Ele27_WP80_v10", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17"), 

    GeneratorLevelTag         = cms.string("generator"),
    LHEEventProductTag        = cms.InputTag("externalLHEProducer"),
    GenJets                   = cms.string("ak4GenJets"),
    akGenJets                 = cms.string("ak4GenJets"),
    GenJetMinPt               = cms.double(5),
    SimTracks                 = cms.string("g4SimHits"),
    DiLeptonMinMass           = cms.double(5),
    out                       = cms.string("./test_mc.root"),
    open                      = cms.string("recreate"),
    pdf                       = cms.string("cteq66.LHgrid"),
    subset                    = cms.int32(0)
    )




process.p = cms.Path(
#process.noscraping*
                     process.primaryVertexFilter*
#                     process.HBHENoiseFilter*
                     process.goodVertices * process.trackingFailureFilter *
                     process.trkPOGFilters*
                     process.myPartons*
		     process.flavourByRefPF*process.flavourByValPF*
		     process.flavourByRefCalo*process.flavourByValCalo*
                     process.flavourByRefGenJet*
                     process.flavourByValGenJet*
                     process.pfPileUpAllChargedParticlesClone*process.kt6PFJetsForCh*process.kt6PFJetsForCh2p4*
                     process.kt6PFJetsForIso*
                     process.ak4CaloJetsL1FastL2L3*
                     #process.pfiso*
#                     process.type0PFMEtCorrection*
#                     process.producePFMETCorrections*
#		     process.metMuonJESCorAK4*
#                     process.mymets*
                     process.superClusters*
                     process.myBTaggers*process.myCaloBTaggers*
 #                   process.QuarkGluonTagger*	
 #                    process.recoPuJetId * process.recoPuJetMva*
#                     process.eleRegressionEnergy * process.calibratedElectrons*
#                     process.mvaTrigV0  * process.mvaNonTrigV0*
                     process.egmGsfElectronIDSequence*
                     process.analyzer)






process.options = cms.untracked.PSet(
 wantSummary = cms.untracked.bool(True)
)



