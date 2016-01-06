import FWCore.ParameterSet.Config as cms
import string

#######################################################################
#
#  global control flag
#
#######################################################################
useAOD   = True
isData   = False # True data; False MC


if (isData == True) : 
    PFCHSJEC = "ak4PFCHSL1FastL2L3Residual"
    CaloJEC  = "ak4CaloL1FastL2L3Residual"
    JPTJEC   = "ak4JPTL1FastL2L3Residual"
    PFJEC    = "ak4PFL1FastL2L3Residual"
    HLT_SINGLEMU1 = "HLT_IsoTkMu20_vVERSION"
    HLT_SINGLEMU2 = "HLT_IsoMu20_vVERSION"
    HLT_DOUBLEMU1 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_vVERSION"
    HLT_DOUBLEMU2 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_vVERSION"
    HLT_MUE       = "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_vVERSION"
    HLT_EMU       = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_vVERSION"
    HLT_SINGLEELE1= "HLT_Ele23_WPLoose_Gsf_vVERSION"
    HLT_SINGLEELE2= "HLT_Ele23_WPLoose_Gsf_vVERSION"
    HLT_DOUBLEELE1= "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_vVERSION"
else :
    PFCHSJEC = "ak4PFCHSL1FastL2L3"
    CaloJEC  = "ak4CaloL1FastL2L3"
    JPTJEC   = "ak4JPTL1FastL2L3"
    PFJEC    = "ak4PFL1FastL2L3"
    HLT_SINGLEMU1 = "HLT_IsoTkMu20_v1"
    HLT_SINGLEMU2 = "HLT_IsoMu20_v1"
    HLT_DOUBLEMU1 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1"
    HLT_DOUBLEMU2 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"
    HLT_MUE       = "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"
    HLT_EMU       = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1"
    HLT_SINGLEELE1= "HLT_Ele27_WP85_Gsf_v1"
    HLT_SINGLEELE2= "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1"
    HLT_DOUBLEELE1= "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"


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

if isData == True :
    process.GlobalTag.globaltag = '74X_dataRun2_v5'
#    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'
else :
#    process.GlobalTag.globaltag = 'MCRUN2_74_V9'
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'


process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")




process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource", 
       	fileNames = cms.untracked.vstring(
       # '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root'
#'file:/uscms_data/d2/ptan/work/sl6/development/CMSSW_7_4_14/src/Analysis_RunII/WZEdmAnalyzer/test/002F7FDD-BA13-E511-AA63-0026189437F5.root'
'file:/uscms_data/d2/ptan/work/sl6/development/CMSSW_7_4_14/src/Analysis_RunII/WZEdmAnalyzer/plugins/00392203-FF0B-E511-AC4F-002590593872.root'
#'/store/data/Run2015D/SingleMuon/AOD/PromptReco-v4/000/258/159/00000/0C2C8F20-246C-E511-B27C-02163E0143D6.root'

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
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff', 
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff', 
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']

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

from RecoJets.JetPlusTracks.JetPlusTrackCorrections_cff import *
process.myJPTeidTight = process.JPTeidTight.clone()
process.myak4JetTracksAssociatorAtVertexJPT = process.ak4JetTracksAssociatorAtVertexJPT.clone()
process.myJetPlusTrackZSPCorJetAntiKt4 = process.JetPlusTrackZSPCorJetAntiKt4.clone()

process.myJetPlusTrackCorrectionsAntiKt4 = cms.Sequence(
     process.JPTeidTight*
     process.ak4JetTracksAssociatorAtVertexJPT*
     process.myJetPlusTrackZSPCorJetAntiKt4
     )


###################################################################
#
#    setup b-tagging MC truth, using hadron-based flavor identification
#
#
###################################################################
process.printList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(1)
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


#for jpt
process.flavourByRefJPT = cms.EDProducer("JetPartonMatcher",
                                          jets = cms.InputTag("myJetPlusTrackZSPCorJetAntiKt4"),
                                          coneSizeToAssociate = cms.double(0.3),
                                          partons = cms.InputTag("myPartons")
                                          )
process.flavourByValJPT = cms.EDProducer("JetFlavourIdentifier",
                                          srcByReference = cms.InputTag("flavourByRefJPT"),
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


from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.myak4PFJetFlavourInfos = ak4JetFlavourInfos.clone(jets = cms.InputTag("ak4PFJets"))

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.myak4PFJetCHSFlavourInfos = ak4JetFlavourInfos.clone(jets = cms.InputTag("ak4PFJetsCHS"))

#from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
#process.myak4CaloJetFlavourInfos = ak4JetFlavourInfos.clone(jets = cms.InputTag("ak4CaloJets"))

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.myak4GenJetFlavourInfos = ak4JetFlavourInfos.clone(jets = cms.InputTag("ak4GenJets"))






###################################################################
#
#   set up b-tagging sequence
#
#
###################################################################
#ak4pfCHS jet
process.MyPFCHSImpactParameterTagInfos = process.pfImpactParameterTagInfos.clone(
    jets = cms.InputTag("ak4PFJetsCHS") # use ak4PFJetsCHS stored in AOD as input
)
process.MyPFCHSSecondaryVertexTagInfos = process.pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("MyPFCHSImpactParameterTagInfos") # use the above IP TagInfos as input
)

process.MyPFCHSTrackCountingHighPurBJetTags = process.pfTrackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFCHSImpactParameterTagInfos"))
)

process.MyPFCHSCombinedSecondaryVertexV2BJetTags = process.pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFCHSImpactParameterTagInfos"), 
                           cms.InputTag("MyPFCHSSecondaryVertexTagInfos")			
	)
)
process.MyPFCHSJetProbabilityBJetTags = process.pfJetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFCHSImpactParameterTagInfos"))
)

process.myPFCHSBTaggers = cms.Sequence(
  process.MyPFCHSImpactParameterTagInfos *
  (
     process.MyPFCHSTrackCountingHighPurBJetTags +
     process.MyPFCHSJetProbabilityBJetTags + 
     (  process.MyPFCHSSecondaryVertexTagInfos *
        process.MyPFCHSCombinedSecondaryVertexV2BJetTags
        )
     )
)



#pf jets
process.MyPFImpactParameterTagInfos = process.pfImpactParameterTagInfos.clone(
    jets = cms.InputTag("ak4PFJets") # use ak4PFJetsCHS stored in AOD as input
)
process.MyPFSecondaryVertexTagInfos = process.pfInclusiveSecondaryVertexFinderTagInfos.clone(
    trackIPTagInfos = cms.InputTag("MyPFImpactParameterTagInfos") # use the above IP TagInfos as input
)

process.MyPFTrackCountingHighPurBJetTags = process.pfTrackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFImpactParameterTagInfos"))
)

process.MyPFCombinedSecondaryVertexV2BJetTags = process.pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFImpactParameterTagInfos"), 
                           cms.InputTag("MyPFSecondaryVertexTagInfos")			
	)
)
process.MyPFJetProbabilityBJetTags = process.pfJetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyPFImpactParameterTagInfos"))
)

process.myPFBTaggers = cms.Sequence(
  process.MyPFImpactParameterTagInfos *
  (
     process.MyPFTrackCountingHighPurBJetTags +
     process.MyPFJetProbabilityBJetTags + 
     (  process.MyPFSecondaryVertexTagInfos *
        process.MyPFCombinedSecondaryVertexV2BJetTags
        )
     )
)



#btagging with calojets
process.MyAk4CaloJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
   process.j2tParametersVX,
   jets = cms.InputTag("ak4CaloJets")
)
process.MyCaloImpactParameterTagInfos = process.impactParameterTagInfos.clone(
  jetTracks = "MyAk4CaloJetTracksAssociatorAtVertex"
)
# SV tag info, use IP tag info as input
process.MyCaloSecondaryVertexTagInfos = process.inclusiveSecondaryVertexFinderTagInfos.clone(
  trackIPTagInfos = cms.InputTag("MyCaloImpactParameterTagInfos"),
)


process.MyCaloTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyCaloImpactParameterTagInfos"))
)

process.MyCaloCombinedSecondaryVertexV2BJetTags = process.combinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("MyCaloImpactParameterTagInfos"), 
                             cms.InputTag("MyCaloSecondaryVertexTagInfos")			
                             )
)
process.MyCaloJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyCaloImpactParameterTagInfos"))
)

process.myCaloBTaggers = cms.Sequence(
  process.MyAk4CaloJetTracksAssociatorAtVertex *
  process.MyCaloImpactParameterTagInfos *
  (
     process.MyCaloTrackCountingHighPurBJetTags +
     process.MyCaloJetProbabilityBJetTags + 
     (  process.MyCaloSecondaryVertexTagInfos *
        process.MyCaloCombinedSecondaryVertexV2BJetTags
        )
     )
)







#btagging with JPTjets
process.MyAk4JPTJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
   process.j2tParametersVX,
   jets = cms.InputTag("myJetPlusTrackZSPCorJetAntiKt4")
)
process.MyJPTImpactParameterTagInfos = process.impactParameterTagInfos.clone(
  jetTracks = "MyAk4JPTJetTracksAssociatorAtVertex"
)
# SV tag info, use IP tag info as input
process.MyJPTSecondaryVertexTagInfos = process.inclusiveSecondaryVertexFinderTagInfos.clone(
  trackIPTagInfos = cms.InputTag("MyJPTImpactParameterTagInfos"),
)


process.MyJPTTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyJPTImpactParameterTagInfos"))
)

process.MyJPTCombinedSecondaryVertexV2BJetTags = process.combinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("MyJPTImpactParameterTagInfos"), 
                             cms.InputTag("MyJPTSecondaryVertexTagInfos")			
                             )
)
process.MyJPTJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("MyJPTImpactParameterTagInfos"))
)

process.myJPTBTaggers = cms.Sequence(
  process.MyAk4JPTJetTracksAssociatorAtVertex *
  process.MyJPTImpactParameterTagInfos *
  (
     process.MyJPTTrackCountingHighPurBJetTags +
     process.MyJPTJetProbabilityBJetTags + 
     (  process.MyJPTSecondaryVertexTagInfos *
        process.MyJPTCombinedSecondaryVertexV2BJetTags
        )
     )
)



process.analyzer = cms.EDAnalyzer(
    "WZEdmAnalyzer",
    #parameters
    DEBUG                     = cms.bool(False),
    DATA                      = cms.bool( isData ),
    GEN_ONLY                  = cms.bool(False),
    SAVE_ALLEVENTS            = cms.bool(True),
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
    EffAreasConfigFile        = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
    EleLooseIdMap             = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    EleMediumIdMap            = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    EleTightIdMap             = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    ElectronEcalPFClusterIsolationProducer = cms.InputTag("electronEcalPFClusterIsolationProducer"),
    ElectronHcalPFClusterIsolationProducer = cms.InputTag("electronHcalPFClusterIsolationProducer"),
    TrigMvaValuesMap          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
    TrigMvaCategoriesMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
    TrigMvaMediumIdMaps       = cms.InputTag( "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"), 
    TrigMvaTightIdMaps        = cms.InputTag( "egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"), 
    NonTrigMvaValuesMap       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    NonTrigMvaCategoriesMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
    NonTrigMvaMediumIdMaps    = cms.InputTag( "egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"), 
    NonTrigMvaTightIdMaps     = cms.InputTag( "egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"), 
    PFCHSJets                 = cms.string(  "ak4PFJetsCHS"),
    PFCHSJetFlavourInfos      = cms.InputTag("myak4PFJetCHSFlavourInfos"),
    PFCHSJetTagInfos          = cms.vstring( "MyPFCHSTrackCountingHighPurBJetTags", "MyPFCHSJetProbabilityBJetTags", "MyPFCHSCombinedSecondaryVertexV2BJetTags"),
    CaloJets                  = cms.string(  "ak4CaloJets"),
    CaloJetFlavourInfos       = cms.InputTag("flavourByValCalo"), 
    CaloJetTagInfos           = cms.vstring( "MyCaloTrackCountingHighPurBJetTags",  "MyCaloJetProbabilityBJetTags", "MyCaloCombinedSecondaryVertexV2BJetTags"),
    JPTJets                   = cms.string(  "myJetPlusTrackZSPCorJetAntiKt4"),
    JPTJetFlavourInfos        = cms.InputTag("flavourByValJPT"), 
    JPTJetTagInfos            = cms.vstring( "MyJPTTrackCountingHighPurBJetTags",   "MyJPTJetProbabilityBJetTags", "MyJPTCombinedSecondaryVertexV2BJetTags"),
    PFJets                    = cms.string(  "ak4PFJets"),
    PFJetFlavourInfos         = cms.InputTag("myak4PFJetFlavourInfos"),
    PFJetTagInfos             = cms.vstring( "MyPFTrackCountingHighPurBJetTags",    "MyPFJetProbabilityBJetTags", "MyPFCombinedSecondaryVertexV2BJetTags"),
    JetTagCollections         = cms.vstring( "trackCountingHighEffBJetTags", "pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    JetMinPt                  = cms.double(10),    
    LeptonThreshold           = cms.double(10),    
    InputJetIDValueMap        = cms.InputTag("ak4JetID"), 
    PFCHSJetCorrectionService = cms.string(  PFCHSJEC ), 
    #'ak4PFCHSL1FastL2L3'),
    CaloJetCorrectionService  = cms.string(  CaloJEC), 
    #'ak4CaloL1FastL2L3'),
    JPTJetCorrectionService   = cms.string(  JPTJEC), 
    #'ak4JPTL1FastL2L3'),
    PFJetCorrectionService    = cms.string(  PFJEC), 
    #'ak4PFL1FastL2L3'),
    FixGridRho                = cms.InputTag('fixedGridRhoFastjetAll'),
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
    PhotonCollection          = cms.string("gedPhotons"),
    L1ParticleMapCollection   = cms.string("l1extraParticleMap"),
    L1GTReadoutRecordLabel    = cms.InputTag('gtDigis'),
    HLTL1GTObjectMapLabel     = cms.InputTag('hltL1GtObjectMap'), 
    TriggerResultsLabel       = cms.InputTag("TriggerResults","","HLT"),
    TriggerSummaryLabel       = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    TriggerRefPath            = cms.string("HLTriggerFinalPath"), 
    #HLTTriggerMuons           = cms.vstring("HLT_IsoTkMu20_v1", "HLT_IsoMu20_v1", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"),
    HLTTriggerMuons           = cms.vstring(HLT_SINGLEMU1, HLT_SINGLEMU2, HLT_DOUBLEMU1, HLT_DOUBLEMU2),
    #HLTTriggerMuonElectrons   = cms.vstring("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1"),
    HLTTriggerMuonElectrons   = cms.vstring(HLT_MUE, HLT_EMU),
    #HLTTriggerElectrons       = cms.vstring("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1", "HLT_Ele27_WP85_Gsf_v1", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"), 
    HLTTriggerElectrons       = cms.vstring(HLT_SINGLEELE1, HLT_SINGLEELE2, HLT_DOUBLEELE1), 
    GeneratorLevelTag         = cms.string("generator"),
    LHEEventProductTag        = cms.InputTag("externalLHEProducer"),
    GenJets                   = cms.string(  "ak4GenJets"),
    akGenJets                 = cms.string(  "ak4GenJets"),
    akGenJetFlavourInfos      = cms.InputTag("myak4GenJetFlavourInfos"),
    GenJetMinPt               = cms.double(5),
    SimTracks                 = cms.string("g4SimHits"),
    DiLeptonMinMass           = cms.double(5)
    )


# rename output file
process.TFileService = cms.Service("TFileService",
                                   fileName = 
                                   cms.string('testme.root' ),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



if isData == True :

    process.p = cms.Path(
        #process.noscraping*
        process.primaryVertexFilter*
        #                     process.HBHENoiseFilter*
        process.goodVertices * process.trackingFailureFilter *
        process.trkPOGFilters*
        process.pfPileUpAllChargedParticlesClone*process.kt6PFJetsForCh*process.kt6PFJetsForCh2p4*
        process.kt6PFJetsForIso*
        process.ak4CaloJetsL1FastL2L3Residual*
        #process.pfiso*
        #                    process.type0PFMEtCorrection*
        #                    process.producePFMETCorrections*
        #		     process.metMuonJESCorAK4*
        #                    process.mymets*
        process.superClusters*
        process.myJetPlusTrackCorrectionsAntiKt4*
        process.myPFCHSBTaggers*process.myPFBTaggers*process.myCaloBTaggers*process.myJPTBTaggers*
        #                    process.QuarkGluonTagger*	
        #                    process.recoPuJetId * process.recoPuJetMva*
        #                    process.eleRegressionEnergy * process.calibratedElectrons*
        #                    process.mvaTrigV0  * process.mvaNonTrigV0*
        process.egmGsfElectronIDSequence*
        process.analyzer)
    
else : 

    process.p = cms.Path(
        #process.noscraping*
        process.printList*
        process.primaryVertexFilter*
        #                     process.HBHENoiseFilter*
        process.goodVertices * process.trackingFailureFilter *
        process.trkPOGFilters*
        process.selectedHadronsAndPartons*
        process.myak4PFJetFlavourInfos*process.myak4PFJetCHSFlavourInfos*process.myak4GenJetFlavourInfos*
        process.myPartons*
        process.flavourByRefPF*process.flavourByValPF*
        process.flavourByRefCalo*process.flavourByValCalo*
        process.flavourByRefGenJet*process.flavourByValGenJet*
        process.pfPileUpAllChargedParticlesClone*process.kt6PFJetsForCh*process.kt6PFJetsForCh2p4*
        process.kt6PFJetsForIso*
        process.ak4CaloJetsL1FastL2L3*
        #process.pfiso*
        #                    process.type0PFMEtCorrection*
        #                    process.producePFMETCorrections*
        #		     process.metMuonJESCorAK4*
        #                    process.mymets*
        process.superClusters*
        process.myJetPlusTrackCorrectionsAntiKt4*
        process.flavourByRefJPT*process.flavourByValJPT*
        process.myPFCHSBTaggers*process.myPFBTaggers*process.myCaloBTaggers*process.myJPTBTaggers*
        #                    process.QuarkGluonTagger*	
        #                    process.recoPuJetId * process.recoPuJetMva*
        #                    process.eleRegressionEnergy * process.calibratedElectrons*
        #                    process.mvaTrigV0  * process.mvaNonTrigV0*
        process.egmGsfElectronIDSequence*
        #process.pfWeightedIsoSeq *
        process.analyzer)
    




process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )



