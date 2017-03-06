import FWCore.ParameterSet.Config as cms

from EgammaAnalysis.ElectronTools.regressionModifier_cfi import regressionModifier

regressionModifier.ecalrechitsEB = cms.InputTag("reducedEcalRecHitsEB")
regressionModifier.ecalrechitsEE = cms.InputTag("reducedEcalRecHitsEE")

egamma_modifications = cms.VPSet( )
egamma_modifications.append( regressionModifier )

gedGsfElectronsClone = cms.EDProducer(
    "ModifiedGsfElectronProducer",
    src = cms.InputTag("gedGsfElectrons",processName=cms.InputTag.skipCurrentProcess()),
    modifierConfig = cms.PSet( modifications = egamma_modifications )
)

gedPhotonsClone = cms.EDProducer(
    "ModifiedGedPhotonProducer",
    src = cms.InputTag("gedPhotons",processName=cms.InputTag.skipCurrentProcess()),
    modifierConfig = cms.PSet( modifications = egamma_modifications )
)

regressionApplication = cms.Sequence( gedGsfElectronsClone * gedPhotonsClone )
