import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaElectronProducers.gedGsfElectronFinalizer_cfi import gedGsfElectrons
gedGsfElectronsClone = gedGsfElectrons.clone()


from EgammaAnalysis.ElectronTools.regressionModifier_cfi import regressionModifier
regressionModifier.ecalrechitsEB = cms.InputTag("reducedEgamma:reducedEBRecHits")
regressionModifier.ecalrechitsEE = cms.InputTag("reducedEgamma:reducedEERecHits")
regressionModifier.useLocalFile  = cms.bool(False)

egamma_modifications = cms.VPSet( )
egamma_modifications.append( regressionModifier )

gedGsfElectronsClone.regressionConfig = regressionModifier

#slimmedElectrons.modifierConfig.modifications = egamma_modifications
#slimmedPhotons.modifierConfig.modifications   = egamma_modifications

regressionApplication = cms.Sequence( gedGsfElectronsClone  )
