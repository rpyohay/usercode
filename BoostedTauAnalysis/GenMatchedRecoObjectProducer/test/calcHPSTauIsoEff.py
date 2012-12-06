import FWCore.ParameterSet.Config as cms

process = cms.Process("EFFANALYSIS")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
W_PDGID = 24
TAU_PDGID = 15
MU_PDGID = 13
NUMU_PDGID = 14
D_PDGID = 1
U_PDGID = 2
S_PDGID = 3
C_PDGID = 4
B_PDGID = 5
T_PDGID = 6
G_PDGID = 21
ANY_PDGID = 0

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

#tau hadronic decay types
TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

#no consideration of pT rank
ANY_PT_RANK = -1

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    INPUT_FILES
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#for mu-less jets
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('BoostedTauAnalysis/CleanJets/cleanjets_cfi')

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(MOM_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

#define a parameter set for the W-->munu selector
WMuNuPSet = commonGenTauDecayIDPSet.clone()
WMuNuPSet.momPDGID = cms.vint32(W_PDGID)

#only proceed if event is a true W-->munu event
process.genWMuNuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(NUMU_PDGID),
    genTauDecayIDPSet = WMuNuPSet,
    primaryTauDecayType = cms.uint32(TAU_ALL),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#require event to fire IsoMu24_eta2p1
process.IsoMu24eta2p1Selector = process.hltHighLevel.clone()
process.IsoMu24eta2p1Selector.HLTPaths = cms.vstring('HLT_IsoMu24_eta2p1_v11')

#search for a muon with pT > 20 GeV as in WHbb CMS AN-2012/349 and proceed if one can be found
#this will produce a ref to the original muon collection
process.WMuonPTSelector = cms.EDFilter('MuonRefSelector',
                                       src = cms.InputTag('muons'),
                                       cut = cms.string('pt > 25.0'),
                                       filter = cms.bool(True)
                                       )

#search for a loose PF isolated tight muon in |eta| < 2.4 with pT > 20 GeV
#(see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1 for
#isolation definition;CMS AN-2012/349 uses loose isolation working point for WHbb muon selection)
#this will produce a ref to the original muon collection
process.WIsoMuonSelector = cms.EDFilter('CustomMuonSelector',
                                        baseMuonTag = cms.InputTag('muons'),
                                        muonTag = cms.InputTag('WMuonPTSelector'),
                                        vtxTag = cms.InputTag('offlinePrimaryVertices'),
                                        muonID = cms.string('tight'),
                                        PFIsoMax = cms.double(0.2),
                                        detectorIsoMax = cms.double(-1.0),
                                        PUSubtractionCoeff = cms.double(0.5),
                                        usePFIso = cms.bool(True),
                                        passIso = cms.bool(True),
                                        etaMax = cms.double(2.4),
                                        minNumObjsToPassFilter = cms.uint32(1)
                                        )

#produce gen tau mu collection
process.genTauMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    primaryTauDecayType = cms.uint32(TAU_MU),
    sisterTauDecayType = cms.uint32(TAU_HAD),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#produce muons
process.recoMuonSelector = cms.EDFilter('MuonRefSelector',
                                        src = cms.InputTag('muons'),
                                        cut = cms.string('pt >= 0.0'),
                                        filter = cms.bool(True)
                                        )

#produce gen-tau-mu-matched muons
process.genTauMuMatchedRecoMuonSelector = cms.EDFilter(
    'GenMatchedMuonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauMuSelector'),
    recoObjTag = cms.InputTag('recoMuonSelector'),
    baseRecoObjTag = cms.InputTag('muons'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(1)
    )

#produce gen-tau-mu-matched muons passing soft ID
process.genTauMuMatchedRecoSoftMuonSelector = cms.EDFilter(
    'CustomMuonSelector',
    baseMuonTag = cms.InputTag('muons'),
    muonTag = cms.InputTag('genTauMuMatchedRecoMuonSelector'),
    vtxTag = cms.InputTag('offlinePrimaryVertices'),
    vetoMuonTag = cms.InputTag('WIsoMuonSelector'),
    muonID = cms.string('soft'),
    PFIsoMax = cms.double(0.2),
    detectorIsoMax = cms.double(-1.0),
    PUSubtractionCoeff = cms.double(0.5),
    usePFIso = cms.bool(True),
    passIso = cms.bool(True),
    etaMax = cms.double(-1.0),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#clean the jets of soft muons, then rebuild the taus
process.CleanJets.muonSrc = cms.InputTag('genTauMuMatchedRecoSoftMuonSelector')
process.CleanJets.cutOnGenMatches = cms.bool(False)
process.CleanJets.outFileName = cms.string('NMSSMSignal_MuProperties.root')
process.recoTauAK5PFJets08Region.src = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.ak5PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJets", "ak5PFJetsNoMu",
                                                             "EFFANALYSIS")
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu",
                                                        "EFFANALYSIS")
process.recoTauCommonSequence = cms.Sequence(process.CleanJets*
                                             process.ak5PFJetTracksAssociatorAtVertex*
                                             process.recoTauAK5PFJets08Region*
                                             process.recoTauPileUpVertices*
                                             process.pfRecoTauTagInfoProducer
                                             )
process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence)

#produce gen-tau-mu-matched muons passing soft ID that are in jets
#will get an output map between new jet and removed muon refs
#see how often the input soft muons are in the output map

#compute efficiency of soft ID for gen-tau-mu-matched muons
process.softMuEffAnalyzer = cms.EDAnalyzer(
    'MuonEfficiencyAnalyzer',
    outFileName = cms.string('/data1/yohay/Wh1_Medium/softMuEff_v1.root'),
    denominatorTag = cms.InputTag('genTauMuMatchedRecoMuonSelector'),
    numeratorTag = cms.InputTag('genTauMuMatchedRecoSoftMuonSelector'),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    decayModeColors = cms.vuint32(1, 2, 4, 6, 8),
    decayModeStyles = cms.vuint32(20, 21, 22, 23, 24),
    decayModeEntries = cms.vstring('#tau_{#mu}', '#tau_{had}, 1 prong',
                                   '#tau_{had}, 1 prong + 1 #pi^{0}',
                                   '#tau_{had}, 1 prong + 2 #pi^{0}', '#tau_{had}, 3 prong')
    )

#produce HPS taus in |eta| < 2.4
process.recoTauSelector = cms.EDFilter(
    'CustomTauSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'ANALYSIS'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'ANALYSIS')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'ANALYSIS'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#produce gen-tau-mu-matched HPS taus in |eta| < 2.4
process.genTauMuMatchedRecoTauSelector = cms.EDFilter(
    'GenMatchedTauProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauSelector'),
    recoObjTag = cms.InputTag('recoTauSelector'),
    baseRecoObjTag = cms.InputTag('hpsPFTauProducer', '', 'ANALYSIS'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(0)
    )

#produce gen-matched HPS taus in |eta| < 2.4 passing MediumCombinedIsolationDBSumPtCorr
process.genMatchedIsoRecoTauSelector = process.recoTauSelector.clone()
process.genMatchedIsoRecoTauSelector.tauTag = cms.InputTag('genMatchedRecoTauSelector')
process.genMatchedIsoRecoTauSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'ANALYSIS')
    )

#compute efficiency of MediumCombinedIsolationDBSumPtCorr for gen-matched HPS taus in |eta| < 2.4
process.isoEffGenMatchedRecoTauAnalyzer = cms.EDAnalyzer(
    'ANALYZER',
    outFileName = cms.string(
    '/data1/yohay/SAMPLE_MCTruthMuonRemoval_MediumCombinedIsolationDBSumPtCorr.root'
    ),
    TAG1 = cms.InputTag('genMatchedRecoTauSelector'),
    TAG2 = cms.InputTag('genMatchedIsoRecoTauSelector'),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    decayModeColors = cms.vuint32(1, 2, 4, 6, 8),
    decayModeStyles = cms.vuint32(20, 21, 22, 23, 24),
    decayModeEntries = cms.vstring('#tau_{#mu}', '#tau_{had}, 1 prong',
                                   '#tau_{had}, 1 prong + 1 #pi^{0}',
                                   '#tau_{had}, 1 prong + 2 #pi^{0}', '#tau_{had}, 3 prong')
    )

#sequences
process.softMuEffSeq = cms.Sequence(process.genWMuNuSelector*process.IsoMu24eta2p1Selector*
                                    process.WMuonPTSelector*process.WIsoMuonSelector*
                                    process.genTauMuSelector*process.recoMuonSelector*
                                    process.genTauMuMatchedRecoMuonSelector*
                                    process.genTauMuMatchedRecoSoftMuonSelector*## process.PFTau*
                                    process.softMuEffAnalyzer)

#path
process.p = cms.Path(process.softMuEffSeq)
