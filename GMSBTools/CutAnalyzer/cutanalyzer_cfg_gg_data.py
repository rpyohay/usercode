#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 35X/36X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV8
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.29 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
    annotation = cms.untracked.string('SUSY pattuple definition')
)

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_100_1_bMI.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_101_1_C2I.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_102_1_HFu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_103_1_prt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_104_1_jJX.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_105_1_XHL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_106_1_biX.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_107_1_Ogv.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_108_1_npS.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_109_1_4hs.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_10_1_BBE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_110_1_7m1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_111_1_1rK.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_112_1_VdR.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_113_1_CZf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_114_1_7NR.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_115_1_U4x.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_116_1_7ri.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_117_1_3lG.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_118_1_1fE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_119_1_gt0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_11_1_guI.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_120_1_uAU.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_121_1_gDi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_122_1_9B0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_123_1_H1u.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_124_1_pGY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_125_1_gQ7.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_126_1_N3k.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_127_1_0d0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_128_1_exe.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_129_1_7V8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_12_1_3f5.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_130_1_kJf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_131_1_hID.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_132_1_9Ua.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_133_1_2Bm.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_134_1_AdI.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_135_1_ZJk.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_136_1_R8x.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_137_1_AiB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_138_1_syq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_139_1_sgy.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_13_1_6MM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_140_1_Ms6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_141_1_R0o.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_142_1_Keb.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_143_1_3ER.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_144_1_a5I.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_145_1_iXw.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_146_1_1Vf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_147_1_UrO.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_148_1_Qle.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_149_1_RSs.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_14_1_Mjy.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_150_1_ED6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_151_1_ixP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_152_1_NA3.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_153_1_X6f.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_154_1_Smd.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_155_1_Tib.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_156_1_Xpp.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_157_1_aFs.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_158_1_65d.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_159_1_ube.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_15_1_pGP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_160_1_Nf3.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_161_1_xek.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_162_1_pqY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_163_1_PTg.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_164_1_1gn.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_165_1_0HY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_166_1_5Se.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_167_1_pgW.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_168_1_LY8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_169_1_Hp8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_16_1_eFo.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_170_1_aBE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_171_1_U2L.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_172_1_MaP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_173_1_YKS.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_174_1_FSu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_175_1_yoM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_176_1_C9p.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_177_1_d6N.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_178_1_6M3.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_179_1_wBE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_17_1_sro.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_180_1_TSQ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_181_1_Nk1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_18_1_cii.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_19_1_zRF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_1_1_YWN.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_20_1_VZD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_21_1_34F.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_22_1_gEv.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_23_1_24a.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_24_1_1yv.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_25_1_MXT.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_26_1_NQK.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_27_1_vPN.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_28_1_sQ4.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_29_1_Ej4.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_2_1_v9N.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_30_1_JkE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_31_2_Wdw.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_32_1_QuZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_33_1_peV.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_34_2_Hh6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_35_1_PL9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_36_1_QSU.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_37_1_Kh1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_38_1_0Kz.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_39_1_MeQ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_3_1_QEC.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_40_1_J22.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_41_1_tdm.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_42_1_FLy.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_43_1_1Yc.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_44_1_m9Y.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_45_1_ScK.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_46_1_aB0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_47_1_Qqy.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_48_1_zSi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_49_1_IIH.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_4_1_eCe.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_50_1_yW1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_51_1_olG.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_52_1_7LJ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_53_1_x1s.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_54_1_wrT.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_55_1_baf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_56_1_lgL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_57_1_8DD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_58_1_69E.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_59_1_d61.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_5_1_rvO.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_60_1_nK8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_61_1_gNM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_62_1_fYm.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_63_1_N8g.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_64_1_see.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_65_1_r97.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_66_1_Sz2.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_67_1_4Kp.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_68_1_xRM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_69_1_7Es.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_6_1_VQt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_70_1_6JZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_71_1_dxU.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_72_1_B9m.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_73_1_LAt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_74_1_aex.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_75_1_zWi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_76_1_BHu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_77_1_Q6S.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_78_1_UHc.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_79_1_lFJ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_7_1_6Vo.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_80_1_Bfc.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_81_1_tad.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_82_1_0yW.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_83_1_6uP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_84_1_eXG.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_85_1_kPZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_86_1_GaF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_87_1_W8j.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_88_1_Fc7.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_89_1_Jwu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_8_1_3xW.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_90_1_zzT.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_91_1_pL2.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_92_1_DPM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_93_1_KIt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_94_1_R1v.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_95_1_YF1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_96_1_TcV.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_97_1_XYN.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_98_1_xYS.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_99_1_88g.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/gg_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_9_1_mPC.root'
    ]
process.maxEvents.input = -1
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'GR_R_36X_V11A::All'

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, mcVersion ('35x' for 35x samples, empty string for 36X samples),JetCollections
addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5JPT']) 
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################

#-- Output module configuration -----------------------------------------------
process.out.fileName = '/data/yohay/SUSYPAT_gg_data.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

#to avoid running out of memory
process.options = cms.untracked.PSet(
    fileMode = cms.untracked.string('NOMERGE')
    )

#load the sample analyzer
process.load("GMSBTools.CutAnalyzer.cutanalyzer_cfi_gg_data")

#load the calo geometry
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

#show TrigReport
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.susyPatDefaultSequence*process.CutAnalyzer )
