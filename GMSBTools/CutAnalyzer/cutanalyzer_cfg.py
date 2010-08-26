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

#keep the logging output to a nice level
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(1) )

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_100_1_OX6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_101_1_bXh.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_102_1_3pR.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_103_1_YtJ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_104_1_a9c.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_105_1_xqL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_106_1_5hw.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_107_1_ZUq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_108_1_5N9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_109_1_R4v.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_10_1_dEU.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_110_1_0tX.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_111_1_t9H.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_112_1_cOi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_113_1_7LY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_114_1_rSj.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_115_1_UaL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_116_1_4ne.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_117_1_l2L.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_118_1_Gzd.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_119_1_Wnt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_11_1_G4k.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_120_1_b7I.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_121_1_8oD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_122_1_Ake.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_123_1_DyS.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_124_1_eNQ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_125_1_nQx.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_126_1_iyS.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_127_1_hFY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_128_1_Os0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_129_1_9VC.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_12_3_iMe.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_130_1_5w6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_131_1_mFf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_132_1_u3n.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_133_1_7UD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_134_1_qSX.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_135_1_34I.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_136_1_rTR.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_137_1_BBW.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_138_1_ObE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_139_1_BmP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_13_1_vs9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_140_1_3S7.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_141_1_Cwm.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_142_1_Scq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_143_1_l0m.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_144_1_oHG.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_145_1_oPo.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_146_1_FkB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_147_1_hOK.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_148_1_yif.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_149_1_8dt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_14_1_bz9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_150_1_IZp.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_151_1_nzV.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_152_1_Yyx.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_153_1_M6X.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_154_1_RKp.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_155_1_gql.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_156_1_bpt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_157_1_aii.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_158_1_jqf.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_159_1_EzE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_15_1_WFD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_160_1_V1v.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_161_1_LWx.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_162_1_h3N.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_163_1_tPz.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_164_1_K3e.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_165_1_Qr7.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_166_1_tux.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_167_1_GSk.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_168_1_Qud.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_169_1_SGD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_16_1_3pm.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_170_1_h5q.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_171_1_gkK.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_172_1_NyF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_173_1_ff0.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_174_1_fKj.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_175_1_rg5.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_176_1_JCR.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_177_1_p7t.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_178_1_5cC.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_179_1_vpi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_17_1_FBt.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_180_1_0x1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_181_1_FGA.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_18_1_WTq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_19_1_4jZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_1_1_zSH.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_20_1_zd2.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_21_1_Cly.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_22_1_Af2.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_23_1_3y4.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_24_1_45B.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_25_1_wQG.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_26_1_AZq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_27_1_6Tp.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_28_1_dre.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_29_1_SYx.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_2_1_Qxk.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_30_1_azi.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_31_1_rFu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_32_1_cHL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_33_1_1r8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_34_1_XL8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_35_1_GXh.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_36_1_DM8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_37_1_kEA.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_38_1_GQP.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_39_1_smz.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_3_1_h1g.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_40_1_RGZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_41_1_yy5.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_42_1_GpB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_43_1_Ita.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_44_1_jfZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_45_1_H9b.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_46_1_R2j.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_47_1_dkD.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_48_1_3vv.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_49_1_79R.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_4_1_0in.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_50_1_7ok.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_51_1_D9z.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_52_1_riM.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_53_1_mlu.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_54_1_ZbE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_55_1_mN9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_56_1_J6s.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_57_1_4lF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_58_1_RuZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_59_1_7b4.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_5_1_cBq.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_60_1_Xoh.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_61_1_cm7.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_62_1_oDC.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_63_1_o6I.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_64_1_gnY.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_65_1_zr5.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_66_1_h4n.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_67_1_ksX.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_68_1_jgv.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_69_1_JlB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_6_2_Tu1.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_6_3_EIU.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_70_1_1zO.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_71_1_40Y.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_72_1_wQh.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_73_1_eDF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_74_1_h2j.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_75_1_1LB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_76_1_KN5.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_77_1_eJC.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_78_1_mcz.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_79_1_8T8.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_7_1_By9.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_80_1_OPO.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_81_1_MUb.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_82_1_wGE.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_83_1_QuQ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_84_1_ut6.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_85_1_mam.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_86_1_LKo.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_87_1_mMA.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_88_1_y8r.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_89_1_w2a.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_8_1_O0A.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_90_1_0Tw.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_91_1_Cj3.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_92_1_Tcl.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_93_1_prL.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_94_1_hzZ.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_95_1_tJB.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_96_1_TRF.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_97_1_Ciz.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_98_1_0La.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_99_1_ADb.root',
    'rfio:/castor/cern.ch/user/y/yohay/361p4/ff_361p4_EG_PD_prompt_reco_runs138560-139790_HLT15_withTimingCut_9_1_eKP.root'
    ]
#process.source.skipEvents = cms.untracked.uint32(49)
process.maxEvents.input = -1
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #)
#)

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
process.out.fileName = '/data/yohay/SUSYPAT_jet_study_RECO.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

#process.out = cms.OutputModule("PoolOutputModule",
                               #fileName = cms.untracked.string("ff_361p4_EG_SD_re-reco_runs132440-135735_HLT15_withoutTimingCut.root"),
                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring( 'p' ) )
                               #)

#to avoid running out of memory
process.options = cms.untracked.PSet(
    fileMode = cms.untracked.string('NOMERGE')
    )
#process.source.inputCommands = cms.untracked.vstring("keep *",  
#                                                     "drop *_MEtoEDMConverter_*_*")

#load the HLT filter
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT")
process.hltHighLevel.HLTPaths = cms.vstring("HLT_Photon15_L1R")

#load the sample maker
process.load("GMSBTools.SampleMaker.samplemaker_cfi")

#load the sample analyzer
process.load("GMSBTools.CutAnalyzer.cutanalyzer_cfi")

#load the GMSB analyzer
#process.load("GMSBTools.GMSBAnalyzer.gmsbanalyzer_cfi")

#load the calo geometry
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_36X_V11A::All' #not the latest global tag, but the one used by Michael to create the skim
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

#show TrigReport
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.susyPatDefaultSequence*process.CutAnalyzer )
#process.p = cms.Path( process.CutAnalyzer )
#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()

#make the sample and check the distributions
#process.p = cms.Path(process.hltHighLevel*process.SampleMaker*process.CutAnalyzer*process.GMSBAnalyzer)
#process.p = cms.Path(process.hltHighLevel*process.SampleMaker*process.CutAnalyzer)
#process.p = cms.Path(process.SampleMaker*process.CutAnalyzer)
#process.p = cms.Path(process.SampleMaker*process.GMSBAnalyzer)
#process.e = cms.EndPath(process.out)
