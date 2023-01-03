# example pipeline 
import os,re
import gzip
import sys, getopt
import megalodon
from megalodon import megalodon_helper as mh
from megalodon import mods
from tqdm import tqdm
import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import ks_2samp
import scipy.stats as stats
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
import pybedtools
import pysam
from scipy.stats import pearsonr
from scipy.stats import spearmanr
# import swifter
from scipy.stats import mannwhitneyu
from scipy.stats import norm
pd.options.mode.chained_assignment = None  # default='warn'
import time
from configparser import ConfigParser
# Config file version 
HelpInfo = '''
usage: python MAGNIFIER_ExtractDat.py [Options] -c <ConfigFile> 

MAGNIFIER_ExtractDat.py is a tool for extracting modification data from megalodon dabase into multiple windows. All output data would be stored in bedformat for further analysis
Required arguments:
-c, --config\tPath to configure file recording parameters required in this example program
optional arguments:
-h, --help\tShow this help message and exit
'''
class MegaDat():
    def __init__(self, dbDir):
        # get per read methylation information
        ## write into the specific path
        self.dbDir = dbDir
        self.mods_db = mods.ModsDb(mh.get_megalodon_fn(dbDir, mh.PR_MOD_NAME), in_mem_dbid_to_uuid=True)
        # mods_txt_fp.write("\t".join(mods_db.text_field_names) + "\n")
        self.rec_tmplt = "\t".join("{}" for _ in self.mods_db.text_field_names) + "\n"
        self.bar = tqdm(
                desc="Processing Per-read Data",
                unit="per-read calls",
                total=self.mods_db.get_num_uniq_stats(),
                smoothing=0,
                dynamic_ncols=True,
            )
        # Check if bam File is ready for analysis
        if not os.path.exists(os.path.join(dbDir, "mappings.bam.bai")):
            print('%s: Prepare bamFile' % time.ctime())
            os.system('samtools sort -@ 40 {bamFile} -o {bamFile};samtools index {bamFile}'.format(bamFile=os.path.join(dbDir, "mappings.bam")))
        self.Bam = os.path.join(dbDir, 'mappings.bam')
    def ExtractDat(self, reference, start, end):
        # No classify m6A or CpG or GpC
        mod_out_text = ""
        for (chrm, strand, pos), pos_lps in self.mods_db.iter_pos_scores(convert_pos=True, pos_range=(reference, start, end)):
            self.bar.update(len(pos_lps))
            str_strand = mh.int_strand_to_str(strand)
            prev_dbid = None
            mod_bs, r_lps = [], []
            for read_dbid, mod_dbid, lp in sorted(pos_lps):
                if prev_dbid != read_dbid and prev_dbid is not None:
                    uuid = self.mods_db.get_uuid(prev_dbid)
                    # compute and store log likelihood ratios
                    with np.errstate(divide="ignore"):
                        can_lp = np.log1p(-np.exp(r_lps).sum())
                    for mod_b, r_lp in zip(mod_bs, r_lps):
                        mod_out_text += self.rec_tmplt.format(
                            uuid, chrm, str_strand, pos, r_lp, can_lp, mod_b
                        )
                    mod_bs, r_lps = [], []
                prev_dbid = read_dbid
                mod_bs.append(self.mods_db.get_mod_base(mod_dbid))
                r_lps.append(lp)
            uuid = self.mods_db.get_uuid(prev_dbid)
            # compute and store log likelihood ratios
            with np.errstate(divide="ignore"):
                can_lp = np.log1p(-np.exp(r_lps).sum())
            for mod_b, r_lp in zip(mod_bs, r_lps):
                mod_out_text += self.rec_tmplt.format(
                    uuid, chrm, str_strand, pos, r_lp, can_lp, mod_b
                )
            # mods_txt_fp = open("per_read_mods_text_%s-%s-%s" % (reference, start+100, end-100),"w")
            # mods_txt_fp.write(mod_out_text)
            # mods_txt_fp.close()
        return(mod_out_text)

def GetInfoBed3(AimDBDir, Region, ModBase, AimBase, OutDir):
    AimDB = MegaDat(AimDBDir)
    # Input Database Text file and region, get the annotated bed filer
    print("%s:Start Data Extracting" % time.ctime())
    # Content = open("{AimDBDir}/per_read_mods_text".format(AimDBDir=AimDBDir), 'r').read()
    Content = AimDB.ExtractDat(chrom, start, end)
    df_select = pd.DataFrame.from_records(list(map(lambda x:x.split("\t"), Content.split("\n")))[:-1], columns=['read_id', 'chrm', 'strand', 'pos', 'mod_log_prob', 'can_log_prob', 'mod_base'])
    df_select['pos'] = np.array(df_select['pos'], dtype=int)
    df_select['mod_log_prob'] = np.array(df_select['mod_log_prob'], dtype=float)
    df_select['can_log_prob'] = np.array(df_select['can_log_prob'], dtype=float)
    df_select['LR'] = df_select['mod_log_prob'] - df_select['can_log_prob']
    # Change The BaseName
    for i in range(len(ModBase)):
        Mod = ModBase[i]
        Aim = AimBase[i]
        df_select.loc[df_select['mod_base']==Mod, 'mod_base'] = Aim
    ## Build Bed format
    C = df_select.groupby(['read_id'])['chrm'].first()
    S = df_select.groupby(['read_id'])['pos'].apply(np.min)
    E = df_select.groupby(['read_id'])['pos'].apply(np.max)
    readN = df_select.groupby(['read_id'])['read_id'].first()
    strand = df_select.groupby(['read_id'])['strand'].first()
    BlockSite = df_select.groupby(['read_id'])['pos'].apply(lambda x:",".join([str(a) for a in x]))
    BlockLR = df_select.groupby(['read_id'])['LR'].apply(lambda x:",".join([str(a) for a in x]))
    BlockBase = df_select.groupby(['read_id'])['mod_base'].apply(lambda x:",".join([str(a) for a in x]))
    bed_Control = pd.concat([C, S, E, readN, strand, BlockSite, BlockLR, BlockBase], axis=1)
    bed_Control.columns = ['chrom', 'start','end','readN', 'strand', 'BlockSite', 'BlockLR', 'BlockBase']
    bed_Control.to_csv("{OutDir}/TMP.{Region}.bed".format(OutDir=OutDir, Region=Region), sep="\t")

def LoadData(ParaM):
    # DBdir, RegionList, ModBase, AimBase, OutDir
    DBdir, RegionFile, ModBase, AimBase, OutDir = ParaM
    AimDB = DBdir
    df_Stat = pd.read_csv(RegionFile, header=None, sep="\t")
    df_Stat.columns = ['chrom', 'start', 'end']
    df_Stat['start'] = list(map(str, list(df_Stat['start'])))
    df_Stat['end'] = list(map(str, list(df_Stat['end'])))
    df_Stat['AimRegion'] = df_Stat['chrom'] + "_" + df_Stat['start'] + "_" + df_Stat['end']
    RegionList = list(df_Stat['AimRegion'])
    for R in RegionList:
        print('%s: Region %s for %s Start' % (time.ctime(), R, OutDir))
        GetInfoBed3(AimDB, R, ModBase, AimBase, OutDir)
        print('%s: Region %s for %s End' % (time.ctime(), R, OutDir))

def BatchSelectDat(ReadName,BlockSite,ALN_df, Ref):
    # Return Whether the BlockSite in reads is the same as Reference Base
    # 1: Ref = Query
    # 0: Ref != Query
    R = list(ALN_df.loc[ALN_df['readN']==ReadName, 'ALN'])[0]
    MatchPairs = np.array(R.get_aligned_pairs(matches_only=False))
    RefPos = np.array([X[1] for X in MatchPairs])
    QPos = np.array([X[0] for X in MatchPairs])
    QSeq = R.query_sequence
    BlockSiteList = np.array(BlockSite.split(","), dtype=int)
    DecisionList = np.array([0] * len(BlockSiteList))
    BlockInRefPos_IDX = np.where(np.in1d(BlockSiteList, RefPos))[0]
    RefInBlock_IDX = np.where(np.in1d(RefPos, BlockSiteList))[0]
    RefSequence = np.array([Ref.fetch(R.reference_name, RefPos[RefInBlock_IDX[I]], RefPos[RefInBlock_IDX[I]]+1).upper() for I in range(len(RefInBlock_IDX))])
    QPos_Associated = QPos[RefInBlock_IDX]
    QSeq_Associated = np.array([None] * len(QPos_Associated))
    # Consider Query gap & MisMatch
    QSeq_Associated[np.where(QPos_Associated!=None)[0]] = np.array([QSeq[I] for I in QPos_Associated[np.where(QPos_Associated!=None)[0]]])
    DecisionList[np.where(QSeq_Associated==RefSequence)[0]] = 1
    # Save the DecisionList
    return(DecisionList)

def BGScore3(dfBed):
    # Get LR value for each single position
    PosList = np.concatenate(dfBed['BlockSite'].apply(lambda x: np.array(x.split(","), dtype=int)))
    LRList = np.concatenate(dfBed['BlockLR'].apply(lambda x: np.array(x.split(","), dtype=float)))
    ReadList = np.concatenate(dfBed.apply(lambda x: [x['readN']]*len(x['BlockSite'].split(",")), axis=1))
    Strand = np.concatenate(dfBed.apply(lambda x: [x['strand']]*len(x['BlockSite'].split(",")), axis=1))
    mod_base = np.concatenate(dfBed['BlockBase'].apply(lambda x: np.array(x.split(","))))
    Label = np.concatenate(dfBed.apply(lambda x: [a for a in x['Quality']], axis=1))
    df_Pos_Control = pd.DataFrame.from_dict({'PosList':PosList, 'LRList':LRList, 'readN':ReadList, 'Strand':Strand, 'mod_base':mod_base, 'Qual':Label})
    return(df_Pos_Control)

def AddQuality(bedFile, bamFile, Region, RefFile):
    bed_Control = pd.read_csv(bedFile, sep="\t", index_col=0)
    chrom,start,end = Region.split("_")[0], int(Region.split("_")[1]), int(Region.split("_")[2])
    ## Introduce the Quality Data
    Control_bam = pysam.AlignmentFile(bamFile, 'rb')
    ReadAlnInfo = []
    for R in Control_bam.fetch(chrom,start,end):
        ReadAlnInfo.append((R.qname, R))
    ALN_df = pd.DataFrame.from_records(ReadAlnInfo)
    if ALN_df.shape[0] == 0:
        return(pd.DataFrame())
    else:
        ALN_df.columns = ['readN', 'ALN']
        Ref = pysam.FastaFile(RefFile)
        bed_Control = bed_Control.loc[set(bed_Control.index).intersection(set(ALN_df['readN']))]
        if bed_Control.shape[0] >0:
            bed_Control['Quality'] = bed_Control.apply(lambda x: BatchSelectDat(x['readN'], x['BlockSite'], ALN_df, Ref), axis=1)
            Pos_Case = BGScore3(bed_Control)
            # Data selection
            Pos_Case_sub = pd.concat([pd.DataFrame(Pos_Case.loc[Pos_Case['Qual']==1].groupby(['PosList'])['LRList'].apply(lambda x: [a for a in x])),
                                      pd.DataFrame(Pos_Case.loc[Pos_Case['Qual']==1].groupby(['PosList'])['readN'].apply(lambda x: [a for a in x])),
                                      pd.DataFrame(Pos_Case.loc[Pos_Case['Qual']==1].groupby(['PosList'])['mod_base'].apply(lambda x: [a for a in x][0]))], axis=1)
            Pos_Case_sub['pos'] = Pos_Case_sub.index
            return(Pos_Case_sub)
        else:
            return(pd.DataFrame())

def SlidingMeanGeneral(Pos, df_summary, ReadN1='CaseReadN', DataL1="Case", ReadN2='ControlReadN', DataL2='Control'):
   # ModBase should be selected before calculation
   S,E = max(0,Pos-50), Pos+50
   mod_base = list(df_summary.loc[df_summary['Pos']==Pos, 'mod_base'])[0]
   df_tmp = df_summary.loc[(df_summary['Pos'].isin(range(S,E)))&(df_summary['mod_base']==mod_base)]
   Pos_Df1 = pd.DataFrame.from_dict({'TMP_Read':np.concatenate(np.array(df_tmp[ReadN1])), 'LRScore':np.concatenate(np.array(df_tmp[DataL1]))})
   LRDf1 = pd.DataFrame(Pos_Df1.groupby(['TMP_Read'])['LRScore'].mean())
   LRMean1 = list(LRDf1['LRScore'])
   LRRead1 = list(LRDf1.index)
   Pos_Df2 = pd.DataFrame.from_dict({'TMP_Read':np.concatenate(np.array(df_tmp[ReadN2])), 'LRScore':np.concatenate(np.array(df_tmp[DataL2]))})
   LRDf2 = pd.DataFrame(Pos_Df2.groupby(['TMP_Read'])['LRScore'].mean())
   LRMean2 = list(LRDf2['LRScore'])
   LRRead2 = list(LRDf2.index)
   return(LRMean1, LRRead1, LRMean2, LRRead2)

def DistWeight(D, T=1/60 * np.log(1/3)):
   # half weight value when distance = 50
   W = 2*np.exp(T*D) / (1+np.exp(T*D))
   return(W)

def EmpricialBetaBase(Pos, df_summary, Success='Case.m', UnSuccess='Case.um'):
   # Estimate Beta parameter from neighborhood sites
   S,E = max(0,Pos-500), Pos+500
   df_tmp = df_summary.loc[df_summary['Pos'].isin(range(S,E))]
   PosSuccess, PosFailed = list(df_tmp.loc[df_tmp['Pos']==Pos, Success])[0], list(df_tmp.loc[df_tmp['Pos']==Pos, UnSuccess])[0]
   df_bg = df_tmp.loc[df_tmp['Pos']!=Pos]
   WeightList = DistWeight(np.array(df_bg['Pos'])-Pos)
   N_success = np.sum(WeightList * np.array(df_bg[Success]))/np.sum(WeightList) + PosSuccess
   N_Failed = np.sum(WeightList * np.array(df_bg[UnSuccess]))/np.sum(WeightList) + PosFailed
   Theta_EXP = N_success / (N_success+N_Failed)
   return(Theta_EXP)

def CollectDat(Region, RawDirDict, BamFileDict, RefFile, ResDir):
    # RawDirDict: {"Case":[CG_A_dir, GC_dir], "Control":[CG_A_dir, GC_dir]} order must be the same
    # BamFileDict: {"Case":[CG_A_bam, GC_bam], "Control":[CG_A_dir, GC_dir]} order must be the same
    chrom, start, end = Region.split("_")[0], int(Region.split("_")[1]), int(Region.split("_")[2])
    print('{Time}: Region {chrom}:{start}-{end} Start'.format(Time=time.ctime(), chrom=chrom, start=start,end=end))
    # Collect Database and mod base information
    # Whether the param is right ?
    ControlList, CaseList = [], []
    for I in range(len(RawDirDict['Case'])):
        Control_bed, Case_bed = os.path.join(RawDirDict['Control'][I], 'TMP.%s.bed' % Region), os.path.join(RawDirDict['Case'][I], 'TMP.%s.bed' % Region)
        Control_bam, Case_bam = BamFileDict['Control'][I], BamFileDict['Case'][I]
        df_Control, df_Case = AddQuality(Control_bed, Control_bam, Region, RefFile), AddQuality(Case_bed, Case_bam, Region, RefFile)
        if df_Control.shape[0] * df_Case.shape[0] > 0:
            ControlList.append(df_Control)
            CaseList.append(df_Case)
    print('{Time}: Region {chrom}:{start}-{end} Quality Adding finished !'.format(Time=time.ctime(), chrom=chrom, start=start,end=end))
    if len(ControlList) * len(CaseList) > 0:
        Pos_Control_sub, Pos_Case_sub = pd.concat(ControlList, axis=0), pd.concat(CaseList, axis=0)
        Pos_Control_sub['Label'] = 'Control'
        Pos_Control_sub.index = range(Pos_Control_sub.shape[0])
        Pos_Case_sub['Label'] = 'Case'
        Pos_Case_sub.index = range(Pos_Case_sub.shape[0])
        PosIDX = list(set(Pos_Case_sub.index[np.where(Pos_Case_sub['LRList'].apply(len)>=5)]).intersection(set(Pos_Control_sub.index[np.where(Pos_Control_sub['LRList'].apply(len)>=5)])))
        if len(PosIDX) > 0:
            df_summary = pd.concat([Pos_Case_sub.loc[PosIDX, 'pos'],
                                    Pos_Case_sub.loc[PosIDX, 'mod_base'],
                                    Pos_Case_sub.loc[PosIDX,'LRList'],
                                    Pos_Case_sub.loc[PosIDX, 'readN'],
                                    Pos_Control_sub.loc[PosIDX,'LRList'],
                                    Pos_Control_sub.loc[PosIDX,'readN']], axis=1)
            df_summary.columns = ['Pos', 'mod_base', 'Case', 'CaseReadN', 'Control', 'ControlReadN']
            # Remove GpC pos that exists in CpG list
            CpGPos = list(df_summary.loc[df_summary['mod_base']=='CpG', 'Pos'])
            df_summary_sub = df_summary.drop(df_summary.loc[(df_summary['mod_base']=='GpC')&(df_summary['Pos'].isin(CpGPos))].index)
            # df_summary_AllMod.to_csv("AllMod.{chrom}_{start}_{end}.tsv".format(chrom=chrom,start=start,end=end),sep="\t")
            # Sliding Mean
            df_summary_sub['CaseLRscore'], df_summary_sub['CaseReadN2'], df_summary_sub['ControlLRscore'], df_summary_sub['ControlReadN2'] = zip(*df_summary_sub['Pos'].apply(lambda x:SlidingMeanGeneral(x, df_summary_sub)))
            # Make FisherTest
            df_summary_AllMod = df_summary_sub
            df_summary_AllMod.sort_values(['Pos'], inplace=True)
            df_summary_AllMod['cutoff2'] = df_summary_AllMod['ControlLRscore'].apply(lambda x: max(np.mean(x)+1.96*np.std(x), np.log(0.75) - np.log(0.25)))
            df_summary_AllMod['Case.m2'] = df_summary_AllMod.apply(lambda x: len(np.where(np.array(x['CaseLRscore'])>x['cutoff2'])[0]), axis=1)
            df_summary_AllMod['Case.um2'] = df_summary_AllMod.apply(lambda x:len(np.where(np.array(x['CaseLRscore'])<=x['cutoff2'])[0]), axis=1)
            df_summary_AllMod['Control.m2'] = df_summary_AllMod.apply(lambda x:len(np.where(np.array(x['ControlLRscore'])>x['cutoff2'])[0]), axis=1)
            df_summary_AllMod['Control.um2'] = df_summary_AllMod.apply(lambda x:len(np.where(np.array(x['ControlLRscore'])<=x['cutoff2'])[0]), axis=1)
            df_summary_AllMod['Pval'] = 1
            df_summary_AllMod['Pval'] = df_summary_AllMod.apply(lambda x: fisher_exact([[x['Case.m2'], x['Control.m2']], [x['Case.um2'], x['Control.um2']]], alternative='greater')[1], axis=1)
            df_summary_AllMod.to_csv("{ResDir}/AllMod.{chrom}_{start}_{end}.Calculated.tsv".format(ResDir=ResDir, chrom=chrom,start=start,end=end),sep="\t")
            print('{Time}: Region {chrom}:{start}-{end} Calculation Finished!'.format(Time=time.ctime(), chrom=chrom, start=start,end=end))
            return(df_summary_AllMod)
        else:
            return(pd.DataFrame())
    else:
        return(pd.DataFrame())

def CutBins_Old(df_summary, Step=50, FisherRes='Pval'):
    # Cut DataFrame into bins
    S,E = np.min(df_summary['Pos']), np.max(df_summary['Pos'])
    bins = [x for x in range(S-1, E, Step)]
    if bins[-1] < E:
        bins.append(E+1)
    df_summary['BinID'] = pd.cut(df_summary['Pos'], bins)
    df_summary['BinStart'] = df_summary['BinID'].apply(lambda x: int(str(x).split(',')[0][1:]))
    df_Count = pd.concat([df_summary.groupby(['BinStart'])[FisherRes].apply(lambda x: len(np.where(np.array(x)<=0.01)[0])),
                          df_summary.groupby(['BinStart'])[FisherRes].apply(lambda x: len(np.where(np.array(x)>0.01)[0]))], axis=1)
    df_Count.columns = ['Bin_Success', 'Bin_UnSuccess']
    df_Count['Pos'] = np.array(df_Count.index, dtype=int)
    df_Count['PosteriorP'] = df_Count['Pos'].apply(lambda x: EmpricialBetaBase(x,df_Count, 'Bin_Success', 'Bin_UnSuccess'))
    df_Count['Summit'] = df_Count['Pos'] + 0.5*Step
    return(df_Count)

def Runs(ParMList):
    Region, RawDirDict, BamFileDict, RefFile, ResDir, PeakDir = ParMList
    df_summary_AllMod = CollectDat(Region, RawDirDict, BamFileDict, RefFile, ResDir)
    df_Count_Old = CutBins_Old(df_summary_AllMod, Step=50)
    df_Count_Old['chrom'] = Region.split("_")[0]
    df_Count_Old['End'] = df_Count_Old['Pos']+50
    df_Count_Old[['chrom', 'Pos', 'End', 'PosteriorP']].to_csv("{PeakDir}/AccessibleScore.{Region}.bedgraph".format(PeakDir=PeakDir, Region=Region), index=False, header=None, sep="\t")
    print('{Time}: Region {Region} CallPeak Finished!'.format(Time=time.ctime(), Region=Region))
    
def main(argv):
    Parameter = {"RegionFile":None,
                 "CaseDBList":None,
                 "ControlDBList":None,
                 "ModBaseList":[],
                 "AimBaseList":[],
                 "RefFile":None, 
                 "CaseBedDir":None,
                 "ControlBedDir":None,
                 'PerReadOutDir':None,
                 'SignalOutDir':None
                }
    RequireParm1 = ['RegionFile', 'RefFile']
    RequireParm2 = ['CaseDBList', 'ControlDBList']
    RequireParm3 = ["ModBaseList", "AimBaseList"]
    RequireParm4 = ["PerReadOutDir", 'SignalOutDir']
    try:
        opts,args = getopt.getopt(argv, "hc:", ['--config'])
    except getopt.GetoptError:
        print(HelpInfo)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(HelpInfo)
            sys.exit(2)
        elif opt in ('-c', '--config'):
            configFile = arg
    # assign parameters into ParameterDict
    config = ConfigParser()
    config.read(configFile)
    Parameter['RegionFile'] = config['AttachFiles']['RegionFile']
    Parameter['RefFile'] = config['AttachFiles']['RefFile']
    Parameter['CaseDBList'] = config['AimDBList']['CaseDBList'].split(",")
    Parameter['ControlDBList'] = config['AimDBList']['ControlDBList'].split(",")
    Parameter['ModBaseList'] = [x.split(",") for x in config['BaseLabel']['ModBaseList'].split(";")]
    Parameter['AimBaseList'] = [x.split(",") for x in config['BaseLabel']['AimBaseList'].split(";")]
    Parameter["CaseBedDir"] = config['outDir']['CaseBedDir'].split(",")
    Parameter["ControlBedDir"] = config['outDir']['ControlBedDir'].split(",")
    Parameter["PerReadOutDir"] = config['outDir']['PerReadOutDir']
    Parameter["SignalOutDir"] = config['outDir']['SignalOutDir']
    # Check required Parameters
    for P in RequireParm1:
        if not os.path.exists(Parameter[P]):
            print("The {P} File {Value} not exists".format(P=P, Value=Parameter[P]))
            sys.exit(2)
    for P in RequireParm2:
        if not Parameter[P]:
            print("No {P} detected Please check python MAGNIFIER.py -h".format(P=P))
        else:
            for Dirs in Parameter[P]:
                if not os.path.exists(Dirs):
                    print("The {P} File {Value} not exists".format(P=P, Value=Parameter[P]))
                    sys.exit(2)
                    break
    for P in RequireParm3:
        if not Parameter[P]:
            print("No {P} detected Please check python MAGNIFIER.py -h".format(P=P))
        else:
            if not len(Parameter[P])==len(Parameter['CaseDBList']):
                print("The length of {P} Parameter is unequal to the case database list check python MAGNIFIER_ExtractDat.py -h".format(P=P))
                sys.exit(2)
                break
    for P in RequireParm4:
        if not Parameter[P]:
            print("No {P} detected Please check python MAGNIFIER.py -h".format(P=P))
            sys.exit(2)
            break
        else:
            if not os.path.exists(Parameter[P]):
                os.makedirs(Parameter[P])
    # Step 1: Extract per read information as bed format 
    ParaMList = [[Parameter["ControlDBList"][0], Parameter['RegionFile'], Parameter['ModBaseList'][0], Parameter['AimBaseList'][0], Parameter["ControlBedDir"][0]],
                 [Parameter["ControlDBList"][1], Parameter['RegionFile'], Parameter['ModBaseList'][1], Parameter['AimBaseList'][1], Parameter["ControlBedDir"][1]],
                 [Parameter["CaseDBList"][0], Parameter['RegionFile'], Parameter['ModBaseList'][0], Parameter['AimBaseList'][0], Parameter["CaseBedDir"][0]],
                 [Parameter["CaseDBList"][1], Parameter['RegionFile'], Parameter['ModBaseList'][1], Parameter['AimBaseList'][1], Parameter["CaseBedDir"][1]]]
    for P in ParaMList:
        OutDir = P[-1]
        if not os.path.exists(OutDir):
            os.makedirs(OutDir)
        LoadData(P)
    # Step2: collect data for MAGNIFIER calculation 
    RawDirDict = {'Case': Parameter["CaseBedDir"], 'Control':Parameter["ControlBedDir"]}
    BamFileDict = {'Case': [x+"/mappings.bam" for x in Parameter["CaseDBList"]], 'Control': [x+"/mappings.bam" for x in Parameter['ControlDBList']]}
    df_Stat = pd.read_csv(Parameter['RegionFile'], header=None, sep="\t")
    df_Stat.columns = ['chrom', 'start', 'end']
    df_Stat['start'] = list(map(str, list(df_Stat['start'])))
    df_Stat['end'] = list(map(str, list(df_Stat['end'])))
    df_Stat['AimRegion'] = df_Stat['chrom'] + "_" + df_Stat['start'] + "_" + df_Stat['end']
    RegionList = list(df_Stat['AimRegion'])
    ParMList = []
    for TestRegion in RegionList:
        ParMList.append([TestRegion, RawDirDict, BamFileDict, Parameter['RefFile'], Parameter["PerReadOutDir"], Parameter["SignalOutDir"]])
    for P in ParMList:
        Runs(P)

if __name__ == '__main__':
    main(sys.argv[1:])

