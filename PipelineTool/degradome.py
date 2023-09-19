"""
For Degradome analysis

Processing PAREsnip2 tool output
and other downstream analysis

"""

from   typing import Dict, List, Tuple, Union
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from   matplotlib.ticker import MaxNLocator
from   scipy.stats import mannwhitneyu

class DegradomeAlign(object):
    """docstring for DegradomeAlign.
        Useful atrribute for downstream aly (after object init):
            - self.tissue
            - self.transcript
            - self.position_to_check
            - self.data
            - self.select_data
            - self.formatTargetAlignTable
    """
    def __init__(self, pdfPipelineResults: pd.DataFrame , dictAlginFilePdfData: Dict[str, pd.DataFrame], tissue: str, transcript: str, position_to_check: int) -> None:
        """
            input : `pdfPipelineResults`, `{"Tissue" : AlginFileDataFrame, ...}`,
            check_tissue, check_transcript, check_position(cleavage position)
            
            - pdfPipelineResults (DataFrame) : include
                - Tissue
                - Transcript_ID
                - Cleavage_Position
            
            will process : one transcript cleavage position info
        """
        self.tissue = tissue
        self.dictAlginFilePdfData = dictAlginFilePdfData
        self.pdfPipelineResults = pdfPipelineResults
        self.transcript = transcript
        self.position_to_check = position_to_check
        try:
            data = dictAlginFilePdfData[tissue]
        except:
            raise Exception("Error: tissue not in dictAlginFilePdfData or other error")

        self.data = data
        self.select_data = data[data["transcript"] == transcript]

    def _MapPredictedMiRNACleavagePosition(self, pdfPiplineResults: pd.DataFrame, transcript_id: str, cleavage_position: int) -> bool:
        """
        Map PredictedMiRNACleavagePosition to data
        if `position` have `miRNA cleavage position` in `pipelineResults` then return `True` or `False`
        """
        mask = (pdfPiplineResults["Transcript_ID"] == transcript_id) & (pdfPiplineResults["Cleavage_Position"] == cleavage_position)
        ResultCounts = len(pdfPiplineResults[mask])
        if ResultCounts >= 1:
            return True
        else:
            return False

    @property
    def formatTargetAlignTable(self):
        """format degradome align table:
            format target align table:
                Get transcript position abundance from dataframe of this object
                * Not sum repeat position abundance
                output (dataframe; example):  
                    | transcript   | start     | abundance |
                    |--------------|-----------|-----------|
                    | transcriptA  | 1         | 3         |
                    | transcriptA  | 1         | 5         |
                    | transcriptA  | 4         | 5         |
                    .....
                    p.s. `start` is potential `cleavage position` 
        """
        data = self.data
        select_data = self.select_data
        transcript = self.transcript

        # 從 sequence 列中提取 abundance
        format_data = select_data[select_data["transcript"] == transcript]
        format_data["abundance"] = format_data["sequence"].str.extract(r"\((\d+)\)").astype(int)
        format_data = format_data[["transcript", "start", "abundance"]]
        
        return format_data
    def getTargetPositionInfo(self, dictTranscriptome_SeqLen: Dict[str, int]) -> pd.DataFrame:
        """get degradome `cleavage position info`:
            1. Sum repeat position `abundance` and get `category`
            2. output (return):
                - each position abundance and category (DataFrame)
                - set `category` = 99 if abundance = 0 
            ---

            input: 
                >>> dictTranscriptome_SeqLen = {'transcriptA' : len(sequence), ..., 'transcriptsB ...}
                >>> TargetPositionInfo = DegradomeAlign.getTargetPositionInfo(dictTranscriptome_SeqLen)
            output: 
                # if target is transcriptA
                >>> pdfTargetPositionInfo = AlyDegradomeAlign.getTargetPositionInfo(dict_transcriptome_seqLen)
                >>> pdfTargetPositionInfo
                transcript	position	abundance	category
                0	TranscriptA	1	0	99
                1	TranscriptA	2	0	99
                2	TranscriptA	3	15	2
                3	TranscriptA	4	2	3
                ...
        """
        target_transcript     = self.transcript
        format_data           = self.formatTargetAlignTable
        pdfPipelineResults    = self.pdfPipelineResults
        target_transcript_len = dictTranscriptome_SeqLen[target_transcript]
        self.target_transcript_len = target_transcript_len

        # 使用 dict 存儲 target 每個位置的 abundance 數總和
        # 先宣告將 abundance 0 覆蓋在每個位置，省去後續補 0 的步驟
        # 後面引用其他 class 處理 category 時也會用類似的邏輯只是 0 換成 ''
        dictTargetPositionAbundance = {i+1:0 for i in range(target_transcript_len)}

        for index, row in format_data.iterrows():
            # for position in range(row["start"], row["end"] + 1):
            position = row["start"]
            dictTargetPositionAbundance[position] += row["abundance"]
        degradomeCategory = PARESnip2Category(dictTargetPositionAbundance)
        dictTargetPositionCategory = degradomeCategory.CalculateCategory

        pdfTargetPositionInfo               = pd.DataFrame(list(dictTargetPositionAbundance.items()), columns=["position", "abundance"])
        pdfTargetPositionInfo["transcript"] = target_transcript
        pdfTargetPositionInfo["category"]   = list(dictTargetPositionCategory.values())
        pdfTargetPositionInfo["category"]   = pdfTargetPositionInfo["category"].apply(lambda x: x if x != '' else 99)

        pdfTargetPositionInfo["PredictedMiRNACleavagePosition"] = pdfTargetPositionInfo.apply(
            lambda row: self._MapPredictedMiRNACleavagePosition(
                pdfPipelineResults, row["transcript"], row["position"]
            ), axis=1
        )
        pdfTargetPositionInfo               = pdfTargetPositionInfo[["transcript", "position", "abundance", "category", "PredictedMiRNACleavagePosition"]]
        
        self.avarage_abundance           = degradomeCategory.avarage_abundance
        self.dictTargetPositionAbundance = dictTargetPositionAbundance
        self.dictTargetPositionCategory  = dictTargetPositionCategory
        self.pdfTargetPositionInfo       = pdfTargetPositionInfo
        
        return pdfTargetPositionInfo

class PARESnip2Category(object):
    """get degradome `category` of each potential `cleavage position`"""
    def __init__(self, dictTargetPositionAbundance: Dict[int, int]) -> Dict[int, int]:
        """object init:
            please input:
                >>> dictTargetPositionAbundance
                {1 : TranscriptA_abundance1, 2: TranscriptA_abundance2, 3: ...ce3, ...}
                >>> dictTargetPositionCategory = PARESnip2Category(dictTargetPositionAbundance)
            output:
                >>> dictTargetPositionCategory
                {1 : TranscriptA_category1, 2: TranscriptA_category2, 3: ...ry3, ...}
        """
        self.dictTargetPositionAbundance = dictTargetPositionAbundance

        # NOTE: 排除 abundance == 0 or 1 的情況
        # 分出 dict.._not_0_1 來處理	
        # NOTE: 後續計算每個切位的 平均值和 target 是否有 category 0 或 1 時會用到，
        # 皆僅考慮 abundance > 1 的 cleavage position1
        dictTargetPositionAbundance_not_0_1 = {
            cleavage_position:abundance for cleavage_position, abundance in dictTargetPositionAbundance.items() if abundance > 1
        }
        self.dictTargetPositionAbundance_not_0_1 = dictTargetPositionAbundance_not_0_1

        arrayTargetPositionAbundance = np.array(
            list(dictTargetPositionAbundance_not_0_1.values())
        )
        self.arrayTargetPositionAbundance = arrayTargetPositionAbundance
    @property
    def _CheckCategory01Exsists(self):
        """
            check category 0 or 1 exists in potential `cleavage position` \n
            (Not include abundance = 1 and 0)
        """
        # NOTE: 如果找到唯一最大 abundance (> 1) 的 cleavage position 則 category 0
        # 若超過一個則 category 1
        # 若沒有則 none
        min_Category01_or_none = None
        max_abundance = 0

        arrayTargetPositionAbundance = self.arrayTargetPositionAbundance
        if len(arrayTargetPositionAbundance) != 0:
            max_abundance = arrayTargetPositionAbundance.max()
            max_indices_size = len(np.where(arrayTargetPositionAbundance == max_abundance)[0])
            if max_indices_size == 1:
                min_Category01_or_none = 0
            elif max_indices_size > 1:
                min_Category01_or_none = 1
            else:
                raise Exception("Error: max_indices_size < 1")
        self.min_Category = min_Category01_or_none
        self.max_abundance = max_abundance
        return min_Category01_or_none, max_abundance
    @property
    def _CalculateAverageAbundance(self):
        """
            get `average` abundance of potential `cleavage position` \n
            (Not include abundance = 1 and 0)
        """
        arrayTargetPositionAbundance = self.arrayTargetPositionAbundance
        avarage_abundance = 0
        if len(arrayTargetPositionAbundance) != 0:
            avarage_abundance = arrayTargetPositionAbundance.mean()
        self.avarage_abundance = avarage_abundance
        return avarage_abundance
    @property
    def CalculateCategory(self):
        """
            get `category` of potential `cleavage position` \n
            ---
            After `object init`, please `run this function` to get result:  \n
            example:
                >>> dictTargetPositionAbundance
                {1 : TranscriptA_abundance1, 2: TranscriptA_abundance2, 3: ...ce3, ...}
                >>> PARESnip2Category = PARESnip2Category(dictTargetPositionAbundance)
                >>> dictTargetPositionCategory = PARESnip2Category.CalculateCategory
                >>> dictTargetPositionCategory
                {1 : TranscriptA_category1, 2: TranscriptA_category2, 3: ...ry3, ...}
        """
        min_category, max_abundance = self._CheckCategory01Exsists
        avarage_abundance           = self._CalculateAverageAbundance
        dictTargetPositionAbundance = self.dictTargetPositionAbundance
        dictTargetPositionCategory  = dict()
        for CleavagePosition, Abundance in dictTargetPositionAbundance.items():
            if Abundance == 1:
                category = 4
            elif Abundance == 0:
                category = ''
            else:
                # NOTE: 現在的 else 為 Abundance > 1 的情況
                if min_category is None:
                    # NOTE: 前面已經分析出結果僅可能剩 Abundance == 0 or 1 (Category '' or 4) 
                    # 的情況，卻還有 Abundance > 1 的情況
                    raise Exception("Error: Abundance > 1 but not in before step check")
                elif Abundance == max_abundance:
                    # NOTE: 當前 Abundance 最大，則套用前面算出的最小 category 0 or 1
                    category = min_category
                else:
                    # NOTE: 當前 Abundance 不是最大，則確認 Category 是否為 2, 3
                    # NOTE: 前面判斷已經分好了不可能為 '', 0, 1, 4
                    if Abundance > avarage_abundance and Abundance < max_abundance:
                        category = 2
                    elif Abundance <= avarage_abundance and Abundance > 1:
                        category = 3
                    else:
                        raise Exception("Error: if this info display, dev please check code (this is impossible))")
            dictTargetPositionCategory[CleavagePosition] = category
        self.dictTargetPositionCategory = dictTargetPositionCategory
        return dictTargetPositionCategory

class BasePlot:
    """
        Base class for PARESnip2Tplot and PARESnip2FlankTplot
    """
    def __init__(self, data: pd.DataFrame, cleavage_position: int):
        self.data              = data
        self.start             = 1
        self.end               = len(data)
        self.target_length     = len(data)
        self.cleavage_position = cleavage_position
        self.flankL            = cleavage_position - 1
        self.flankR	           = len(data) - cleavage_position
        self.target_range      = self.target_length

        self.target_abundance, \
        self.target_category  = self._get_target_cleavage_info(cleavage_position)

    def _get_target_cleavage_info(self, position):
        abundance = self.data[self.data["position"] == position]["abundance"].values[0]
        category = self.data[self.data["position"] == position]["category"].values[0]
        return abundance, category
    def _get_format_data_position_info(self, position):
        filtered_data = self.data[self.data["position"] == position]
        if filtered_data.empty:
            return "{}(abd:N/A; catgy:N/A)".format(position)
        return "{}(abd:{}; catgy:{})".format(
            position,
            filtered_data["abundance"].values[0],
            filtered_data["category"].values[0]
        )

    def print_position_info(self, flankL=None, flankR=None):
        dictPrintData = {
            "fl": flankL,
            "fr": flankR,
            "s": getattr(self, 'start', None),
            "e": getattr(self, 'end', None),
            "c": self.cleavage_position,
            "sInfo": self._get_format_data_position_info(getattr(self, 'start', None)),
            "eInfo": self._get_format_data_position_info(getattr(self, 'end', None)),
            "cInfo": self._get_format_data_position_info(self.cleavage_position),
            "r": getattr(self, 'end', 0) - getattr(self, 'start', 0)
        }

        print("CheckPositionFlank: {fl} ---- {fr}".format(**dictPrintData))
        print("CheckPositionAdjust(range:{r}) ---- {sInfo} -- {cInfo} -- {eInfo} ----".format(**dictPrintData))

    def PlotTplot(self, CleavageColor: str = "red", OtherColor: str = "skyblue", MinNormalResizeRange: list = [401, 400], resizeFigureSize: Tuple = (10, 6)):
        listColors = [CleavageColor if pos == self.cleavage_position else OtherColor for pos in self.data["position"]]
        # listColors = [
        # 	CleavageColor if pos == self.cleavage_position 
        # 	else "black" if self.data[self.data["position"] == pos]["abundance"].iloc[0] >= self.target_abundance * 0.7
        # 	else OtherColor 
        # 	for pos in self.data["position"]
        # ]
        MinResizeRange, NormalResizeRange = MinNormalResizeRange

        # 繪製柱狀圖 (fix : 401, v : len / 400) adjust figure size (fix : 10 / v : v * 5 * 1.5  , 6)
        if self.target_range <= MinResizeRange:
            plt.figure(figsize=resizeFigureSize)
        else:
            width = ((self.target_range + 1) / NormalResizeRange) * 5 * 1.5
            plt.figure(figsize=(width, 6))
        
        bars = plt.bar(self.data["position"], self.data["abundance"], color=listColors)
        plt.xlabel('Position')
        plt.ylabel('Abundance')
        plt.title('Degraded Abundance at Each Position')
        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tight_layout()
        plt.show()

    def _filter_groups_by_category(self, data: pd.DataFrame, category_threshold: int = 1) -> pd.DataFrame:
        """
            篩選掉包含 category <= threshold (default:1) 的 position_in_group (先前的 dataframe 要分好，先前都用 24 nt)
            input : datafame
            output : datafame
                - 將回傳篩選過不包含 category <= threshold (default:1) 的 group 的 dataframe
        """
        groups_to_exclude = data[data["category"] <= category_threshold]['position_in_group'].unique()

        # 從原始 DataFrame 中排除這些 position_in_group
        filtered_df = data[~data['position_in_group'].isin(groups_to_exclude)]

        return filtered_df

    class GroupCluster(object):
        """
            GroupCluster: \n
                - group `degradome target positions (have abundance)` by max_size (define_cluster_max_size; default: 24)
                - input: dataframe, max_size (default: 24)
                - output : dataframe (add column `cluster_AfterReGroup`)
        """
        def __init__(self, data: pd.DataFrame, max_size: int = 24):
            self.data                 = data.copy()
            self.max_size             = max_size
            self.group                = 0
            self.group_size           = 0
            self.group_start_position = 0
        def _group_positions(self, position: int):
            if self.group_size == self.max_size:
                self.group += 1
                self.group_size = 1
            else:
                # self.group_size += 1
                group_size = self.group_size
                group_size += position - self.group_start_position
                if group_size >= self.max_size:
                    self.group += 1
                    self.group_size = 1
            if self.group_size == 1:
                self.group_start_position = position

            return self.group
        @property
        def regroup_positions_for_cluster_data(self) -> pd.DataFrame:
            """
                regroup degradome target positions define group to cluster: \n
                    - last group dataframe `is testing simple group`
                    - `regroup` is `more powerful` than group (for `define cluster` in degradome target positions)
            """
            self.group      = 0
            self.group_size = 0
            data = self.data.sort_values('position')

            self.group_start_position = data["position"].values[0]
            data['cluster_AfterReGroup'] = data.apply(lambda x: self._group_positions(x["position"]), axis=1)
            return data

class PARESnip2Tplot(BasePlot):
    """
        PARESnip2Tplot: \n
            1. PlotTplot
            2. runTplotFilterPipeline :
                1. Background filter
                2. Cleavage cluster filter
                3. Accurate cleavage filter
    """
    def __init__(self, data: pd.DataFrame, cleavage_position: int, average_abundance: int = 0):
        super().__init__(data, cleavage_position)
        self.average_abundance = average_abundance
        flankL = self.flankL
        flankR =  self.flankR
        self.print_position_info(flankL, flankR)
        print("---------------------------")
        print(self.data[self.data["position"] == self.cleavage_position])
        print("---------------------------")
        print(self.data)
    def runTplotFilterPipeline(self) -> bool:
        """
        Executes the T-plot filtering pipeline on the data: \n
            1. Background filter
            2. Cleavage cluster filter
            3. Accurate cleavage filter
        - example : \n
            >>> objPARESnip2Tplot = PARESnip2Tplot(pdfTargetPositionInfo, Cleavage_Position, average_abundance)
            >>> boolTplotFilterOrPass = objPARESnip2Tplot.runTplotFilterPipeline()
            >>> boolTplotPassOrFilter
            True
        - Output: True/False (pass or filter)
        """
        data              = self.data
        cleavage_position = self.cleavage_position
        target_abundance  = self.target_abundance
        target_category   = self.target_category

        self.background_filter_log        = ''
        self.cleavage_cluster_filter_log  = ''
        self.accurate_cleavage_filter_log = ''
        background_filter_pass = self._background_filter(data)
        cleavage_cluster_filter_pass = self._cleavage_cluster_filter(
            data=data,
            rate=0.5, 
            define_cluster_size=24,
            clusters_threshold=2,
            ignore_miRNA_target_flank=(13, 10)
        )
        accurate_cleavage_filter_pass = self._accurate_cleavage_filter(data, slide_size=5)

        boolPassOrFilter = True
        if background_filter_pass and \
            cleavage_cluster_filter_pass and \
            accurate_cleavage_filter_pass:
            print("pass")
            boolPassOrFilter = True
        else:
            print("filter")
            boolPassOrFilter = False

        listFilterLog = [
            "background_filter : {}".format(background_filter_pass),
            "cleavage_cluster_filter : {}".format(cleavage_cluster_filter_pass),
            "accurate_cleavage_filter : {}".format(accurate_cleavage_filter_pass)
        ]
        strFilterLog = "\n".join(listFilterLog)
        print("-------------------")
        print(strFilterLog)
        print("-------------------")
        return boolPassOrFilter 

    def _background_filter(self, data: pd.DataFrame) -> bool:
        """
            使用 Mann-Whitney U 檢驗，當前 miRNA target cleavage position 的 abundance 
            是否顯著高於其他 position 的 abundance : \n
            - 檢定範圍  current miRNA cleavage position abudnace v.s. other position abundance (category 0~4)
            - alternative : greater (單尾)
            - p_value < 0.05 : True
            - 不判斷樣本量過小 or 當前要判定的 miRNA 切位 category = 4 的情況: \n
                這些情況都會回傳 True
        """
        
        target_abundance  = self.target_abundance
        target_category   = self.target_category
        Cleavage_Position = self.cleavage_position

        mask = data["category"] <= 4
        data = data[mask]

        group1 = data[data['position'] == Cleavage_Position]['abundance']
        group2 = data[
            (data['position'] != Cleavage_Position) & (data['category'] > 1)
        ]['abundance']
        group2_not_category4 = data[(data['position'] != Cleavage_Position) & (data['category'] > 1) & (data['category'] < 4)]
        # NOTE: 如果 category 為 4 則不做背景計算，或者`全部統計樣本 - category4 的樣本`的數量 < 20 也不計算
        # 即要統計樣本中 category < 4 的 data 要在 20 以上 (含)
        if target_category == 4 or len(group2_not_category4) < 20:
            if target_category == 4:
                self.background_filter_log = "True (category == 4)"
            if len(group2_not_category4) < 20:
                self.background_filter_log = "True (statics sample < 20)"
            return True

        # 使用 Mann-Whitney U 檢驗
        stat, p_value = mannwhitneyu(group1, group2, alternative='greater')
        # print(p_value)
        if p_value < 0.05:
            self.background_filter_log = f"True (p-value: {p_value})"
            return True
        else:
            self.background_filter_log = f"False (p-value: {p_value})"
            return False

    def _cleavage_cluster_filter(self, data: pd.DataFrame, rate: float = 0.5, define_cluster_size: int = 24, clusters_threshold: int = 2, ignore_miRNA_target_flank: tuple = (13, 10)) -> bool:
        """
            filter mutiple clusters result of cleavage position : \n
                - default : rate = 0.5
                    - define cluster high >= target_abundance * rate
                - define_cluster_size = 24, clusters_threshold = 2, ignore_miRNA_target_flank = (13, 10)
        """
        ignoreFlankL, ignoreFlankR = ignore_miRNA_target_flank
        target_abundance  = self.target_abundance
        target_category   = self.target_category
        cleavage_position = self.cleavage_position
        # default: 13
        ignore_range_start = cleavage_position - ignoreFlankL
        # default: 10
        ignore_range_end = cleavage_position + ignoreFlankR
        # default: 0.5
        high_degrede_signore_abundance_threshold = target_abundance * rate
        signore_check_mask = (data["position"] < ignore_range_start) | (data["position"] > ignore_range_end)
        signore_check = data[signore_check_mask].copy()
        # default: 24
        # NOTE: 將 miRNA target 的 24 nt 範圍外的 position 分群 
        # (可用 // 24+1 ，將每 24 nt 分成 1 群)
        # signore_check['position_in_group'] = signore_check['position'] // (define_cluster_size+1)
        signore_check.loc[:, 'position_in_group'] = signore_check['position'] // (define_cluster_size+1)
        mask = signore_check["abundance"] >= high_degrede_signore_abundance_threshold
        signore_check = signore_check[mask]
        # NOTE: 將 miRNA target 的 24 nt 範圍外的 position 分群後，篩選掉包含 category <= 1 的 group
        # NOTE: 由於 category <= 1，可能為重要強烈的切位，因此不排除 (不考慮於計算)
        signore_check = self._filter_groups_by_category(data=signore_check, category_threshold=1)
        if len(signore_check) > 0:
            objGroupCluster = self.GroupCluster(data=signore_check, max_size=define_cluster_size)
            # NOTE: 將 group 重新 group，目的是將 24 nt 範圍外的 position 重新分群
            signore_check = objGroupCluster.regroup_positions_for_cluster_data
            # NOTE: 過濾掉只有一個 position 的 cluster (只有一個 position 該稱為 cluster)
            signore_check = signore_check.groupby('cluster_AfterReGroup').filter(lambda x: len(x) > 1)
        print("--------------------")
        print(signore_check)
        print("--------------------")
        self.signore_check = signore_check
        self.high_degrede_signore_abundance_threshold = high_degrede_signore_abundance_threshold
        # print(miRNA_target_out_range_above_target70per)

        # grouped = signore_check.groupby('position_in_group')
        if len(signore_check) > 0:
            grouped = signore_check.groupby('cluster_AfterReGroup')
        else:
            # NOTE: 沒有 cluster 則值接回傳 True
            self.cleavage_cluster_filter_log = "True (no cluster)"
            return True
        # NOTE: 計算被當作計算的 clusters 數量
        # 當成 cluster 必須滿足與 cleavage_position_abundance 接近，且不能是 category <= 1
        # category <1 的排除在前面 _filter_groups_by_category() 已處利
        # default : 當前 cleavage_position_abundance * 0.5
        result_counts = grouped['abundance'].apply(lambda x: (x >= high_degrede_signore_abundance_threshold).sum())
        groups_above_threshold = (result_counts > 0).sum()
        # print(groups_above_threshold)

        # NOTE: 檢查計算的 clusters 有無 > 允許值
        # default : 2
        if groups_above_threshold > clusters_threshold:
            # print("filter")
            self.cleavage_cluster_filter_log = f"False ({groups_above_threshold} clusters > {clusters_threshold})"
            return False
        else:
            # print("pass")
            self.cleavage_cluster_filter_log = f"True ({groups_above_threshold} clusters <= {clusters_threshold})"
            return True


    def _accurate_cleavage_filter(self, data: pd.DataFrame, slide_size: int = 5) -> bool:
        """
            get accurate_cleavage_filter_results :
                set `max abundance` of target position `+- 5 nt (default)` \n
            input : `datafame`, slide_size (default:5) \n
            - slide_size = flank beforeAfter cleavage position `+- 5 nt (default)` \n
            output : accurate_cleavage True/False
        """
        cleavage_position = self.cleavage_position
        target_abundance  = self.target_abundance
        target_category   = self.target_category

        # NOTE: 計算 cleavage position 前後的 (default: 5) abundance 最大值，
        # 要確認其他位置是否可能才是 miRNA family 主要切位 (將其用作代表)
        # 但這些切位必須滿足上面的 miRNA 被我們先前 pipeline 找到

        # default slide_size = 5
        mask_within_range_of_PredictedMiRNACleavagePosition = \
            (data["position"] >= cleavage_position - slide_size) & \
            (data["position"] <= cleavage_position + slide_size) & \
            (data["position"] != cleavage_position) & \
            (data["PredictedMiRNACleavagePosition"] == True)

        mask_within_range_of_CleavagePosition = \
            (data["position"] >= cleavage_position - slide_size) & \
            (data["position"] != cleavage_position) & \
            (data["position"] <= cleavage_position + slide_size) & \
            (data["PredictedMiRNACleavagePosition"] == False)

        abundance_notPredict_have_max_target_abundance_BeforeAffter_SlideSize = data[mask_within_range_of_CleavagePosition]["abundance"].max()
        category_notPredict_have_max_target_abundance_BeforeAffter_SlideSize  = data[mask_within_range_of_CleavagePosition]["category"].min()

        abundance_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize = data[mask_within_range_of_PredictedMiRNACleavagePosition]["abundance"].max()
        category_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize  = data[mask_within_range_of_PredictedMiRNACleavagePosition]["category"].min()
        if target_abundance < abundance_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize:
            self.accurate_cleavage_filter_log = f"False (larger abundace on +-{slide_size}nt; categroy {category_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize})"
            return False
        else:
            if target_category <= 1:
                self.accurate_cleavage_filter_log = f"True (category == {target_category})"
                # NOTE: 僅有可能找到 == 1，因為 0 和 1 互斥
                # 若當前為 0 則其他不可能為 0 or 1
                # 當前為 1 則其他不可能為 0
                if category_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize <= 1:
                    self.accurate_cleavage_filter_log = f"True (category == {target_category})"
            else:
                if target_abundance == abundance_onAllPipelinePredictCleavagePosition_have_max_target_abundance_BeforeAffter_SlideSize:
                    self.accurate_cleavage_filter_log = f"True (same abundance on +-{slide_size}nt [Predicted])"
                else:
                    self.accurate_cleavage_filter_log = f"True (largest abundance than +-{slide_size}nt [Predicted])"
            return True

class PARESnip2FlankTplot(BasePlot):
    """
        PARESnip2FlankTplot: \n
        * PlotTplot (Only Plot Setting Range)
    """
    def __init__(self, parent, flankL: int = 10, flankR: int = 10):
        super().__init__(parent.data, parent.cleavage_position)
        self.flankL = flankL
        self.flankR = flankR
        self.start, self.end = self._adjusted_positions(flankL, flankR)
        mask_flank_range = (self.data["position"] >= self.start) & (self.data["position"] <= self.end)
        self.data = self.data[mask_flank_range]
        self.print_position_info(flankL, flankR)
        print("---------------------------")
        print(self.data[self.data["position"] == self.cleavage_position])
        print("---------------------------")
        print(self.data)

    def _adjusted_positions(self, flankL, flankR):
        start = self.cleavage_position - flankL
        end = self.cleavage_position + flankR
        start = max(start, 1)
        end = min(end, self.target_length)
        return start, end
