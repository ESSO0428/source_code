"""
    Transcriptome Analysis Tool
        since 230919 (Andy6)
"""

from __future__ import annotations
from pandas.core.frame import DataFrame
import pandas as pd
from   matplotlib import pyplot as plt
import seaborn as sns
from   typing import List, Dict, Tuple, Union, Optional, Any

class GetGeneIsoformCount(object):
    """
        GetGeneIsoformCount
            proceess transcriptome dataframe:
            ---
            - get `gene_id` and `isoform_count` from `transcript_id` \n
            - input: `DataFrame` (include transcript_id) \n
            - example: \n
                >>> pdfTranscriptome
                    transcript_id
                0	PggD10000_c0_g1_i1
                1	PggD10002_c0_g1_i1
                2	PggD10002_c0_g1_i2
                ...	...
                257743 rows × 1 columns
                >>> objGetGeneIsoformCount = GetGeneIsoformCount(pdfTranscriptome)
                >>> objGetGeneIsoformCount.CountGeneIsoform()
                {'PggD10000_c0_g1': 2,
                'PggD10002_c0_g1': 4,
                'PggD10002_c10_g1': 2,
                ...}
                >>> objGetGeneIsoformCount.GetResults(how="default")
                    transcript_id	gene_id	isoform_count
                0	PggD10000_c0_g1_i1	PggD10000_c0_g1	1
                1	PggD10002_c0_g1_i1	PggD10002_c0_g1	2
                2	PggD10002_c0_g1_i2	PggD10002_c0_g1	2
                ...	...	...	...
                257743 rows × 3 columns
                >>> objGetGeneIsoformCount.GetResults(how="gene_and_isoform_count")
                    gene_id	isoform_count
                0	PggD10000_c0_g1	1
                1	PggD10002_c0_g1	2
                ...	...	...
                221328 rows × 2 columns
    """
    def __init__(self, pdfTranscriptome: DataFrame):
        """
            GetGeneIsoformCount:
                - input: DataFrame (include transcript_id) \n
                - useful attributes :
                    self.pdfTranscriptome : 
                        DataFrame (include transcript_id and gene_id)
        """
        self.dictGeneIsoFormCount = dict()
        self.pdfTranscriptome = pdfTranscriptome
        self.pdfTranscriptome['gene_id'] = self.pdfTranscriptome.apply(lambda row: self._get_gene_id(row), axis=1)
    def _get_gene_id(self, row: pd.Series) -> List:
        """
            Get gene_id from transcript_id
        """
        transcript_id     = row['transcript_id']
        gene_id = '_'.join(transcript_id.split('_')[:-1])
        return gene_id
    def CountGeneIsoform(self) -> Dict[str, int]:
        """
            Count gene isoform
            input : DataFrame (include gene_id)
        """
        for index, row in self.pdfTranscriptome.iterrows():
            gene_id = row['gene_id']
            if gene_id not in self.dictGeneIsoFormCount:
                self.dictGeneIsoFormCount[gene_id] = 0
            self.dictGeneIsoFormCount[gene_id] += 1
        return self.dictGeneIsoFormCount
    def _add_isoform_count(self, row: pd.Series) -> List:
        """
            Get gene isofrom
            input :
                row : pd.Series (include gene_id)
            output :
                row : list (list of gene isoform counts for each gene_id)
        """
        gene_id = row['gene_id']
        isoform_count = self.dictGeneIsoFormCount[gene_id]
        return isoform_count
    def GetResults(self, how: str = "default") -> DataFrame:
        """
            Get result  \n
                return : `DataFrame`  \n
                Args :
                    `default` : return all columns  \n
                    `gene_and_isoform_count` : return gene_id and isoform_count
        """
        self.pdfTranscriptome['isoform_count'] = self.pdfTranscriptome.apply(
            lambda row: self._add_isoform_count(row),
            axis=1
        )
        if how != "gene_and_isoform_count":
            Result = self.pdfTranscriptome
        else:
            Result = self.pdfTranscriptome[['gene_id', 'isoform_count']].drop_duplicates()
        return Result

class AnalayGeneIsoForm(object):
    def __init__(self, pdfTranscriptomeGeneLevel: DataFrame) -> None:
        self.pdfTranscriptomeGeneLevel = pdfTranscriptomeGeneLevel
        self.dictIsoFormCountGeneDistribution = dict()
    @property
    def GetIsoFormCountGeneDistribution(self) -> Dict:
        """
            GetIsoFormCountGeneDistribution:
            ---
                output: 
            ---
                    dict (key: isoform_count, value: gene_count)
                    
            ---
                example:
            ---
                >>> pdfTranscriptomeGeneLevel
                    gene_id	isoform_count
                171509	PggD4975_c16_g1	24
                89481	PggD14842_c13_g1	22
                ...	...	...
                221328 rows × 2 columns
                >>> objAnalayGeneIsoForm = AnalayGeneIsoForm(pdfTranscriptomeGeneLevel)
                >>> objAnalayGeneIsoForm.GetIsoFormCountGeneDistribution
                {24: 2,
                22: 2,
                ...}
                >>> objAnalayGeneIsoForm.GetFormatDistribution(how="DataFrame")
                    isoform_count	gene_count
                0	24	1
                1	22	1
                2	21	1
                ...
                >>> objAnalayGeneIsoForm.GetFormatDistribution(how="Dict")
                {24: 2,
                22: 2,
                ...}
        """
        for index, row in self.pdfTranscriptomeGeneLevel.iterrows():
            isoform_count = row['isoform_count']
            if isoform_count not in self.dictIsoFormCountGeneDistribution:
                self.dictIsoFormCountGeneDistribution[isoform_count] = 0
            self.dictIsoFormCountGeneDistribution[isoform_count] += 1
        return self.dictIsoFormCountGeneDistribution
    def GetFormatDistribution(self, how="DataFrame") -> DataFrame | Dict:
        """
            GetFormatDistribution:
            ---
                - After `GetIsoFormCountGeneDistribution()`
                - Args: how (str) : "DataFrame" or "Dict"
                    - for your return result type
                - output: 
                    - DataFrame (isoform_count, gene_count)
                    - or dict (key: isoform_count, value: gene_count)
        """
        if how != "Dict":
            data = pd.DataFrame(
                {
                    'isoform_count' : self.dictIsoFormCountGeneDistribution.keys(),
                    'gene_count'    : self.dictIsoFormCountGeneDistribution.values()
                }
            )
            return data
        else:
            return self.dictIsoFormCountGeneDistribution
