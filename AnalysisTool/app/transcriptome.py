"""
    Transcriptome Analysis Tool
        since 230919 (Andy6)
"""

from __future__ import annotations
from pandas.core.frame import DataFrame, Series
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
            - p.s. if dataframe include `gene_id`, set `gene_id_include` to `True` for `avoid` use `transcript_id convert gene_id`\n
            - example (gene_id_include=False): \n
                >>> data
                    transcript_id
                0	PggD10000_c0_g1_i1
                1	PggD10002_c0_g1_i1
                2	PggD10002_c0_g1_i2
                ...	...
                257743 rows × 1 columns
                >>> objGetGeneIsoformCount = GetGeneIsoformCount(data=data)
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
    def __init__(self, data: DataFrame, gene_id_include: bool = False) -> None:
        """
            GetGeneIsoformCount:
                - input: DataFrame (include transcript_id) \n
                - useful attributes :
                    self.data_add_gene_id : 
                        DataFrame (include transcript_id and gene_id)
        """
        self.dictGeneIsoFormCount         = dict()
        self.data_add_gene_id             = data
        if gene_id_include:
            pass
        else:
            self.data_add_gene_id['gene_id']  = self.data_add_gene_id.apply(lambda row: self._get_gene_id(row), axis=1)
    def _get_gene_id(self, row: Series) -> List:
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
        for index, row in self.data_add_gene_id.iterrows():
            gene_id = row['gene_id']
            if gene_id not in self.dictGeneIsoFormCount:
                self.dictGeneIsoFormCount[gene_id] = 0
            self.dictGeneIsoFormCount[gene_id] += 1
        return self.dictGeneIsoFormCount
    def _add_isoform_count(self, row: Series) -> List:
        """
            Get gene isofrom
            input :
                row : Series (include gene_id)
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
        self.data_add_gene_id['isoform_count'] = self.data_add_gene_id.apply(
            lambda row: self._add_isoform_count(row),
            axis=1
        )
        if how != "gene_and_isoform_count":
            Result = self.data_add_gene_id
        else:
            Result = self.data_add_gene_id[['gene_id', 'isoform_count']].drop_duplicates()
        return Result

class AnalayGeneIsoForm(object):
    """
        AnalayGeneIsoForm:
            Analay `genes distribution` of `difference isoform count`
        ---
        - input: DataFrame (include gene_id and isoform_count)
        - output:  `gene_id` and `isoform_count` to dict or DataFrame (GeneIsoFormCount)
    """
    def __init__(self, data: DataFrame) -> None:
        self.data = data
        self.dictIsoFormCountGeneDistribution = dict()
    @property
    def GetIsoFormCountGeneDistribution(self) -> Dict:
        """
            GetIsoFormCountGeneDistribution:
            ---
                input:
            ---
                    DataFrame (include gene_id and isoform_count)

            ---
                output: 
            ---
                    dict (key: isoform_count, value: gene_count)
                    
            ---
                example:
            ---
                >>> data
                    gene_id	isoform_count
                171509	PggD4975_c16_g1	24
                89481	PggD14842_c13_g1	22
                ...	...	...
                221328 rows × 2 columns
                >>> objAnalayGeneIsoForm = AnalayGeneIsoForm(data)
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
        for index, row in self.data.iterrows():
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

class StaticAnnotationGene(object):
    """
        StaticAnnotationGene:
        ---
            static annotation gene
        ---
                `input`: 
        ---
                    DataFrame (include gene_id and annotation)

        ---
                `format` and `return results`:
        ---
                    dict (key: annotation, value: gene_count; Static) \n
                    or DataFrame (annotation, gene_count; to_DataFrame)
    """
    def __init__(self, data: DataFrame) -> None:
        self.data = data
        self.dictAnnotationGene = dict()
    @property
    def Static(self) -> Dict[str, int]:
        """
            StaticAnnotationGene:
            ---
                - static annotation gene
                - input: DataFrame (include gene_id and annotation)
                - output: 
                    - dict (key: annotation, value: gene_count)
            ---

            example
            ---
                >>> data
                    gene_id	annotation
                0	PggD13791_c37_g1	XM_020715060.1 PREDICTED: Phalaenopsis equestr...
                1	PggD3380_c11_g1	XM_020724321.1 PREDICTED: Phalaenopsis equestr...
                2	PggD7976_c8_g1	XM_020720713.1 PREDICTED: Phalaenopsis equestr...
                ...	...	...
                87959 rows × 2 columns
                >>> objStaticAnnotationGene = StaticAnnotationGene(data)
                >>> objStaticAnnotationGene.Static
                {'XM_020715060.1 PREDICTED: Phalaenopsis equestris midasin-like (LOC110017904), mRNA': 2,
                'XM_020724321.1 PREDICTED: Phalaenopsis equestris auxin transport protein BIG (LOC110024389), mRNA': 2,
                'XM_020720713.1 PREDICTED: Phalaenopsis equestris uncharacterized LOC110021972 (LOC110021972), transcript variant X1, mRNA': 1,
                ...}
        """
        for index, row in self.data.iterrows():
            annotation = row['annotation']
            if annotation not in self.dictAnnotationGene:
                self.dictAnnotationGene[annotation] = 0
            self.dictAnnotationGene[annotation] += 1
        return self.dictAnnotationGene
    @property
    def to_DataFrame(self):
        """
            FormatResults to DataFrame:
            ---
                - After `Static`
                - output: DataFrame (annotation, gene_count)
            ---

            example
            ---
                >>> objStaticAnnotationGene = StaticAnnotationGene(data)
                >>> objStaticAnnotationGene.Static
                ...
                >>> objStaticAnnotationGene.to_DataFrame
                    annotation	gene_count
                0	XM_020715060.1 PREDICTED: Phalaenopsis equestr...	2
                1	XM_020724321.1 PREDICTED: Phalaenopsis equestr...	2
                2	XM_020720713.1 PREDICTED: Phalaenopsis equestr...	1
                ...	...	...
                30511 rows × 2 columns


        """
        data = pd.DataFrame(
            {
                'annotation' : self.dictAnnotationGene.keys(),
                'gene_count' : self.dictAnnotationGene.values()
            }
        )
        return data

