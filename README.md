---
title: My Analysis Tools
author: Andy6
date: 230922
---

# My Analysis Tools

## Introduction

### AnalysisTool 

- fasta Aly
    + AnalysisTool.fasta
- miRNA aly (template)
    + AnalysisTool.template
        * miRNA fasta process
    + AnalysisTool.VFold
        * predict good miRNA secondary structure by tool-Vfold (process RNAfold data)
- transcriptome Aly
    + AnalysisTool.transcriptome
        * transcriptome fasta process

### PipelineTool

- degradome data process
    + PipelineTool.degradome
        * degradome data process
            - process PARESnip2 output

### PlotTool

- create simple model for python plot tool
    + intergrate `matplotlib`, `seaborn`, `plotly`
    + please import/from PlotTool
        * use `format argument` to user friendly plot
        * now support plot include :
            - `barplot`
            - `pie`

