# scMetabolism
`scMetabolism` is a R package for quantifying metabolism activity at the single-cell resolution
![Screenshot](https://github.com/wu-yc/scMetabolism/raw/main/logo.jpg)

## Requirements
    install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
    devtools::install_github("YosefLab/VISION")
    

## Install
    devtools::install_github("wu-yc/scMetabolism")

## Quick Start
`scMetabolism` generally supports the quantification and visualization of metabolism at the single-cell resolution. 

`scMetabolism` currently supports human scRNA-seq data.


### 1. Load packages and demo data
The demo data is the dataset of Peripheral Blood Mononuclear Cells (PBMC) from 10X Genomics open access dataset (~2,700 single cells, also used by Seurat tutorial). The demo Seurat object can be downloaded from [here](https://figshare.com/articles/dataset/scMetabolism_-_pbmc_demo_rda/13670038).


    load(file = "pbmc_demo.rda")
    
    library(scMetabolism)
    library(ggplot2)
    library(rsvd)


### 2. Quantify single-cell metabolism with Seurat (Recommended)
    countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")

`obj` is a Seurat object containing the UMI count matrix. 

`method` supports `VISION`, `AUCell`, `ssgsea`, and `gsva`, which VISION is the default method.

`imputation` allows users to choose whether impute their data before metabolism scoring.

`ncores` is the number of threads of parallel computation.

`metabolism.type` supports `KEGG` and `REACTOME`, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.

To extract the metabolism score, just run `metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score`, where `metabolism.matrix` is the matrix.

### 3. Visualize 
#### Dimplot

    DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)

`countexp.Seurat` is a Seurat object containing the UMI count matrix. 

`pathway` is the pathway of interest to visualize. 

`dimention.reduction.type` supports `umap` and `tsne`.

`dimention.reduction.run` allows users to choose whether re-run the dimention reduction of the given Seurat object.

`size` is the dot size in the plot.

This function returns a ggplot object, which can be DIY by users.

![Screenshot](https://github.com/wu-yc/scMetabolism/raw/main/scmetab_dim.png)

#### Dot plot

    input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
    DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", norm = "y")

`obj` is a Seurat object containing the UMI count matrix. 

`pathway` is the pathway of interest to visualize. 

`phenotype` is the one of the features contained in the metadata in the Seurat object.

`norm` refers to scale the value according to row or column. Users can choose "x", "y", and "na".

This function returns a ggplot object, which can be DIY by users.

![Screenshot](https://github.com/wu-yc/scMetabolism/raw/main/scmetab_dot.png)

#### Box plot

    BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", ncol = 1)

`obj` is a Seurat object containing the UMI count matrix. 

`pathway` is the pathway of interest to visualize. 

`phenotype` is the one of the features contained in the metadata in the Seurat object.

`ncol` refers to the column number per row.

This function returns a ggplot object, which can be DIY by users.

![Screenshot](https://github.com/wu-yc/scMetabolism/raw/main/scmetab_box.png)

### 4. Quantify single-cell metabolism WITHOUT Seurat (Not recommended)
scMetabolism also supports quantifying metabolism independent of Seurat. 

    metabolism.matrix<-sc.metabolism(countexp = countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")

`countexp` is a data frame of UMI count matrix (col is cell ID, row is gene name). 

`method` supports `VISION`, `AUCell`, `ssgsea`, and `gsva`, which VISION is the default method.

`imputation` allows users to choose whether impute their data before metabolism scoring.

`ncores` is the number of threads of parallel computation.

`metabolism.type` supports `KEGG` and `REACTOME`, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.

## Citations
**_scMetabolism_**

- Yingcheng Wu, Shuaixi Yang, Jiaqiang Ma, Zechuan Chen, Guohe Song, Dongning Rao, Yifei Cheng, Siyuan Huang, Yifei Liu, Shan Jiang, Jinxia Liu, Xiaowu Huang, Xiaoying Wang, Shuangjian Qiu, Jianmin Xu, Ruibin Xi, Fan Bai, Jian Zhou, Jia Fan, Xiaoming Zhang, and Qiang Gao. Spatiotemporal Immune Landscape of Colorectal Cancer Liver Metastasis at Single-Cell Level. Cancer Discovery. 2021.

**_Genesets and algorithms_**
1. DeTomaso D, et al. Nat Commun. 2019 Sep 26;10(1):4376.
2. Aibar S, et al. Nat Methods. 2017 Nov;14(11):1083-1086.
3. Xiao Z, et al. Nat Commun. 2019 Aug 21;10(1):3763.
4. Hänzelmann S, et al. BMC Bioinformatics. 2013 Jan 16;14:7.
5. George C. Linderman, et al. bioRxiv 2019.


## Online version of scMetabolism
http://cancerdiversity.asia/scMetabolism/


## Contact

Qiang Gao, MD, PhD

Department of Liver Surgery and Transplantation, Liver Cancer Institute, Zhongshan Hospital, Fudan University, Shanghai, China

gaoqiang@fudan.edu.cn


Any technical question please contact Yingcheng Wu (wuyc@usa.com).

Copyright (C) 2020-2021 Gao Lab @ Fudan University.



