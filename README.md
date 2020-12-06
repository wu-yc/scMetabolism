# scMetabolism
`scMetabolism` is a R package for quantifying metabolism activity at the single-cell resolution

## Requirements
    install.packages(c("data.table", "wesanderson", "Seurat", "VISION", "AUCell", "GSEABase", "GSVA", "ggplot2"))
    

## Install
    install.packages(devtools)
    install_github("wu-yc/scMetabolism")

## Quick Start
`scMetabolism` generally supports the quantification and visualization of metabolism at the single-cell resolution. 

### 1. Quantify single-cell metabolism with Seurat (Recommended)
    countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")

`obj` is a Seurat object containing the UMI count matrix. 

`method` supports `VISION`, `AUCell`, `ssgsea`, and `gsva`, which VISION is the default method.

`imputation` allows users to choose whether impute their data before metabolism scoring.

`ncores` is the number of threads of parallel computation.

`metabolism.type` supports `KEGG` and `REACTOME`, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.

To extract the metabolism score, just run `metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score`, where `metabolism.matrix` is the matrix.

### 2. Visualize 
    plot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = T, size = 1)

`obj` is a Seurat object containing the UMI count matrix. 

`pathway` is the pathway of interest to visualize. 

`dimention.reduction.type` supports `UMAP` and `tSNE`.

`dimention.reduction.run` allows users to choose whether re-run the dimention reduction of the given Seurat object.

`size` is the dot size in the plot.

This function returns a ggplot object, which can be DIY by users.

![Screenshot](https://github.com/wu-yc/scMetabolism/raw/main/scMetabolism_demo.png)


### 3. Quantify single-cell metabolism WITHOUT Seurat (Not recommended)
scMetabolism also supports quantifying metabolism independent of Seurat. 

    metabolism.matrix<-sc.metabolism(countexp = countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")

`countexp` is a data frame of UMI count matrix (col is cell ID, row is gene name). 

`method` supports `VISION`, `AUCell`, `ssgsea`, and `gsva`, which VISION is the default method.

`imputation` allows users to choose whether impute their data before metabolism scoring.

`ncores` is the number of threads of parallel computation.

`metabolism.type` supports `KEGG` and `REACTOME`, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.

## Citations
1. DeTomaso D, et al. Nat Commun. 2019 Sep 26;10(1):4376.
2. Aibar S, et al. Nat Methods. 2017 Nov;14(11):1083-1086.
3. Xiao Z, et al. Nat Commun. 2019 Aug 21;10(1):3763.
4. Hänzelmann S, et al. BMC Bioinformatics. 2013 Jan 16;14:7.
5. George C. Linderman, et al. bioRxiv 2019.


## Contact

Qiang Gao, MD, PhD

Cheung Kong Scholar Distinguished Professor, Department of Liver Surgery and Transplantation, Liver Cancer Institute, Zhongshan Hospital, Fudan University, Shanghai, China

gaoqiang@fudan.edu.cn

Any technical question please contact wuyc@usa.com (Yingcheng Wu)

Copyright (C) 2020 Gao Lab @ Fudan University.



