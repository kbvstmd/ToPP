# 1. Introduction
ToPP: Tumor online Prognostic analysis Platform <br><br>
ToPP integrated multi-omics features (mutation, CNV, gene fusion, DNA methylation, mRNA, miRNA, lncRNA, and protein expression) and clinical data of 55 tumor types, all related to prognosis. It allows validation of the most molecule’s prognosis value in pan-cancer datasets with survival data available. It also provides multiple ways for customized prognostic studies which include automatic feature optimization and univariate analysis, multivariate analysis, subgroup survival analysis, prognostic modeling, pan-cancer survival analysis and feature combination analysis from public datasets, or user-uploaded specific datasets. More sophisticated researches regarding prognosis feature regulation mechanism are also available.
# 2. univariate_analysis.R
Univariate analysis was performed with log-rank test. For continuous variables such as gene or protein expression, patients were divided into two groups according to the median value of the variable, quantile or the ‘best-cut’ value.
# 3. multivariate_analysis.R
Multivariate analysis is performed by cox regression analyses (or Cox proportional hazards model).
# 4. prognostic_modeling.R
When constructing a prognostic model, we firstly fit a naive Cox model including all covariates. In real world, we should remove redundant or irrelevant variables from the model, which describes the data by reducing variance on the expense of bias to make model more robust. Interaction between covariates and time-dependent covariates should also be considered in the modeling process. Therefore, stepwise regression is performed to select the variables of the model and is evaluated by the Akaike information criterion (AUC) value. 
# 5. combination_analysis.R
Combination analysis is to analyze the synergistic or antagonistic effects for two factors on the prognosis. These two factors can be in the same level of data such as the expression of two genes, or they can be in different levels.
# 6. pancancer.R
Pan-cancer module was designed to investigate the prognostic effects of factors in a variety of tumors. All data were acquired from the Pan-Cancer Atlas with unified normalization and standardization. 
