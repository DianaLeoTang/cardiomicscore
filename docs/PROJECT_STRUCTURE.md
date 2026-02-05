# CardiOmicScore 项目结构与文件说明

本文档描述项目「AI-based multiomics profiling for personalized prediction of cardiovascular disease: A prospective UK Biobank study」的目录结构、各文件功能与数据/模型流程。

---

## 一、项目概述

- **目的**：基于 UK Biobank 多组学（代谢组、蛋白组、基因组 PRS）与临床变量，用深度学习与 Cox 等模型预测 6 种心血管结局（CAD、卒中、心衰、房颤、外周动脉病、静脉血栓）。
- **技术栈**：Python（数据与模型）、R（竞争风险、可视化）。
- **产出**：训练好的 MetNet/ProNet、C-index/ΔC-index、SHAP 重要性、生存分析与论文图表。

---

## 二、根目录与顶层文件

| 路径 | 类型 | 功能说明 |
|------|------|----------|
| `README.md` | 文档 | 项目说明、结构指引、环境与依赖、联系方式 |
| `LICENSE` | 文档 | MIT 许可证 |
| `src/` | 目录 | 论文插图与漫画解读（中英文 PDF、封面图） |

---

## 三、配置目录 `config/`

存放各深度学习模型（OmicsNet）的**最终超参数**，以 JSON 形式保存。

| 路径 | 功能 |
|------|------|
| `config/Metabolomics/final_param.json` | 代谢组学模型（MetNet）最终参数：输入维度、共享/跳跃/预测头 MLP 结构、dropout、激活函数、数据路径、训练轮数等 |
| `config/Metabolomics_no_statins/final_param.json` | 排除他汀使用者后的代谢组学模型参数 |
| `config/Proteomics/final_param.json` | 蛋白组学模型（Proteomics/ProNet）最终参数 |

这些 JSON 被 `model/training/` 与 `model/hyperparameter_tuning/` 中的脚本读取，用于训练或调参。

---

## 四、数据目录 `data/`

### 4.1 `data/gwas_summary_statistics/`

用于**多基因风险评分（PRS）**的 GWAS 汇总统计与评分文件。

| 文件 | 功能 |
|------|------|
| `donwload_summary_statistics.sh` | 从 PGS Catalog 下载各结局的 PRS 评分文件（CAD、卒中、房颤、心衰、PAD、VTE），并解压 |
| `PGS00xxxxx.txt` | 原始 PGS 评分文件（按结局） |
| `PGS00xxxxx_formatted.txt` | 经 `data_loader/DataPreparation/OmicsGenerator/PRSProcessing/` 预处理后的格式，供 PRS 计算使用 |

### 4.2 其他数据目录（README 中提及，当前仓库未包含）

- `data/processed/`：中间处理结果（组学、协变量、结局等）
- `data/split_seed-XXX/`：按种子划分并完成插补后的 train/val/test 数据集（模型训练主数据）
- `data/split_seed_raw-250901/`：未做缺失插补的原始划分，用于基线特征表等描述性分析

---

## 五、数据准备模块 `data_loader/`

负责从原始 UK Biobank 数据生成**特征（X）**与**结局（y, e）**，并划分/插补得到最终数据集。

### 5.1 顶层入口

| 文件 | 功能 |
|------|------|
| `data_loader.py` | 定义 `UKBiobankDataset`（PyTorch Dataset）和 `UKBiobankDataLoader`：按 `data_dir`、`dataset_type`（train/val/internal_test/external_test）、`predictor_set` 读取 X/y/e 的 feather 文件，供训练与评估使用 |
| `data_preparation.py` | `UKBiobankData()` 类：合并临床、基因组、代谢组等预测因子，定义排除规则、划分策略（train/val/test）、缺失值插补（如 KNN/SimpleImputer）、标准化，并写出各 `predictor_set` 的 X/y/e feather 到 `data/split_seed-{seed}/` |
| `data_preparation_raw.py` | 基于原始数据生成**未插补**的划分（如 `split_seed_raw-250901`），用于基线特征表等 |

### 5.2 组学生成 `DataPreparation/OmicsGenerator/`

#### 代谢组学 `NMRDataProcessing/`

| 文件 | 功能 |
|------|------|
| `S0_ReadRawNMRData.r` | R 脚本：读取原始 NMR 代谢组学数据，做初步整理与输出 |
| `S1_MetabolomicsQC.py` | 代谢组学 QC：读取上步输出，剔除高缺失代谢物等，输出 `Metabolomics.csv` 等至 `data/processed/omics/` |

#### 蛋白组学 `OlinkDataProcessing/`

| 文件 | 功能 |
|------|------|
| `S0_ReadRawOlinkData.py` | 读取原始 Olink 蛋白组学数据 |
| `S1_ProteomicsQC.py` | 蛋白组学 QC（缺失率过滤等），输出至 processed |

#### PRS 处理 `PRSProcessing/`

| 文件 | 功能 |
|------|------|
| `S0_summary_statistics_preprocessing.py` | 将 GWAS 汇总统计预处理为 PRS 计算所需格式（如生成 `*_formatted.txt`） |
| `S1_ukbrapr_prs_calculation.r` | 在 UK Biobank 基因型上计算 PRS（调用预处理后的评分文件） |
| `S2_merge_prs.py` | 将各结局的 PRS 合并到统一表，供 `data_preparation.py` 合并进 X |

### 5.3 结局生成 `DataPreparation/OutcomesGenerator/`

| 文件 | 功能 |
|------|------|
| `S0_OutcomesInfo.py` | 从 UK Biobank 诊断与随访信息中提取结局相关字段与编码（如 ICD-10、日期），定义 6 种心血管结局的识别规则 |
| `S1_OutcomesGenerator.py` | 基于上述规则计算各结局的随访时间、事件发生（y、e），与 `data_preparation.py` 使用的列名（如 `bl2{cad}_yrs`、`cad`）对应 |

### 5.4 预测因子生成 `DataPreparation/PredictorsGenerator/`

每个脚本从 UK Biobank 原始字段生成一类临床/生活方式预测变量，供 `data_preparation.py` 在 `get_predictors_dict()` 中引用并合并。

| 文件 | 功能 |
|------|------|
| `Demographic.py` | 人口学：年龄、性别、种族等 |
| `PhysicalMeasurements.py` | 体格：身高、体重、BMI、腰围、腰臀比、血压等 |
| `DiseaseHistory.py` | 疾病史：从诊断编码与日期中提取各病种病史（含 `diag_min_dates` 等工具函数） |
| `FamilyHistory.py` | 家族史：心血管、卒中、高血压、糖尿病等 |
| `Lifestyle.py` | 生活方式：吸烟、饮酒、运动等 |
| `LifestyleDiet.py` | 饮食相关变量 |
| `LifestyleSocial.py` | 社会活动等 |
| `MedicationHistory.py` | 用药史：降压、降脂等 |
| `Biofluids.py` | 血液指标：血常规、生化等（对应 `blood_count_*`、`blood_biochem_*`） |
| `GenotypeQC.py` | 基因型 QC，供 PRS 前筛选样本/位点 |

`data_preparation.py` 中 `predictor_set_mapping` 将上述模块组合成：`AgeSex`、`Clinical`、`PANEL`、`ASCVD`、`SCORE2`、`Genomics`、`Metabolomics`、`Proteomics`、`Metabolomics_no_statins` 等。

---

## 六、模型模块 `model/`

### 6.1 核心模型与训练器

| 文件 | 功能 |
|------|------|
| `model.py` | **OmicsNet**：多任务 MLP 网络。含共享 MLP、每个结局一个 TaskSpecificMLP（跳跃连接 MLP + 预测头 MLP）。支持从 legacy kwargs 构建 MLP；输入为组学特征，输出各结局的 logit。 |
| `trainer.py` | **BaseTrainer**：通用训练循环（epoch、早停、验证 C-index、保存 checkpoint）、测试集评估（内部/外部）、AUC/C-index/PRAUC/混淆矩阵等。**WeightedMLPTrainer** / **MLPTrainer**：多任务 BCEWithLogitsLoss，支持按结局加权，与 OmicsNet 配合使用。 |

### 6.2 超参数调优 `model/hyperparameter_tuning/`

| 文件 | 功能 |
|------|------|
| `metabolomics.py` | 使用 Optuna 对代谢组学 OmicsNet（MetNet）进行超参搜索（网络结构、学习率、dropout 等），以验证集 C-index 为目标 |
| `proteomics.py` | 对蛋白组学 OmicsNet（ProNet）进行 Optuna 超参搜索 |

### 6.3 正式训练 `model/training/`

| 文件 | 功能 |
|------|------|
| `metabolomics_5times.py` | 使用 `config/Metabolomics/final_param.json`，多随机种子（如 5 次）训练 MetNet；每轮训练+验证+测试，计算 SHAP、保存预测分数与模型 |
| `metabolomics_no_statins_5times.py` | 同上，但使用排除他汀人群的数据与对应 config |
| `proteomics_5times.py` | 多种子训练 ProNet，并产出分数与 SHAP |

### 6.4 基线模型 `model/baseline/`

| 文件 | 功能 |
|------|------|
| `clinical_scores.py` | 临床评分（如 ASCVD、SCORE2）在外部测试集上的点估计 C-index 与 Bootstrap 置信区间；`ClinicalScoreEvaluator` 基类与具体评分子类 |
| `machine_learning_metabolomics.py` | 代谢组学上的传统 ML（如 XGBoost、随机森林）训练与评估 |
| `machine_learning_proteomics.py` | 蛋白组学上的传统 ML 训练与评估 |

### 6.5 Cox 比例风险模型 `model/coxph/`

在组合好的预测因子（临床、PRS、MetNet/ProNet 分数等）上拟合 Cox 模型，评估 C-index、ΔC-index 与发病概率，支持 Bootstrap。

| 文件 | 功能 |
|------|------|
| `coxph.py` | **PointEstimateEvaluator**：加载各 predictor 组合的 X 与 y/e，拟合 Cox，计算 C-index，保存发病概率等 |
| `coxph_bootstrap.py` | 在 Cox 点估计基础上做 Bootstrap，得到 C-index、ΔC-index 的置信区间 |
| `coxph_event_2yrs.py` | 2 年事件发生相关的 Cox/二分类分析 |
| `coxph_event_2yrs_bootstrap.py` | 上述 2 年分析的 Bootstrap 版本 |
| `coxph_individual_biomarker.py` | 单一生化标志物（或单变量）Cox 分析 |
| `coxph_individal_biomarker_metabolite_bootstrap.py` | 代谢物单变量 Cox 的 Bootstrap（注意文件名拼写 individual） |
| `coxph_individal_biomarker_protein_bootstrap.py` | 蛋白单变量 Cox 的 Bootstrap |
| `coxph_LDL_biomarker.py` | LDL 相关标志物组合的 Cox |
| `coxph_LDL_biomarker_bootstrap.py` | 同上，Bootstrap |
| `coxph_LDL_biomarker_no_statins.py` | 排除他汀人群的 LDL 标志物 Cox |
| `coxph_LDL_biomarker_no_statins_bootstrap.py` | 同上，Bootstrap |
| `coxph_LDL_biomarker_no_statins_panel.py` | 排除他汀 + PANEL 的 LDL 相关分析 |
| `coxph_LDL_biomarker_no_statins_panel_bootstrap.py` | 同上，Bootstrap |
| `coxph_ml.py` | 将机器学习（如 MetNet/ProNet）输出分数作为 Cox 的输入进行拟合与比较 |
| `coxph_ml_bootstrap.py` | 同上，Bootstrap |
| `coxph_single_omics_vs_multiomics.py` | 单组学 vs 多组学（如 Clinical+Met+Pro）的 Cox 比较 |
| `coxph_single_omics_vs_multiomics_bootstrap.py` | 同上，Bootstrap |
| `coxph_subgroup.py` | 亚组分析（年龄、性别、用药等）的 Cox |
| `coxph_subgroup_bootstrap_*.py` | 各亚组的 Bootstrap（age/sex/meds_ht/meds_lipid） |
| `coxph_train_val_internal_test.py` | 在 train/val/internal test 上的 Cox 评估流程 |
| `coxph_train_val_internal_test_bootstrap.py` | 同上，Bootstrap |
| `coxph_website.py` | 为在线计算器或网站提供的简化 Cox/评分接口 |

### 6.6 竞争风险（Fine-Gray）`model/fine_gray/`

在竞争风险框架下比较不同预测因子组合的效应（考虑死亡竞争事件）。

| 文件 | 功能 |
|------|------|
| `fine_gray.R` | 加载 y/e/死亡时间与各 predictor 数据，拟合 Fine-Gray 模型，输出子分布 HR 等 |
| `fine_gray_bootstrap_cad_stroke.R` | CAD、卒中的 Fine-Gray Bootstrap |
| `fine_gray_bootstrap_hf_af.R` | 心衰、房颤的 Fine-Gray Bootstrap |
| `fine_gray_bootstrap_pad_vte.R` | PAD、VTE 的 Fine-Gray Bootstrap |

---

## 七、分析脚本 `analysis/`

将模型与 Cox 输出汇总为表格与中间结果，供 `visualization/` 画图制表。

| 路径 | 功能 |
|------|------|
| `MetricMerge/MetricMerge.ipynb` | 读取 C-index 点估计与 Bootstrap CI（如 `cindex_summary.csv`、`cindex_bootstrap_ci_summary.csv`），按 outcome、baseline_model、comparison_model 合并，整理为 c_index 与 delta_c_index 的长表，保存为 `cindex_final.csv` 等 |
| `OmicsScoresMerge/OmicsScoresMerge.ipynb` | 读取多种子（如 1234–1238）的 OmicsNet 预测分数，按 eid 平均，得到 Final 分数文件（train/val/internal_test/external_test），供 Cox 与可视化使用 |
| `SHAPImportanceMerge/SHAPImportanceMerge.ipynb` | 将多轮运行的 SHAP 值按样本/特征取平均，得到 MetNet/ProNet 的最终 SHAP 表，供热图与 beeswarm 等解释性图表使用 |
| `SurvivalAnalysis/SurvivalAnalysis.ipynb` | 加载外部测试集临床数据与 PRS/metscore/proscore，做 Kaplan-Meier、Cox、logrank 等生存分析，输出分位数分组与风险比等，为森林图与生存曲线提供数据 |

---

## 八、可视化 `visualization/`

| 文件 | 功能 |
|------|------|
| `visualization.Rmd` | **主控 R Markdown**：读取 `saved/results/` 下的汇总表与 C-index/SHAP 等，生成论文用图表与表格。包含：GWAS 文件格式转换、基线特征表、随访与发病描述、Kaplan-Meier 与 HR 森林图、C-index/ΔC-index 森林图（主分析、地理验证、单 vs 多组学、亚组与敏感性）、SHAP 热图与 beeswarm 等。输出到 `results/Figure/` 与 `results/Table/`（路径需在 Rmd 内配置）。 |

---

## 九、典型数据与输出目录（README 约定）

- **输入**：`data/processed/`、`data/split_seed-XXX/`、`data/gwas_summary_statistics/`
- **输出**：
  - `saved/log/`：数据准备、训练、分析的日志
  - `saved/models/`：训练好的 OmicsNet 等模型权重
  - `saved/results/`：预测分数（Scores）、C-index（Cindex）、SHAP（SHAP）、发病概率（Incident_Probabilities）、生存分析结果等
- **图表**：由 `visualization.Rmd` 写入 `results/Figure/`、`results/Table/`（路径以 Rmd 中设置为准）

---

## 十、流程概览（从数据到图表）

1. **数据准备**：`OmicsGenerator` + `OutcomesGenerator` + `PredictorsGenerator` → 原始与 processed 数据；`data_preparation.py` → 划分与插补 → `split_seed-XXX` 的 X/y/e。
2. **组学模型**：`hyperparameter_tuning` 得到 `final_param.json` → `training/*_5times.py` 训练 MetNet/ProNet → 分数与 SHAP 写入 `saved/results/`。
3. **基线与传统模型**：`baseline/` 计算临床评分与 ML；`coxph/` 与 `fine_gray/` 在多种 predictor 组合上做 Cox/Fine-Gray 与 Bootstrap。
4. **结果整合**：`analysis/*.ipynb` 合并 C-index、分数、SHAP、生存分析结果。
5. **论文图表**：`visualization.Rmd` 读取 `saved/results/` 与合并结果，生成 Figure/Table。

---

## 十一、文件索引（按目录）

```
config/
  Metabolomics/final_param.json
  Metabolomics_no_statins/final_param.json
  Proteomics/final_param.json

data/
  gwas_summary_statistics/
    donwload_summary_statistics.sh
    PGS*.txt, PGS*_formatted.txt

data_loader/
  data_loader.py
  data_preparation.py
  data_preparation_raw.py
  DataPreparation/
    OmicsGenerator/
      NMRDataProcessing/S0_ReadRawNMRData.r, S1_MetabolomicsQC.py
      OlinkDataProcessing/S0_ReadRawOlinkData.py, S1_ProteomicsQC.py
      PRSProcessing/S0_..., S1_..., S2_merge_prs.py
    OutcomesGenerator/S0_OutcomesInfo.py, S1_OutcomesGenerator.py
    PredictorsGenerator/
      Biofluids, Demographic, DiseaseHistory, FamilyHistory, GenotypeQC,
      Lifestyle, LifestyleDiet, LifestyleSocial, MedicationHistory, PhysicalMeasurements

model/
  model.py
  trainer.py
  baseline/clinical_scores.py, machine_learning_metabolomics.py, machine_learning_proteomics.py
  hyperparameter_tuning/metabolomics.py, proteomics.py
  training/metabolomics_5times.py, metabolomics_no_statins_5times.py, proteomics_5times.py
  coxph/*.py
  fine_gray/*.R

analysis/
  MetricMerge/MetricMerge.ipynb
  OmicsScoresMerge/OmicsScoresMerge.ipynb
  SHAPImportanceMerge/SHAPImportanceMerge.ipynb
  SurvivalAnalysis/SurvivalAnalysis.ipynb

visualization/
  visualization.Rmd

src/
  comics_*.pdf, cover_with_authors.png
```

以上为 CardiOmicScore 项目结构与各文件功能的说明文档。使用前请根据 README 配置 Python/R 环境与路径（如 `/your path/cardiomicscore`）。
