# ================================
# 0. 设置工作环境
# ================================
rm(list = ls())
set.seed(1615)  # 设置随机数种子以保证结果可重复
library(dplyr)        # 数据处理
library(ggplot2)      # 数据可视化
library(tidyr)        # 数据重塑
library(pROC)         # ROC曲线分析
library(splines)      # 限制性样条函数
library(segmented)    # 分段回归
library(survival)     # 生存分析
library(glmnet)       # LASSO回归
library(rms)          # 校准曲线
library(vcd)          # 可视化关联
library(corrplot)     # 相关性图
library(igraph)       # 网络图
library(patchwork)    # 组合图表
library(mgcv)         # 广义加性模型
library(boot)         # Bootstrap分析
library(pheatmap)     # 热图
library(survminer)    # 生存分析可视化
library(caret)        # 机器学习工具
library(SHAPforxgboost)  # SHAP值计算

# ================================
# 1. 创建模拟数据
# ================================

# 1.1 患者基本信息数据
n_patients <- 150  # 样本量

# 基本人口统计学和临床特征
patients <- data.frame(
  patient_id = 1:n_patients,
  age = round(rnorm(n_patients, mean = 63.8, sd = 15.2)),
  gender = sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.573, 0.427)),
  weight = round(rnorm(n_patients, mean = 67.5, sd = 12.3), 1),
  SOFA_score = round(rnorm(n_patients, mean = 7.8, sd = 3.2)),
  APACHE_II = round(rnorm(n_patients, mean = 18.5, sd = 6.8)),
  hypertension = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.57, 0.43)),
  diabetes = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.693, 0.307)),
  cardiac_disease = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.787, 0.213)),
  stroke = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.867, 0.133)),
  infection_site = sample(c("Pulmonary", "Abdominal", "Urinary", "Bloodstream", "Other"), 
                          n_patients, replace = TRUE, 
                          prob = c(0.453, 0.213, 0.167, 0.107, 0.06)),
  septic_shock = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.75, 0.25)),
  mechanical_ventilation = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.653, 0.347)),
  vasopressors = sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.68, 0.32))
)

# 添加更多合并症和总合并症数
patients$renal_disease <- sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.88, 0.12))
patients$liver_disease <- sample(c(0, 1), n_patients, replace = TRUE, prob = c(0.9, 0.1))
patients$comorbidity_count <- patients$hypertension + patients$diabetes + 
patients$cardiac_disease + patients$stroke + patients$renal_disease + patients$liver_disease

# 1.2 抗菌药物使用数据 (更加详细)
# 定义抗生素类别和具体药物
antibiotic_classes <- list(
  "Penicillins" = c("Piperacillin", "Ampicillin", "Amoxicillin"),
  "Cephalosporins" = c("Ceftazidime", "Cefepime", "Ceftriaxone", "Cefuroxime"),
  "Carbapenems" = c("Imipenem", "Meropenem"),
  "Aminoglycosides" = c("Gentamicin", "Amikacin"),
  "Quinolones" = c("Ciprofloxacin", "Levofloxacin"),
  "Glycopeptides" = c("Vancomycin", "Teicoplanin"),
  "Others" = c("Linezolid", "Piperacillin/Tazobactam", "Trimethoprim/Sulfamethoxazole")
)

# 抗生素使用原因和等级
reasons <- c("Empirical", "Pathogen-directed", "Prophylactic")
levels <- c("Ordinary", "Restricted", "Special")

# 创建抗生素使用数据框 - 每位患者可能使用多种抗生素
antibiotics <- data.frame()
antibiotic_count <- numeric(n_patients)  # 存储每位患者的抗生素种类数

for (i in 1:n_patients) {
  # 确定该患者使用的抗生素种类数
  n_antibiotics <- sample(1:7, 1, prob = c(0.187, 0.233, 0.28, 0.167, 0.08, 0.033, 0.02))
  antibiotic_count[i] <- n_antibiotics
  
  # 为更严重的患者增加抗生素使用数量的概率
  if (patients$SOFA_score[i] > 10 || patients$septic_shock[i] == 1) {
    n_antibiotics <- max(n_antibiotics, sample(2:5, 1))
  }
  
  for (j in 1:n_antibiotics) {
    # 选择抗生素类别和具体药物
    class_index <- sample(1:length(antibiotic_classes), 1)
    class_name <- names(antibiotic_classes)[class_index]
    drug_name <- sample(antibiotic_classes[[class_index]], 1)
    
    # 为每种抗生素添加使用特征
    antibiotics <- rbind(antibiotics, data.frame(
      patient_id = i,
      antibiotic_class = class_name,
      antibiotic_name = drug_name,
      start_time = round(runif(1, 0, 24), 1),  # 开始时间(小时)
      duration = round(runif(1, 3, 14), 1),    # 持续时间(天)
      daily_dose = round(runif(1, 1, 4), 1),   # 日剂量(g)
      reason = sample(reasons, 1),             # 使用原因
      level = sample(levels, 1, prob = c(0.62, 0.28, 0.1)),  # 抗生素等级
      spectrum = sample(c("Narrow", "Medium", "Broad"), 1, 
                        prob = c(0.2, 0.5, 0.3)),  # 抗菌谱
      stringsAsFactors = FALSE
    ))
  }
}

# 更新患者数据框中的抗生素种类数
patients$antibiotic_count <- antibiotic_count

# 特殊级抗生素使用标记
patients$special_antibiotic <- sapply(patients$patient_id, function(pid) {
  any(antibiotics$patient_id == pid & antibiotics$level == "Special")
})

# 1.3 联合用药模式和治疗调整
# 确定每位患者的联合用药模式
determine_combination_mode <- function(pid) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  n_abs <- nrow(patient_abs)
  
  if (n_abs == 1) return("Single")
  
  # 按开始时间排序
  patient_abs <- patient_abs[order(patient_abs$start_time), ]
  
  # 检查是否存在时间重叠
  overlap <- FALSE
  for(i in 1:(n_abs-1)) {
    if((patient_abs$start_time[i] + patient_abs$duration[i]) > 
       patient_abs$start_time[i+1]) {
      overlap <- TRUE
      break
    }
  }
  
  # 判断模式
  if(!overlap) return("Sequential")
  else if(all(patient_abs$start_time < 2)) return("Parallel")
  else return("Mixed")
}

patients$combination_mode <- sapply(patients$patient_id, determine_combination_mode)

# 治疗调整模式
determine_adjustment <- function(pid) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  
  # 根据SOFA评分和临床表现决定调整模式
  if(patients$SOFA_score[patients$patient_id == pid] > 10 || 
     patients$septic_shock[patients$patient_id == pid] == 1) {
    return(sample(c("No adjustment", "Early de-escalation", "Late de-escalation", "Escalation"), 
                  1, prob = c(0.1, 0.1, 0.2, 0.6)))
  } else {
    return(sample(c("No adjustment", "Early de-escalation", "Late de-escalation", "Escalation"), 
                  1, prob = c(0.4, 0.3, 0.1, 0.2)))
  }
}

patients$adjustment_mode <- sapply(patients$patient_id, determine_adjustment)

# 1.4 微生物学和耐药性数据
# 创建微生物学检测结果
microbiology <- data.frame()

# 定义常见病原体
pathogens <- c("Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa", 
               "Staphylococcus aureus", "Enterococcus faecalis", "Acinetobacter baumannii", 
               "Candida albicans")

# 定义耐药表型
resistance_types <- c("None", "ESBL", "CRE", "MRSA", "VRE", "MDR-Acinetobacter")

for (i in 1:n_patients) {
  # 决定是否有病原体检出
  if (runif(1) < 0.8) {  # 80%的患者有病原体检出
    n_pathogens <- sample(1:3, 1, prob = c(0.6, 0.3, 0.1))  # 可能检出多个病原体
    
    for (j in 1:n_pathogens) {
      pathogen <- sample(pathogens, 1)
      
      # 根据抗生素使用情况决定是否出现耐药
      resistance_prob <- 0.2  # 基础耐药概率
      
      # 抗生素种类数越多，耐药概率越高
      if (patients$antibiotic_count[i] > 3) {
        resistance_prob <- resistance_prob + 0.1 * (patients$antibiotic_count[i] - 3)
      }
      
      # 特殊级抗生素使用增加耐药概率
      if (patients$special_antibiotic[i]) {
        resistance_prob <- resistance_prob + 0.15
      }
      
      # 根据病原体类型选择可能的耐药表型
      if (pathogen == "Escherichia coli" || pathogen == "Klebsiella pneumoniae") {
        resistance <- sample(c("None", "ESBL", "CRE"), 1, 
                             prob = c(1 - resistance_prob, resistance_prob * 0.7, resistance_prob * 0.3))
      } else if (pathogen == "Staphylococcus aureus") {
        resistance <- sample(c("None", "MRSA"), 1, 
                             prob = c(1 - resistance_prob, resistance_prob))
      } else if (pathogen == "Enterococcus faecalis") {
        resistance <- sample(c("None", "VRE"), 1, 
                             prob = c(1 - resistance_prob, resistance_prob))
      } else if (pathogen == "Acinetobacter baumannii") {
        resistance <- sample(c("None", "MDR-Acinetobacter"), 1, 
                             prob = c(1 - resistance_prob, resistance_prob))
      } else {
        resistance <- "None"
      }
      
      # 添加到微生物学数据框
      microbiology <- rbind(microbiology, data.frame(
        patient_id = i,
        pathogen = pathogen,
        resistance = resistance,
        culture_time = round(runif(1, 24, 120), 1),  # 培养时间(小时)
        specimen_source = sample(c("Blood", "Sputum", "Urine", "Wound", "Other"), 1),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 为患者添加耐药标记
patients$resistance_detected <- sapply(patients$patient_id, function(pid) {
  any(microbiology$patient_id == pid & microbiology$resistance != "None")
})

# MDR分类（多重耐药、广泛耐药、全耐药）
determine_mdr_class <- function(pid) {
  patient_micro <- microbiology[microbiology$patient_id == pid, ]
  
  if (nrow(patient_micro) == 0 || all(patient_micro$resistance == "None")) {
    return("No resistant")
  }
  
  # 检查耐药类型和数量
  resistance_count <- sum(patient_micro$resistance != "None")
  
  if (resistance_count == 1) return("MDR")
  else if (resistance_count == 2) return("XDR")
  else return("PDR")
}

patients$mdr_class <- sapply(patients$patient_id, determine_mdr_class)

# 1.5 临床结局数据
# 28天死亡率预测函数
predict_mortality <- function(patient) {
  # 基础概率
  prob <- 0.05
  
  # 年龄影响
  if (patient$age > 70) prob <- prob + 0.1
  
  # SOFA评分影响
  prob <- prob + 0.01 * patient$SOFA_score^2
  
  # 药物和治疗因素
  if (patient$antibiotic_count > 4) prob <- prob + 0.05
  if (patient$special_antibiotic) prob <- prob + 0.1
  if (patient$adjustment_mode == "Escalation") prob <- prob + 0.1
  if (patient$combination_mode == "Mixed") prob <- prob + 0.05
  
  # 耐药因素
  if (patient$resistance_detected) {
    if (patient$mdr_class == "MDR") prob <- prob + 0.1
    else if (patient$mdr_class == "XDR") prob <- prob + 0.15
    else if (patient$mdr_class == "PDR") prob <- prob + 0.2
  }
  
  # 器官支持
  if (patient$mechanical_ventilation) prob <- prob + 0.15
  if (patient$vasopressors) prob <- prob + 0.15
  
  # 感染部位
  if (patient$infection_site == "Pulmonary") prob <- prob + 0.05
  if (patient$infection_site == "Bloodstream") prob <- prob + 0.1
  
  return(min(prob, 0.95))  # 限制最高死亡概率为95%
}

# 计算死亡概率并分配结局
patients$mortality_prob <- sapply(1:nrow(patients), function(i) predict_mortality(patients[i, ]))

# 确保死亡率约为10%
target_deaths <- round(0.10 * nrow(patients))  # 目标10%死亡率
ordered_patients <- patients[order(-patients$mortality_prob), ]
ordered_patients$death_28d <- 0
ordered_patients$death_28d[1:target_deaths] <- 1
shuffled_indices <- sample(1:nrow(ordered_patients))
ordered_patients <- ordered_patients[shuffled_indices, ]
patients$death_28d <- ordered_patients$death_28d[match(1:nrow(patients), shuffled_indices)]

# 验证结果
cat("生成的28天死亡率:", mean(patients$death_28d) * 100, "%\n")

# 其他临床结局
patients$icu_los <- round(rexp(n_patients, 1/7) * (1 + patients$SOFA_score/10), 1)  # ICU住院日
patients$hospital_los <- round(patients$icu_los + rexp(n_patients, 1/10), 1)  # 总住院日
patients$mechanical_ventilation_days <- ifelse(patients$mechanical_ventilation == 1, 
                                               round(rexp(n_patients, 1/5), 1), 0)  # 机械通气日
patients$antibiotic_treatment_days <- round(runif(n_patients, 5, 21), 1)  # 抗生素治疗日

# 临床反应 (72小时)
patients$clinical_response_72h <- sapply(1:nrow(patients), function(i) {
  if (patients$death_28d[i] == 1 || patients$SOFA_score[i] > 12) {
    sample(c("Improved", "Stable", "Deteriorated"), 1, prob = c(0.1, 0.3, 0.6))
  } else {
    sample(c("Improved", "Stable", "Deteriorated"), 1, prob = c(0.6, 0.3, 0.1))
  }
})

# 治疗失败 (死亡或7天内无临床改善)
patients$treatment_failure <- ifelse(
  patients$death_28d == 1 | patients$clinical_response_72h == "Deteriorated" | 
    (patients$clinical_response_72h != "Improved" & runif(n_patients) > 0.7),
  1, 0
)

# 1.6 创建Shannon多样性指数
# 为每位患者计算Shannon多样性指数
calculate_shannon_index <- function(pid) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  
  if (nrow(patient_abs) <= 1) return(0)  # 单一抗生素无多样性
  
  # 计算每种抗生素的持续时间比例
  total_days <- sum(patient_abs$duration)
  proportions <- patient_abs$duration / total_days
  
  # Shannon指数计算
  -sum(proportions * log(proportions))
}

patients$shannon_index <- sapply(patients$patient_id, calculate_shannon_index)

# 确保数据合理性
summary(patients)

# ================================
# 2. 数据预处理与亚组划分
# ================================

# 2.1 创建亚组
# 基于年龄分组（老年vs非老年）
patients$age_group <- ifelse(patients$age >= 65, "Elderly (≥65)", "Non-elderly (<65)")
table(patients$age_group)

# 基于感染部位分组
table(patients$infection_site)

# 基于SOFA评分分组（重症vs轻中症）
patients$severity_group <- ifelse(patients$SOFA_score >= 10, "Severe (SOFA≥10)", "Non-severe (SOFA<10)")
table(patients$severity_group)

# 基于休克状态分组
table(patients$septic_shock)
patients$shock_group <- ifelse(patients$septic_shock == 1, "Septic shock", "No shock")

# 基于耐药状态分组
table(patients$resistance_detected)

# 2.2 数据清洗与缺失值处理
# 检查缺失值
missing_values <- colSums(is.na(patients))
print(missing_values[missing_values > 0])

# 对缺失数据进行填充或剔除
patients$shannon_index[is.na(patients$shannon_index)] <- 0  # 单一抗生素的患者

# 2.3 计算关键变量间的相关性
numeric_vars <- c("age", "SOFA_score", "antibiotic_count", "shannon_index", "comorbidity_count")
cor_matrix <- cor(patients[, numeric_vars], use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         title = "关键变量相关性矩阵")

# ================================
# 3. 描述性统计分析
# ================================

# 3.1 基本特征分析
cat("患者基本特征分析\n")
cat("==================\n")
cat("总患者数:", nrow(patients), "\n")
cat("年龄(均值±标准差):", mean(patients$age), "±", sd(patients$age), "\n")
cat("性别(男性比例):", sum(patients$gender == "Male") / nrow(patients), "\n")
cat("SOFA评分(均值±标准差):", mean(patients$SOFA_score), "±", sd(patients$SOFA_score), "\n")
cat("28天死亡率:", mean(patients$death_28d), "\n")

# 按生存状态分组比较基线特征
patients_by_survival <- split(patients, patients$death_28d)
survivors <- patients_by_survival$`0`
non_survivors <- patients_by_survival$`1`

# 创建比较表格函数
compare_groups <- function(var_name, var1, var2, continuous = TRUE) {
  if (continuous) {
    # 连续变量
    test_result <- t.test(var1, var2)
    p_value <- round(test_result$p.value, 4)
    data.frame(
      Variable = var_name,
      Survivors = paste0(round(mean(var1), 2), " ± ", round(sd(var1), 2)),
      Non_survivors = paste0(round(mean(var2), 2), " ± ", round(sd(var2), 2)),
      P_value = p_value,
      stringsAsFactors = FALSE
    )
  } else {
    # 分类变量
    test_result <- prop.test(c(sum(var1), sum(var2)), c(length(var1), length(var2)))
    p_value <- round(test_result$p.value, 4)
    data.frame(
      Variable = var_name,
      Survivors = paste0(sum(var1), " (", round(100 * mean(var1), 1), "%)"),
      Non_survivors = paste0(sum(var2), " (", round(100 * mean(var2), 1), "%)"),
      P_value = p_value,
      stringsAsFactors = FALSE
    )
  }
}

# 构建比较表格
comparison_table <- rbind(
  compare_groups("年龄", survivors$age, non_survivors$age),
  compare_groups("SOFA评分", survivors$SOFA_score, non_survivors$SOFA_score),
  compare_groups("抗生素种类数", survivors$antibiotic_count, non_survivors$antibiotic_count),
  compare_groups("特殊级抗生素使用", survivors$special_antibiotic, non_survivors$special_antibiotic, FALSE),
  compare_groups("耐药菌检出", survivors$resistance_detected, non_survivors$resistance_detected, FALSE),
  compare_groups("机械通气", survivors$mechanical_ventilation, non_survivors$mechanical_ventilation, FALSE),
  compare_groups("血管活性药物使用", survivors$vasopressors, non_survivors$vasopressors, FALSE)
)

print(comparison_table)

# 3.2 亚组描述性统计
# 创建一个函数进行亚组描述性统计
subgroup_stats <- function(data, group_var, outcome_var = "death_28d") {
  # 确保组变量是因子
  data[[group_var]] <- as.factor(data[[group_var]])
  
  # 计算每组的基本统计量
  result <- data %>%
    group_by(.data[[group_var]]) %>%
    dplyr::summarize(
      N = n(),
      Death_rate = mean(.data[[outcome_var]]),
      Resistance_rate = mean(resistance_detected),
      Mean_antibiotic_count = mean(antibiotic_count),
      Special_antibiotic_use = mean(special_antibiotic),
      Mean_SOFA = mean(SOFA_score)
    )
  
  return(result)
}

# 各亚组的基本特征
age_group_stats <- subgroup_stats(patients, "age_group")
infection_site_stats <- subgroup_stats(patients, "infection_site")
severity_group_stats <- subgroup_stats(patients, "severity_group")
shock_group_stats <- subgroup_stats(patients, "shock_group")

print("年龄亚组基本特征:")
print(age_group_stats)
print("感染部位亚组基本特征:")
print(infection_site_stats)
print("疾病严重程度亚组基本特征:")
print(severity_group_stats)
print("休克状态亚组基本特征:")
print(shock_group_stats)

# 3.3 可视化亚组特征
# 年龄亚组的抗生素种类数与死亡率
p_age_group <- ggplot(patients, aes(x = antibiotic_count, fill = factor(death_28d))) +
  geom_bar(position = "fill") +
  facet_wrap(~age_group) +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("Survivors", "Non-survivors")) +
  labs(title = "年龄亚组中抗生素种类数与死亡率的关系",
       x = "抗生素种类数", y = "比例", fill = "死亡") +
  theme_minimal()

# 疾病严重程度亚组的抗生素种类数与死亡率
p_severity_group <- ggplot(patients, aes(x = antibiotic_count, fill = factor(death_28d))) +
  geom_bar(position = "fill") +
  facet_wrap(~severity_group) +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("Survivors", "Non-survivors")) +
  labs(title = "疾病严重程度亚组中抗生素种类数与死亡率的关系",
       x = "抗生素种类数", y = "比例", fill = "死亡") +
  theme_minimal()

# 合并亚组可视化
print(p_age_group)
print(p_severity_group)

# 3.4 各亚组相关性比较
# 计算各亚组中抗生素种类数与耐药风险的相关性
cor_by_subgroup <- function(data, group_var, var1 = "antibiotic_count", var2 = "resistance_detected") {
  groups <- unique(data[[group_var]])
  result <- data.frame(Group = character(), Correlation = numeric(), P_value = numeric())
  
  for(group in groups) {
    sub_data <- data[data[[group_var]] == group, ]
    
    # 计算相关性(点二列相关)
    if(is.numeric(sub_data[[var1]]) && is.logical(sub_data[[var2]])) {
      cor_test <- cor.test(sub_data[[var1]], as.numeric(sub_data[[var2]]))
      result <- rbind(result, data.frame(
        Group = as.character(group), 
        Correlation = round(cor_test$estimate, 3),
        P_value = round(cor_test$p.value, 4)
      ))
    }
  }
  
  return(result)
}

# 计算各亚组相关性
cor_age_group <- cor_by_subgroup(patients, "age_group")
cor_severity_group <- cor_by_subgroup(patients, "severity_group")
cor_shock_group <- cor_by_subgroup(patients, "shock_group")

print("年龄亚组中抗生素种类数与耐药风险的相关性:")
print(cor_age_group)
print("疾病严重程度亚组中抗生素种类数与耐药风险的相关性:")
print(cor_severity_group)
print("休克状态亚组中抗生素种类数与耐药风险的相关性:")
print(cor_shock_group)

# ================================
# 4. 抗生素多样性与耐药风险分析
# ================================

# 4.1 全人群分析
# 抗生素种类数与耐药风险的关系
resistance_by_count <- patients %>%
  group_by(antibiotic_count) %>%
  dplyr::summarize(
    n_patients = n(),
    n_resistant = sum(resistance_detected),
    resistance_rate = mean(resistance_detected)
  )

print(resistance_by_count)

# 绘制抗生素种类数与耐药风险的关系图
p1 <- ggplot(resistance_by_count, aes(x = antibiotic_count, y = resistance_rate)) +
  geom_point(aes(size = n_patients), color = "darkblue", alpha = 0.7) +
  geom_line() +
  labs(title = "全人群: 抗生素种类数与耐药风险的关系", 
       x = "抗生素种类数", 
       y = "耐药菌检出率",
       size = "患者数") +
  theme_minimal() +
  ylim(0, 1)

# 4.2 亚组分析
# 针对多个亚组进行分析
subgroups <- list(
  Age = c("age_group", "年龄亚组"),
  Severity = c("severity_group", "疾病严重程度亚组"),
  Shock = c("shock_group", "休克状态亚组")
)

# 创建各亚组的抗生素种类数与耐药风险关系图
subgroup_plots <- list()
subgroup_models <- list()
subgroup_thresholds <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  sg_name <- subgroups[[sg]][2]
  
  # 计算亚组中抗生素种类数与耐药风险的关系
  sg_resistance_by_count <- patients %>%
    group_by(.data[[sg_var]], antibiotic_count) %>%
    dplyr::summarize(
      n_patients = n(),
      n_resistant = sum(resistance_detected),
      resistance_rate = mean(resistance_detected),
      .groups = "drop"
    )
  
  # 绘制亚组关系图
  p <- ggplot(sg_resistance_by_count, aes(x = antibiotic_count, y = resistance_rate, color = .data[[sg_var]])) +
    geom_point(aes(size = n_patients), alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = paste0(sg_name, ": 抗生素种类数与耐药风险的关系"), 
         x = "抗生素种类数", 
         y = "耐药菌检出率",
         color = sg_name,
         size = "患者数") +
    theme_minimal() +
    ylim(0, 1)
  
  subgroup_plots[[sg]] <- p
  
  # 亚组Logistic回归分析
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_model_name <- paste0(sg, "_", group)
    
    # Logistic回归分析
    sg_model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score, 
                    family = binomial, data = sg_data)
    subgroup_models[[sg_model_name]] <- sg_model
    
    # 尝试进行分段回归分析
    tryCatch({
      sg_segmented <- segmented(sg_model, seg.Z = ~ antibiotic_count, psi = 5)
      sg_threshold <- sg_segmented$psi[2]
      cat(sg_name, "亚组", group, "的抗生素种类数阈值点:", sg_threshold, "\n")
      subgroup_thresholds[[sg_model_name]] <- sg_threshold
    }, error = function(e) {
      cat(sg_name, "亚组", group, "的分段回归分析失败:", e$message, "\n")
      subgroup_thresholds[[sg_model_name]] <- NA
    })
  }
}

# 显示亚组图表
for (sg in names(subgroup_plots)) {
  print(subgroup_plots[[sg]])
}

# 4.3 汇总各亚组的阈值点和风险差异
threshold_summary <- data.frame(
  Subgroup = character(),
  Group = character(),
  Threshold = numeric(),
  OR_below_threshold = numeric(),
  OR_above_threshold = numeric(),
  OR_ratio = numeric(),
  stringsAsFactors = FALSE
)

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_model_name <- paste0(sg, "_", group)
    
    # 提取阈值
    threshold <- subgroup_thresholds[[sg_model_name]]
    
    if (!is.na(threshold)) {
      # 筛选亚组数据
      sg_data <- patients[patients[[sg_var]] == group, ]
      
      # 计算阈值前后风险比
      sg_data$above_threshold <- sg_data$antibiotic_count > threshold
      
      # 阈值前后单独建模
      model_below <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score, 
                         family = binomial, data = sg_data[!sg_data$above_threshold, ])
      
      model_above <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score, 
                         family = binomial, data = sg_data[sg_data$above_threshold, ])
      
      # 提取风险比
      or_below <- exp(coef(model_below)["antibiotic_count"])
      or_above <- exp(coef(model_above)["antibiotic_count"])
      or_ratio <- or_above / or_below
      
      # 添加到汇总表
      threshold_summary <- rbind(threshold_summary, data.frame(
        Subgroup = sg,
        Group = group,
        Threshold = round(threshold, 2),
        OR_below_threshold = round(or_below, 2),
        OR_above_threshold = round(or_above, 2),
        OR_ratio = round(or_ratio, 2),
        stringsAsFactors = FALSE
      ))
    }
  }
}

print("各亚组抗生素种类数阈值及阈值前后风险比:")
print(threshold_summary)

# 4.4 全人群与亚组的Logistic回归比较
# 全人群Logistic回归
logit_model_all <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score + infection_site, 
                       family = binomial, data = patients)
summary(logit_model_all)

# 创建亚组回归系数比较表
coef_comparison <- data.frame(
  Variable = c("(Intercept)", "antibiotic_count", "age", "SOFA_score"),
  All_population = round(coef(logit_model_all)[1:4], 3),
  stringsAsFactors = FALSE
)

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_model_name <- paste0(sg, "_", group)
    sg_model <- subgroup_models[[sg_model_name]]
    
    if (!is.null(sg_model)) {
      # 提取系数
      sg_coefs <- coef(sg_model)[1:4]
      names(sg_coefs) <- c("(Intercept)", "antibiotic_count", "age", "SOFA_score")
      
      # 添加到比较表
      coef_comparison[[paste0(sg, "_", group)]] <- round(sg_coefs, 3)
    }
  }
}

print("全人群与亚组Logistic回归系数比较:")
print(coef_comparison)

# 4.5 全人群分段回归分析
segmented_model_all <- try(
  segmented(logit_model_all, seg.Z = ~ antibiotic_count, psi = 5),
  silent = TRUE
)

if (!inherits(segmented_model_all, "try-error")) {
  # 提取阈值点
  threshold_all <- segmented_model_all$psi[2]
  cat("全人群抗生素种类数阈值点:", threshold_all, "\n")
  
  # 输出阈值前后的斜率
  slopes_all <- slope(segmented_model_all)
  print("全人群阈值前后斜率:")
  print(slopes_all)
  
  # 预测值曲线 - 分段回归模型
  newdata <- data.frame(
    antibiotic_count = seq(1, 7, by = 0.1),
    age = mean(patients$age),
    SOFA_score = mean(patients$SOFA_score),
    infection_site = "Pulmonary"  # 使用最常见的感染部位
  )
  
  newdata$pred_logit <- predict(logit_model_all, newdata)
  newdata$pred_seg <- predict(segmented_model_all, newdata)
  
  newdata$prob_logit <- exp(newdata$pred_logit) / (1 + exp(newdata$pred_logit))
  newdata$prob_seg <- exp(newdata$pred_seg) / (1 + exp(newdata$pred_seg))
  
  # 同时绘制全人群和亚组的阈值效应图
  p2 <- ggplot(newdata, aes(x = antibiotic_count)) +
    geom_line(aes(y = prob_logit, color = "Linear"), size = 1) +
    geom_line(aes(y = prob_seg, color = "Segmented"), size = 1) +
    geom_point(data = resistance_by_count, aes(y = resistance_rate, size = n_patients), 
               color = "black", alpha = 0.7) +
    geom_vline(xintercept = threshold_all, linetype = "dashed", color = "darkred") +
    annotate("text", x = threshold_all + 0.3, y = 0.8, 
             label = paste("全人群阈值:", round(threshold_all, 2)), color = "darkred") +
    scale_color_manual(values = c("Linear" = "blue", "Segmented" = "red"),
                       name = "模型类型") +
    labs(title = "全人群抗生素种类数与耐药风险的非线性关系", 
         subtitle = "分段线性回归显示阈值效应",
         x = "抗生素种类数", 
         y = "预测耐药风险",
         size = "患者数") +
    theme_minimal()
  
  print(p2)
  
  # 4.6 Bootstrap抽样分析阈值点的稳定性
  threshold_boot <- function(data, indices) {
    d <- data[indices, ]
    model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score, 
                 family = binomial, data = d)
    
    # 尝试获取阈值点，如果模型收敛失败则返回NA
    tryCatch({
      seg_model <- segmented(model, seg.Z = ~ antibiotic_count, psi = list(antibiotic_count = c(threshold_all)))
      return(seg_model$psi[2])
    }, error = function(e) {
      return(NA)
    })
  }
  
  # 执行Bootstrap
  set.seed(1615)  # 设置随机种子以确保结果可重复
  boot_results <- boot(data = patients, statistic = threshold_boot, R = 500)
  
  # 计算Bootstrap置信区间
  boot_ci <- boot.ci(boot_results, type = "basic", conf = 0.95)
  cat("全人群阈值点的Bootstrap 95%置信区间：", boot_ci$basic[4], "到", boot_ci$basic[5], "\n")
  
  # 绘制Bootstrap结果的直方图
  valid_results <- boot_results$t[!is.na(boot_results$t)]
  hist(valid_results, breaks = 30, main = "全人群阈值点Bootstrap分布", 
       xlab = "阈值点估计值", ylab = "频率", col = "lightblue", border = "white")
  abline(v = threshold_all, col = "red", lwd = 2)
  abline(v = boot_ci$basic[4], col = "blue", lty = 2)
  abline(v = boot_ci$basic[5], col = "blue", lty = 2)
  legend("topright", legend = c("原始估计", "95%置信区间"), 
         col = c("red", "blue"), lty = c(1, 2), lwd = c(2, 1))
}

# 4.7 亚组与全人群阈值比较
# 创建综合阈值比较表
threshold_all_value <- ifelse(exists("threshold_all"), threshold_all, NA)

threshold_comparison <- data.frame(
  Group = c("全人群", threshold_summary$Group),
  Threshold = c(threshold_all_value, threshold_summary$Threshold),
  stringsAsFactors = FALSE
)

# 可视化阈值比较
p3 <- ggplot(threshold_comparison, aes(x = reorder(Group, Threshold), y = Threshold)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(Threshold, 2)), vjust = -0.5) +
  labs(title = "全人群与亚组抗生素种类数阈值比较",
       x = "分组", y = "阈值点") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# ================================
# 5. 抗生素使用时序特征分析
# ================================

# 5.1 创建抗生素使用时间轴数据
# 增加时间特征
antibiotics$end_time <- antibiotics$start_time + antibiotics$duration * 24  # 结束时间(小时)
antibiotics$time_to_first_antibiotic <- antibiotics$start_time  # 从入院到首次使用的时间

# 计算每位患者的抗生素时序特征
patient_time_features <- data.frame()
for(pid in unique(patients$patient_id)) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  
  # 按开始时间排序
  patient_abs <- patient_abs[order(patient_abs$start_time), ]
  
  # 计算时序特征
  if(nrow(patient_abs) > 0) {
    # 首次抗生素使用时间
    first_use_time <- min(patient_abs$start_time)
    
    # 总暴露时间(考虑重叠)
    exposure_timeline <- numeric(24*30)  # 创建30天的时间线（每小时一个点）
    
    for(i in 1:nrow(patient_abs)) {
      start_hour <- ceiling(patient_abs$start_time[i])
      end_hour <- floor(patient_abs$end_time[i])
      if(start_hour <= end_hour && start_hour <= length(exposure_timeline) && end_hour >= 1) {
        start_idx <- max(1, start_hour)
        end_idx <- min(length(exposure_timeline), end_hour)
        exposure_timeline[start_idx:end_idx] <- 1
      }
    }
    
    total_exposure_hours <- sum(exposure_timeline)
    
    # 药物变更频率 (计算不同抗生素开始使用的时间点数)
    change_points <- unique(patient_abs$start_time)
    change_frequency <- length(change_points)
    
    # 同时使用的最大抗生素数量
    max_concurrent <- 0
    for(hour in 1:length(exposure_timeline)) {
      concurrent_at_hour <- sum(patient_abs$start_time <= hour & patient_abs$end_time >= hour)
      max_concurrent <- max(max_concurrent, concurrent_at_hour)
    }
    
    # 抗生素之间的时间间隔
    time_gaps <- numeric(0)
    if(nrow(patient_abs) > 1) {
      for(i in 1:(nrow(patient_abs)-1)) {
        if(patient_abs$end_time[i] < patient_abs$start_time[i+1]) {
          time_gaps <- c(time_gaps, patient_abs$start_time[i+1] - patient_abs$end_time[i])
        }
      }
    }
    mean_gap <- ifelse(length(time_gaps) > 0, mean(time_gaps), NA)
    
    # 抗生素持续时间的一致性
    duration_consistency <- sd(patient_abs$duration) / mean(patient_abs$duration)
    duration_consistency <- ifelse(is.nan(duration_consistency), 0, duration_consistency)
    
    # 添加到患者时序特征数据框
    patient_time_features <- rbind(patient_time_features, data.frame(
      patient_id = pid,
      first_use_time = first_use_time,
      total_exposure_hours = total_exposure_hours,
      change_frequency = change_frequency,
      max_concurrent = max_concurrent,
      mean_gap = mean_gap,
      duration_consistency = duration_consistency,
      stringsAsFactors = FALSE
    ))
  }
}

# 合并到患者数据
patients <- merge(patients, patient_time_features, by = "patient_id", all.x = TRUE)

# 5.2 时序特征与耐药风险的关系
# 全人群时间特征与耐药风险的关系
time_vs_resistance_all <- glm(
  resistance_detected ~ first_use_time + total_exposure_hours + change_frequency + max_concurrent + 
    # mean_gap +
    age + SOFA_score,
  family = binomial, data = patients
)
summary(time_vs_resistance_all)

# 计算风险比和置信区间
time_or_all <- exp(coef(time_vs_resistance_all))
time_ci_all <- exp(confint(time_vs_resistance_all))
time_results_all <- data.frame(
  Variable = names(coef(time_vs_resistance_all)),
  OR = time_or_all,
  Lower_CI = time_ci_all[,1],
  Upper_CI = time_ci_all[,2],
  P_value = summary(time_vs_resistance_all)$coefficients[,4]
)
print("全人群时序特征与耐药风险:")
print(time_results_all)

# 5.3 亚组时序特征分析
# 为各亚组构建时间特征模型
time_models_by_subgroup <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_model_name <- paste0("time_", sg, "_", group)
    
    # 构建模型
    sg_time_model <- glm(
      resistance_detected ~ first_use_time + total_exposure_hours + change_frequency + 
        max_concurrent + 
        # mean_gap 
      age + SOFA_score,
      family = binomial, data = sg_data
    )
    
    time_models_by_subgroup[[sg_model_name]] <- sg_time_model
  }
}

# 创建亚组时间特征风险比比较表
time_or_comparison <- data.frame(
  Variable = c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent"),
  All_population = round(time_or_all[c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent")], 2),
  stringsAsFactors = FALSE
)

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_model_name <- paste0("time_", sg, "_", group)
    sg_model <- time_models_by_subgroup[[sg_model_name]]
    
    if (!is.null(sg_model)) {
      # 提取风险比
      sg_or <- exp(coef(sg_model)[c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent")])
      
      # 添加到比较表
      time_or_comparison[[paste0(sg, "_", group)]] <- round(sg_or, 2)
    }
  }
}

print("全人群与亚组时序特征风险比比较:")
print(time_or_comparison)

# 5.4 可视化时序特征与耐药风险的关系
# 首次抗生素使用时间与耐药风险
p4 <- ggplot(patients, aes(x = first_use_time, y = as.numeric(resistance_detected))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  facet_wrap(~severity_group) +  # 按严重程度分组
  labs(title = "首次抗生素使用时间与耐药风险 (按严重程度分组)", 
       x = "首次使用时间(小时)", y = "耐药风险") +
  theme_minimal()

# 抗生素暴露时间与耐药风险
p5 <- ggplot(patients, aes(x = total_exposure_hours, y = as.numeric(resistance_detected))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  facet_wrap(~age_group) +  # 按年龄分组
  labs(title = "抗生素暴露时间与耐药风险 (按年龄分组)", 
       x = "总暴露时间(小时)", y = "耐药风险") +
  theme_minimal()

print(p4)
print(p5)

# 5.5 时序特征与患者生存的关系
# 全人群时间依赖性Cox比例风险模型
cox_time_model_all <- coxph(
  Surv(hospital_los, death_28d) ~ first_use_time + total_exposure_hours + change_frequency + 
    max_concurrent + 
    # mean_gap + 
    age + SOFA_score,
  data = patients
)
summary(cox_time_model_all)

# 5.6 按亚组进行生存分析
cox_models_by_subgroup <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_model_name <- paste0("cox_", sg, "_", group)
    
    # 构建Cox模型
    sg_cox_model <- try(
      coxph(
        Surv(hospital_los, death_28d) ~ first_use_time + total_exposure_hours + change_frequency + 
          max_concurrent + mean_gap + age + SOFA_score,
        data = sg_data
      ), 
      silent = TRUE
    )
    
    if (!inherits(sg_cox_model, "try-error")) {
      cox_models_by_subgroup[[sg_model_name]] <- sg_cox_model
    }
  }
}

# 创建亚组Cox模型风险比比较表
cox_hr_comparison <- data.frame(
  Variable = c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent"),
  All_population = round(exp(coef(cox_time_model_all)[c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent")]), 2),
  stringsAsFactors = FALSE
)

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_model_name <- paste0("cox_", sg, "_", group)
    sg_model <- cox_models_by_subgroup[[sg_model_name]]
    
    if (!is.null(sg_model)) {
      # 提取风险比
      sg_hr <- try(exp(coef(sg_model)[c("first_use_time", "total_exposure_hours", "change_frequency", "max_concurrent")]), silent = TRUE)
      
      if (!inherits(sg_hr, "try-error")) {
        # 添加到比较表
        cox_hr_comparison[[paste0(sg, "_", group)]] <- round(sg_hr, 2)
      }
    }
  }
}

print("全人群与亚组时序特征Cox风险比比较:")
print(cox_hr_comparison)

# 5.7 按首次抗生素使用时间分组的生存曲线
# 按首次抗生素使用时间分组(早/中/晚)
patients$first_use_group <- cut(patients$first_use_time, 
                                breaks = c(0, 2, 6, Inf), 
                                labels = c("早期(<2h)", "中期(2-6h)", "晚期(>6h)"))

# 全人群Kaplan-Meier生存曲线
km_all <- survfit(Surv(hospital_los, death_28d) ~ first_use_group, data = patients)
p6 <- ggsurvplot(
  km_all,
  data = patients,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = "Dark2",
  title = "全人群: 首次抗生素使用时间与生存",
  xlab = "住院日",
  ylab = "生存概率",
  legend.title = "首次使用时间"
)

print(p6)

# 按严重程度亚组的生存曲线
# 重症患者
severe_data <- patients[patients$severity_group == "Severe (SOFA≥10)", ]
km_severe <- survfit(Surv(hospital_los, death_28d) ~ first_use_group, data = severe_data)
p7 <- ggsurvplot(
  km_severe,
  data = severe_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = "Dark2",
  title = "重症亚组: 首次抗生素使用时间与生存",
  xlab = "住院日",
  ylab = "生存概率",
  legend.title = "首次使用时间"
)

# 非重症患者
non_severe_data <- patients[patients$severity_group == "Non-severe (SOFA<10)", ]
km_non_severe <- survfit(Surv(hospital_los, death_28d) ~ first_use_group, data = non_severe_data)
p8 <- ggsurvplot(
  km_non_severe,
  data = non_severe_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = "Dark2",
  title = "非重症亚组: 首次抗生素使用时间与生存",
  xlab = "住院日",
  ylab = "生存概率",
  legend.title = "首次使用时间"
)

print(p7)
print(p8)

# ================================
# 6. 抗生素剂量特征与耐药风险分析
# ================================

# 6.1 计算抗生素剂量特征
# WHO DDD标准值(示例)
who_ddd <- data.frame(
  antibiotic_name = c("Piperacillin", "Ampicillin", "Amoxicillin", "Ceftazidime", "Cefepime", "Ceftriaxone", 
                      "Cefuroxime", "Imipenem", "Meropenem", "Gentamicin", "Amikacin", "Ciprofloxacin", 
                      "Levofloxacin", "Vancomycin", "Teicoplanin", "Linezolid", "Piperacillin/Tazobactam", 
                      "Trimethoprim/Sulfamethoxazole"),
  ddd_value = c(14, 6, 3, 4, 2, 2, 3, 2, 2, 0.24, 1, 0.5, 0.5, 2, 0.4, 1.2, 14, 1.92)
)

# 合并DDD值并计算标准化剂量
antibiotics <- merge(antibiotics, who_ddd, by = "antibiotic_name", all.x = TRUE)
antibiotics$ddd_adjusted <- antibiotics$daily_dose / antibiotics$ddd_value

# 计算每位患者的累计DDD和日均DDD
patient_dose_features <- antibiotics %>%
  group_by(patient_id) %>%
  dplyr::summarize(
    total_ddd = sum(ddd_adjusted * duration),
    average_daily_ddd = sum(ddd_adjusted * duration) / max(duration),
    ddd_per_antibiotic = mean(ddd_adjusted),
    ddd_variability = sd(ddd_adjusted) / mean(ddd_adjusted)
  ) %>%
  ungroup()

# 合并到患者数据
patients <- merge(patients, patient_dose_features, by = "patient_id", all.x = TRUE)

# 6.2 全人群剂量特征与耐药风险
# 对DDD进行四分位数分组
patients$total_ddd_group <- cut(patients$total_ddd, 
                                breaks = quantile(patients$total_ddd, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                                labels = c("Q1", "Q2", "Q3", "Q4"),
                                include.lowest = TRUE)

# 计算各组的耐药率
resistance_by_ddd_all <- patients %>%
  group_by(total_ddd_group) %>%
  dplyr::summarize(
    n_patients = n(),
    resistance_rate = mean(resistance_detected)
  )

print("全人群累计DDD与耐药风险:")
print(resistance_by_ddd_all)

# 绘制全人群累计DDD与耐药率的关系
p9 <- ggplot(resistance_by_ddd_all, aes(x = total_ddd_group, y = resistance_rate, fill = total_ddd_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "全人群: 累计DDD与耐药风险", 
       x = "累计DDD四分位组", 
       y = "耐药率") +
  theme_minimal() +
  theme(legend.position = "none")

# 6.3 亚组剂量特征分析
# 计算各亚组的DDD与耐药风险关系
ddd_subgroup_plots <- list()
ddd_models_by_subgroup <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  sg_name <- subgroups[[sg]][2]
  
  # 计算亚组中DDD与耐药风险的关系
  sg_resistance_by_ddd <- patients %>%
    group_by(.data[[sg_var]], total_ddd_group) %>%
    summarize(
      n_patients = n(),
      resistance_rate = mean(resistance_detected),
      .groups = "drop"
    )
  
  # 绘制亚组关系图
  p <- ggplot(sg_resistance_by_ddd, aes(x = total_ddd_group, y = resistance_rate, fill = .data[[sg_var]])) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    labs(title = paste0(sg_name, ": 累计DDD与耐药风险的关系"), 
         x = "累计DDD四分位组", 
         y = "耐药菌检出率",
         fill = sg_name) +
    theme_minimal()
  
  ddd_subgroup_plots[[sg]] <- p
  
  # 亚组剂量特征Logistic回归分析
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_model_name <- paste0("ddd_", sg, "_", group)
    
    # Logistic回归分析
    sg_ddd_model <- glm(
      resistance_detected ~ total_ddd + average_daily_ddd + 
        ddd_per_antibiotic + age + SOFA_score,
      family = binomial, data = sg_data
    )
    ddd_models_by_subgroup[[sg_model_name]] <- sg_ddd_model
  }
}

# 全人群剂量特征的Logistic回归分析
dose_resistance_model_all <- glm(
  resistance_detected ~ total_ddd + average_daily_ddd + 
    ddd_per_antibiotic + age + SOFA_score,
  family = binomial, data = patients
)

summary(dose_resistance_model_all)

# 计算风险比和置信区间
dose_or_all <- exp(coef(dose_resistance_model_all))
dose_ci_all <- exp(confint(dose_resistance_model_all))
dose_results_all <- data.frame(
  Variable = names(coef(dose_resistance_model_all)),
  OR = dose_or_all,
  Lower_CI = dose_ci_all[,1],
  Upper_CI = dose_ci_all[,2],
  P_value = summary(dose_resistance_model_all)$coefficients[,4]
)

print("全人群剂量特征与耐药风险:")
print(dose_results_all)

# 显示亚组图表
for (sg in names(ddd_subgroup_plots)) {
  print(ddd_subgroup_plots[[sg]])
}

# 创建亚组剂量特征风险比比较表
dose_or_comparison <- data.frame(
  Variable = c("total_ddd", "average_daily_ddd", "ddd_per_antibiotic"),
  All_population = round(dose_or_all[c("total_ddd", "average_daily_ddd", "ddd_per_antibiotic")], 2),
  stringsAsFactors = FALSE
)

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_model_name <- paste0("ddd_", sg, "_", group)
    sg_model <- ddd_models_by_subgroup[[sg_model_name]]
    
    if (!is.null(sg_model)) {
      # 提取风险比
      sg_or <- try(exp(coef(sg_model)[c("total_ddd", "average_daily_ddd", "ddd_per_antibiotic")]), silent = TRUE)
      
      if (!inherits(sg_or, "try-error")) {
        # 添加到比较表
        dose_or_comparison[[paste0(sg, "_", group)]] <- round(sg_or, 2)
      }
    }
  }
}

print("全人群与亚组剂量特征风险比比较:")
print(dose_or_comparison)

# 6.4 剂量与临床结局关系的亚组差异
# DDD与治疗效果的关系
p10 <- ggplot(patients, aes(x = total_ddd, y = as.numeric(treatment_failure), color = severity_group)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess") +
  scale_color_manual(values = c("Severe (SOFA≥10)" = "red", "Non-severe (SOFA<10)" = "blue")) +
  labs(title = "累计DDD与治疗失败的关系 (按严重程度分组)", 
       x = "累计DDD", 
       y = "治疗失败(1=是, 0=否)",
       color = "严重程度") +
  theme_minimal()

print(p10)

# DDD与死亡风险的关系
p11 <- ggplot(patients, aes(x = total_ddd, y = as.numeric(death_28d), color = infection_site)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess") +
  labs(title = "累计DDD与死亡风险的关系 (按感染部位分组)", 
       x = "累计DDD", 
       y = "28天死亡(1=是, 0=否)",
       color = "感染部位") +
  theme_minimal()

print(p11)

# ================================
# 7. 抗生素联合方式与耐药模式分析
# ================================

# 7.1 全人群联合方式分析
# 不同联合用药模式对耐药风险和临床结局的影响
combination_vs_outcome_all <- patients %>%
  group_by(combination_mode) %>%
  summarize(
    n_patients = n(),
    resistance_rate = mean(resistance_detected),
    death_rate = mean(death_28d),
    treatment_failure_rate = mean(treatment_failure)
  )

print("全人群联合用药模式与临床结局:")
print(combination_vs_outcome_all)

# 7.2 联合方式详细分类
# 创建更详细的联合方式分类
antibiotics$class_spectrum <- ifelse(antibiotics$spectrum == "Broad", "广谱", 
                                     ifelse(antibiotics$spectrum == "Medium", "中谱", "窄谱"))

# 按照联合抗生素的谱进行分类
determine_spectrum_combination <- function(pid) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  
  if(nrow(patient_abs) == 1) {
    return(paste0("单一", patient_abs$class_spectrum))
  }
  
  # 统计不同谱的抗生素数量
  spectrum_counts <- table(patient_abs$class_spectrum)
  
  # 判断联合类型
  if(length(spectrum_counts) == 1) {
    # 单一谱联合
    return(paste0("同谱联合-", names(spectrum_counts)))
  } else if("广谱" %in% names(spectrum_counts) && "窄谱" %in% names(spectrum_counts)) {
    # 广谱+窄谱
    return("广窄联合")
  } else if("广谱" %in% names(spectrum_counts) && "中谱" %in% names(spectrum_counts)) {
    # 广谱+中谱
    return("广中联合")
  } else if("中谱" %in% names(spectrum_counts) && "窄谱" %in% names(spectrum_counts)) {
    # 中谱+窄谱
    return("中窄联合")
  } else {
    # 三谱联合
    return("全谱联合")
  }
}

patients$spectrum_combination <- sapply(patients$patient_id, determine_spectrum_combination)

# 7.3 亚组联合方式分析
# 计算每个亚组中不同联合方式的临床结局
combination_stats_by_subgroup <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    
    # 计算联合方式与结局
    sg_combination_stats <- sg_data %>%
      group_by(combination_mode) %>%
      summarize(
        n_patients = n(),
        resistance_rate = mean(resistance_detected),
        death_rate = mean(death_28d),
        treatment_failure_rate = mean(treatment_failure)
      )
    
    combination_stats_by_subgroup[[paste0(sg, "_", group)]] <- sg_combination_stats
  }
}

# 显示示例亚组分析结果
print("重症亚组联合用药模式与临床结局:")
print(combination_stats_by_subgroup[["Severity_Severe (SOFA≥10)"]])

print("非重症亚组联合用药模式与临床结局:")
print(combination_stats_by_subgroup[["Severity_Non-severe (SOFA<10)"]])

# 7.4 联合方式与耐药谱系的关系分析
# 计算不同联合方式下各耐药类型的检出率
combination_resistance_types <- data.frame()

for(combo in unique(patients$spectrum_combination)) {
  # 使用该联合方式的患者ID
  patients_with_combo <- patients$patient_id[patients$spectrum_combination == combo]
  
  # 计算各类耐药菌检出率
  for(res_type in c("ESBL", "CRE", "MRSA", "VRE", "MDR-Acinetobacter")) {
    # 计算该耐药类型在这组患者中的检出率
    detection_count <- 0
    for(pid in patients_with_combo) {
      patient_micro <- microbiology[microbiology$patient_id == pid, ]
      if(nrow(patient_micro) > 0 && any(patient_micro$resistance == res_type)) {
        detection_count <- detection_count + 1
      }
    }
    
    detection_rate <- detection_count / max(length(patients_with_combo), 1)
    
    combination_resistance_types <- rbind(combination_resistance_types, data.frame(
      combination = combo,
      resistance_type = res_type,
      detection_rate = detection_rate,
      patient_count = length(patients_with_combo)
    ))
  }
}

# 创建热图数据
combination_resistance_wide <- combination_resistance_types %>%
  pivot_wider(names_from = resistance_type, values_from = detection_rate)

# 可视化联合方式与耐药谱系的关系
if(nrow(combination_resistance_wide) > 1) {
  heatmap_matrix <- as.matrix(combination_resistance_wide[, 3:7])
  rownames(heatmap_matrix) <- combination_resistance_wide$combination
  
  # 如果有足够的数据点，创建热图
  pheatmap(heatmap_matrix, 
           display_numbers = TRUE, 
           number_format = "%.2f",
           main = "抗生素联合方式与耐药谱系的关系",
           fontsize = 10,
           fontsize_number = 10)
}

# 7.5 亚组联合方式对临床结局的可视化比较
# 选择一个亚组(严重程度)进行联合方式比较可视化
severe_data <- patients[patients$severity_group == "Severe (SOFA≥10)", ]
non_severe_data <- patients[patients$severity_group == "Non-severe (SOFA<10)", ]

# 计算各种联合方式的死亡率
severe_combo_death <- severe_data %>%
  group_by(combination_mode) %>%
  summarize(
    n_patients = n(),
    death_rate = mean(death_28d)
  )

non_severe_combo_death <- non_severe_data %>%
  group_by(combination_mode) %>%
  summarize(
    n_patients = n(),
    death_rate = mean(death_28d)
  )

# 将两组结果合并为一个数据框进行可视化
combo_death_comparison <- rbind(
  cbind(severe_combo_death, group = "Severe (SOFA≥10)"),
  cbind(non_severe_combo_death, group = "Non-severe (SOFA<10)")
)

# 可视化比较
p12 <- ggplot(combo_death_comparison, aes(x = combination_mode, y = death_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste0(round(death_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "不同严重程度亚组中抗生素联合方式与死亡风险", 
       x = "抗生素联合方式", 
       y = "28天死亡率",
       fill = "严重程度") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p12)

# Logistic回归分析联合方式与耐药风险的亚组差异
combination_resistance_models <- list()

# 全人群模型
combination_resistance_models[["All"]] <- glm(
  resistance_detected ~ combination_mode + age + SOFA_score,
  family = binomial, data = patients
)

# 亚组模型
for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_model_name <- paste0(sg, "_", group)
    
    combination_resistance_models[[sg_model_name]] <- glm(
      resistance_detected ~ combination_mode + age + SOFA_score,
      family = binomial, data = sg_data
    )
  }
}

# 提取全人群和亚组模型的系数
combination_coef_table <- data.frame(
  Model = character(),
  Parallel_vs_Single = numeric(),
  Sequential_vs_Single = numeric(),
  Mixed_vs_Single = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(combination_resistance_models)) {
  model <- combination_resistance_models[[model_name]]
  
  # 提取系数(如果存在)
  parl_coef <- if("combination_modeParallel" %in% names(coef(model))) exp(coef(model)["combination_modeParallel"]) else NA
  seq_coef <- if("combination_modeSequential" %in% names(coef(model))) exp(coef(model)["combination_modeSequential"]) else NA
  mixed_coef <- if("combination_modeMixed" %in% names(coef(model))) exp(coef(model)["combination_modeMixed"]) else NA
  
  combination_coef_table <- rbind(combination_coef_table, data.frame(
    Model = model_name,
    Parallel_vs_Single = round(parl_coef, 2),
    Sequential_vs_Single = round(seq_coef, 2),
    Mixed_vs_Single = round(mixed_coef, 2),
    stringsAsFactors = FALSE
  ))
}

print("全人群和亚组联合方式的耐药风险比较(OR):")
print(combination_coef_table)

# ================================
# 8. 序贯模式分析
# ================================

# 8.1 创建序贯模式数据
# 提取有序贯治疗的患者
sequential_patients <- patients$patient_id[patients$combination_mode == "Sequential"]

# 创建序贯模式数据框
sequential_patterns <- data.frame()

for(pid in sequential_patients) {
  # 获取该患者的抗生素使用记录并按时间排序
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  patient_abs <- patient_abs[order(patient_abs$start_time), ]
  
  # 提取序贯模式
  if(nrow(patient_abs) >= 2) {
    # 创建序贯对
    for(i in 1:(nrow(patient_abs)-1)) {
      # 检查是否确实是序贯(第一种结束后才开始使用第二种)
      if(patient_abs$end_time[i] <= patient_abs$start_time[i+1]) {
        sequential_patterns <- rbind(sequential_patterns, data.frame(
          patient_id = pid,
          from_class = patient_abs$antibiotic_class[i],
          to_class = patient_abs$antibiotic_class[i+1],
          from_spectrum = patient_abs$spectrum[i],
          to_spectrum = patient_abs$spectrum[i+1],
          gap_time = patient_abs$start_time[i+1] - patient_abs$end_time[i],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# 8.2 序贯方向分类(升级/降级/水平)
if(nrow(sequential_patterns) > 0) {
  sequential_patterns$direction <- ifelse(
    sequential_patterns$from_spectrum == "Narrow" & sequential_patterns$to_spectrum == "Broad" |
      sequential_patterns$from_spectrum == "Medium" & sequential_patterns$to_spectrum == "Broad" |
      sequential_patterns$from_spectrum == "Narrow" & sequential_patterns$to_spectrum == "Medium",
    "升级",
    ifelse(
      sequential_patterns$from_spectrum == "Broad" & sequential_patterns$to_spectrum == "Narrow" |
        sequential_patterns$from_spectrum == "Broad" & sequential_patterns$to_spectrum == "Medium" |
        sequential_patterns$from_spectrum == "Medium" & sequential_patterns$to_spectrum == "Narrow",
      "降级",
      "水平"
    )
  )
  
  # 统计各患者的序贯升级/降级/水平次数
  patient_seq_stats <- sequential_patterns %>%
    group_by(patient_id) %>%
    summarize(
      escalation_count = sum(direction == "升级"),
      de_escalation_count = sum(direction == "降级"),
      lateral_count = sum(direction == "水平"),
      total_seq_count = n()
    )
  
  # 合并到患者数据
  patients <- merge(patients, patient_seq_stats, by = "patient_id", all.x = TRUE)
  
  # 填充NA值为0(没有序贯的患者)
  patients$escalation_count[is.na(patients$escalation_count)] <- 0
  patients$de_escalation_count[is.na(patients$de_escalation_count)] <- 0
  patients$lateral_count[is.na(patients$lateral_count)] <- 0
  patients$total_seq_count[is.na(patients$total_seq_count)] <- 0
  
  # 分析序贯方向与临床结局
  # 计算有升级、降级或水平序贯的患者的结局
  patients$has_escalation <- patients$escalation_count > 0
  patients$has_de_escalation <- patients$de_escalation_count > 0
  patients$has_lateral <- patients$lateral_count > 0
  
  # 统计不同序贯方向的临床结局
  direction_stats_all <- data.frame(
    Direction = c("有升级", "无升级", "有降级", "无降级", "有水平", "无水平"),
    N_patients = c(
      sum(patients$has_escalation),
      sum(!patients$has_escalation),
      sum(patients$has_de_escalation),
      sum(!patients$has_de_escalation),
      sum(patients$has_lateral),
      sum(!patients$has_lateral)
    ),
    Death_rate = c(
      mean(patients$death_28d[patients$has_escalation]),
      mean(patients$death_28d[!patients$has_escalation]),
      mean(patients$death_28d[patients$has_de_escalation]),
      mean(patients$death_28d[!patients$has_de_escalation]),
      mean(patients$death_28d[patients$has_lateral]),
      mean(patients$death_28d[!patients$has_lateral])
    ),
    Resistance_rate = c(
      mean(patients$resistance_detected[patients$has_escalation]),
      mean(patients$resistance_detected[!patients$has_escalation]),
      mean(patients$resistance_detected[patients$has_de_escalation]),
      mean(patients$resistance_detected[!patients$has_de_escalation]),
      mean(patients$resistance_detected[patients$has_lateral]),
      mean(patients$resistance_detected[!patients$has_lateral])
    )
  )
  
  print("全人群序贯方向与临床结局:")
  print(direction_stats_all)
  
  # 8.3 亚组序贯模式分析
  # 创建亚组序贯方向临床结局表
  seq_direction_by_subgroup <- list()
  
  for (sg in names(subgroups)) {
    sg_var <- subgroups[[sg]][1]
    
    for (group in unique(patients[[sg_var]])) {
      sg_data <- patients[patients[[sg_var]] == group, ]
      sg_name <- paste0(sg, "_", group)
      
      # 统计亚组中不同序贯方向的临床结局
      sg_direction_stats <- data.frame(
        Direction = c("有升级", "无升级", "有降级", "无降级"),
        N_patients = c(
          sum(sg_data$has_escalation),
          sum(!sg_data$has_escalation),
          sum(sg_data$has_de_escalation),
          sum(!sg_data$has_de_escalation)
        ),
        Death_rate = c(
          mean(sg_data$death_28d[sg_data$has_escalation]),
          mean(sg_data$death_28d[!sg_data$has_escalation]),
          mean(sg_data$death_28d[sg_data$has_de_escalation]),
          mean(sg_data$death_28d[!sg_data$has_de_escalation])
        ),
        Resistance_rate = c(
          mean(sg_data$resistance_detected[sg_data$has_escalation]),
          mean(sg_data$resistance_detected[!sg_data$has_escalation]),
          mean(sg_data$resistance_detected[sg_data$has_de_escalation]),
          mean(sg_data$resistance_detected[!sg_data$has_de_escalation])
        )
      )
      
      seq_direction_by_subgroup[[sg_name]] <- sg_direction_stats
    }
  }
  
  # 显示重症亚组序贯方向结局
  print("重症亚组序贯方向与临床结局:")
  print(seq_direction_by_subgroup[["Severity_Severe (SOFA≥10)"]])
  
  # 8.4 可视化序贯方向与耐药风险
  p13 <- ggplot(direction_stats_all[1:4,], aes(x = Direction, y = Resistance_rate, fill = Direction)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Resistance_rate * 100, 1), "%")), 
              position = position_dodge(width = 0.9), vjust = -0.5) +
    scale_fill_manual(values = c("有升级" = "salmon", "无升级" = "lightblue", 
                                 "有降级" = "lightgreen", "无降级" = "pink")) +
    labs(title = "全人群: 抗生素序贯方向与耐药风险", 
         x = "序贯方向", 
         y = "耐药菌检出率") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p13)
  
  # 8.5 亚组序贯方向比较可视化
  # 准备比较数据
  comparison_data <- data.frame()
  
  # 添加重症和非重症亚组数据
  if(!is.null(seq_direction_by_subgroup[["Severity_Severe (SOFA≥10)"]]) && 
     !is.null(seq_direction_by_subgroup[["Severity_Non-severe (SOFA<10)"]])) {
    
    severe_data <- seq_direction_by_subgroup[["Severity_Severe (SOFA≥10)"]][1:2,]
    severe_data$Group <- "重症患者"
    
    non_severe_data <- seq_direction_by_subgroup[["Severity_Non-severe (SOFA<10)"]][1:2,]
    non_severe_data$Group <- "非重症患者"
    
    comparison_data <- rbind(severe_data, non_severe_data)
    
    # 可视化比较
    p14 <- ggplot(comparison_data, aes(x = Direction, y = Resistance_rate, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_text(aes(label = paste0(round(Resistance_rate * 100, 1), "%")), 
                position = position_dodge(width = 0.9), vjust = -0.5) +
      labs(title = "不同严重程度亚组中抗生素序贯升级与耐药风险", 
           x = "序贯方向", 
           y = "耐药菌检出率",
           fill = "患者分组") +
      theme_minimal()
    
    print(p14)
  }
  
  # 8.6 分析常见序贯模式及其结局
  common_patterns <- sequential_patterns %>%
    group_by(from_class, to_class) %>%
    summarize(
      count = n(),
      avg_gap = mean(gap_time),
      direction = first(direction)
    ) %>%
    filter(count >= 3)  # 至少出现3次的模式
  
  if(nrow(common_patterns) > 0) {
    print("常见序贯模式:")
    print(common_patterns)
    
    # 为每位患者确定最常用的序贯模式
    patients_with_seq <- unique(sequential_patterns$patient_id)
    patient_most_common_pattern <- data.frame()
    
    for(pid in patients_with_seq) {
      # 获取该患者的序贯模式
      patient_patterns <- sequential_patterns[sequential_patterns$patient_id == pid, ]
      
      # 统计最常见的序贯模式
      pattern_counts <- table(paste(patient_patterns$from_class, "→", patient_patterns$to_class))
      if(length(pattern_counts) > 0) {
        most_common <- names(which.max(pattern_counts))
        count <- max(pattern_counts)
        
        patient_most_common_pattern <- rbind(patient_most_common_pattern, data.frame(
          patient_id = pid,
          most_common_pattern = most_common,
          pattern_count = count,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # 合并到患者数据
    patients <- merge(patients, patient_most_common_pattern, by = "patient_id", all.x = TRUE)
    
    # 分析最常见序贯模式的临床结局
    pattern_outcomes <- patients %>%
      filter(!is.na(most_common_pattern)) %>%
      group_by(most_common_pattern) %>%
      summarize(
        n_patients = n(),
        death_rate = mean(death_28d),
        resistance_rate = mean(resistance_detected),
        treatment_failure_rate = mean(treatment_failure)
      ) %>%
      filter(n_patients >= 3)  # 至少3名患者
    
    print("常见序贯模式的临床结局:")
    print(pattern_outcomes)
  }
}

# ================================
# 9. 死亡预测模型构建与评估
# ================================

# 9.1 准备建模数据
# 选择预测变量
pred_vars <- c("age", "SOFA_score", "antibiotic_count", "special_antibiotic", "resistance_detected", 
               "comorbidity_count", "mechanical_ventilation", "vasopressors", "infection_site",
               "shannon_index")

# 如果时间特征和剂量特征存在，添加到预测变量
if("total_exposure_hours" %in% names(patients)) {
  pred_vars <- c(pred_vars, "total_exposure_hours", "change_frequency", "max_concurrent")
}

if("total_ddd" %in% names(patients)) {
  pred_vars <- c(pred_vars, "total_ddd", "ddd_per_antibiotic")
}

if("escalation_count" %in% names(patients)) {
  pred_vars <- c(pred_vars, "escalation_count", "de_escalation_count")
}

# 准备建模数据
model_data <- patients[, c("death_28d", pred_vars)]
model_data <- model_data[complete.cases(model_data), ]  # 移除含有缺失值的行

# 9.2 全人群LASSO回归进行特征选择
set.seed(1615)  # 设置随机数种子

x <- model.matrix(death_28d ~ ., data = model_data)[, -1]  # 移除截距
y <- model_data$death_28d

# 交叉验证确定最佳λ值
cv.lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
plot(cv.lasso)

# 提取1SE规则λ值下的非零系数
lambda_1se <- cv.lasso$lambda.1se
lasso_coef <- coef(cv.lasso, s = lambda_1se)
selected_features <- which(lasso_coef != 0)
selected_vars <- rownames(lasso_coef)[selected_features]
selected_vars <- selected_vars[selected_vars != "(Intercept)"]
print("LASSO回归选择的特征变量:")
print(selected_vars)

# 9.3 全人群模型构建与验证
set.seed(1615)

# 创建训练和测试集
train_index <- createDataPartition(model_data$death_28d, p = 0.7, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# 构建Logistic回归模型
formula_str <- paste("death_28d ~", paste(selected_vars, collapse = " + "))
logistic_model_all <- glm(as.formula(formula_str), family = binomial, data = train_data)
summary(logistic_model_all)

# 模型性能评估
# 训练集性能
train_pred_prob <- predict(logistic_model_all, type = "response")
train_pred <- ifelse(train_pred_prob > 0.5, 1, 0)
train_accuracy <- mean(train_pred == train_data$death_28d)
train_confusion <- table(Predicted = train_pred, Actual = train_data$death_28d)
train_sensitivity <- ifelse(sum(train_confusion[, 2]) > 0, 
                            train_confusion[2, 2] / sum(train_confusion[, 2]), NA)
train_specificity <- ifelse(sum(train_confusion[, 1]) > 0, 
                            train_confusion[1, 1] / sum(train_confusion[, 1]), NA)

# 测试集性能
test_pred_prob <- predict(logistic_model_all, newdata = test_data, type = "response")
test_pred <- ifelse(test_pred_prob > 0.5, 1, 0)
test_accuracy <- mean(test_pred == test_data$death_28d)
test_confusion <- table(Predicted = test_pred, Actual = test_data$death_28d)
test_sensitivity <- ifelse(sum(test_confusion[, 2]) > 0, 
                           test_confusion[2, 2] / sum(test_confusion[, 2]), NA)
test_specificity <- ifelse(sum(test_confusion[, 1]) > 0, 
                           test_confusion[1, 1] / sum(test_confusion[, 1]), NA)

# 打印模型性能
cat("全人群模型训练集性能:\n")
cat("准确率:", train_accuracy, "\n")
cat("敏感性:", train_sensitivity, "\n")
cat("特异性:", train_specificity, "\n")
cat("\n全人群模型测试集性能:\n")
cat("准确率:", test_accuracy, "\n")
cat("敏感性:", test_sensitivity, "\n")
cat("特异性:", test_specificity, "\n")

# 绘制ROC曲线
train_roc <- roc(train_data$death_28d, train_pred_prob)
test_roc <- roc(test_data$death_28d, test_pred_prob)

p15 <- ggroc(list(Training = train_roc, Testing = test_roc)) +
  scale_color_manual(values = c("Training" = "blue", "Testing" = "red")) +
  labs(title = "全人群死亡预测模型ROC曲线", 
       subtitle = paste("训练集AUC =", round(auc(train_roc), 3), 
                        ", 测试集AUC =", round(auc(test_roc), 3))) +
  theme_minimal()

print(p15)

# 9.4 亚组预测模型构建
# 为特定亚组构建预测模型
subgroup_models <- list()
subgroup_performance <- list()

# 选择特定亚组进行单独建模
key_subgroups <- c("Severity_Severe (SOFA≥10)", "Age_Elderly (≥65)")

for (sg_name in key_subgroups) {
  # 解析亚组名称
  sg_parts <- strsplit(sg_name, "_")[[1]]
  sg_var <- sg_parts[1]
  sg_val <- sg_parts[2]
  
  # 在patients中找到对应的变量
  if (sg_var == "Severity") {
    sg_var <- "severity_group"
  } else if (sg_var == "Age") {
    sg_var <- "age_group"
  }
  
  # 筛选亚组数据
  sg_patients <- patients[patients[[sg_var]] == sg_val, ]
  sg_model_data <- sg_patients[, c("death_28d", pred_vars)]
  sg_model_data <- sg_model_data[complete.cases(sg_model_data), ]
  
  if (nrow(sg_model_data) > 20) {  # 确保有足够的样本进行建模
    # 创建训练和测试集
    set.seed(1615)
    sg_train_index <- createDataPartition(sg_model_data$death_28d, p = 0.7, list = FALSE)
    sg_train_data <- sg_model_data[sg_train_index, ]
    sg_test_data <- sg_model_data[-sg_train_index, ]
    
    # 使用全人群模型的变量构建亚组模型
    sg_formula_str <- formula_str  # 使用与全人群相同的公式
    sg_model <- try(glm(as.formula(sg_formula_str), family = binomial, data = sg_train_data), silent = TRUE)
    
    if (!inherits(sg_model, "try-error")) {
      # 模型验证
      sg_test_pred_prob <- predict(sg_model, newdata = sg_test_data, type = "response")
      sg_test_roc <- try(roc(sg_test_data$death_28d, sg_test_pred_prob), silent = TRUE)
      
      if (!inherits(sg_test_roc, "try-error")) {
        sg_test_auc <- auc(sg_test_roc)
        
        # 存储模型和性能
        subgroup_models[[sg_name]] <- sg_model
        subgroup_performance[[sg_name]] <- list(
          AUC = sg_test_auc,
          ROC = sg_test_roc
        )
      }
    }
  }
}

# 9.5 全人群与亚组模型性能比较
performance_comparison <- data.frame(
  Model = c("全人群"),
  AUC = c(round(auc(test_roc), 3)),
  stringsAsFactors = FALSE
)

for (sg_name in names(subgroup_performance)) {
  performance_comparison <- rbind(performance_comparison, data.frame(
    Model = sg_name,
    AUC = round(subgroup_performance[[sg_name]]$AUC, 3),
    stringsAsFactors = FALSE
  ))
}

print("全人群与亚组模型性能比较:")
print(performance_comparison)

# 9.6 模型变量重要性分析
# 计算风险比
odds_ratios <- exp(coef(logistic_model_all))
ci <- exp(confint(logistic_model_all))
var_importance <- data.frame(
  Variable = names(odds_ratios),
  OR = odds_ratios,
  Lower_CI = ci[, 1],
  Upper_CI = ci[, 2]
)

# 可视化变量重要性
p16 <- ggplot(var_importance[-1, ], aes(x = reorder(Variable, OR), y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "全人群死亡预测模型变量重要性", 
       x = "变量", 
       y = "优势比(OR) 及 95%置信区间") +
  theme_minimal()

print(p16)

# 9.7 校准曲线
# 使用rms包创建校准曲线
lrm_model <- lrm(as.formula(formula_str), data = train_data, x=TRUE, y=TRUE)

# 校准曲线
cal <- calibrate(lrm_model, method="boot", B=200)
plot(cal, main="死亡预测模型校准曲线")

# ================================
# 10. 抗生素时序特征与联合模式的交互影响
# ================================

# 10.1 创建时序+联合模式的组合特征
# 联合抗生素时序特征和联合模式
if("first_use_time" %in% names(patients)) {
  patients$timing_combination <- paste(
    ifelse(patients$first_use_time <= 2, "早期启动", 
           ifelse(patients$first_use_time <= 6, "中期启动", "晚期启动")),
    patients$combination_mode
  )
  
  # 10.2 分析组合特征对耐药风险和死亡率的影响
  # 计算各组合特征下的耐药率和死亡率
  timing_combo_outcomes <- patients %>%
    group_by(timing_combination) %>%
    summarize(
      n_patients = n(),
      resistance_rate = mean(resistance_detected),
      death_rate = mean(death_28d),
      treatment_failure_rate = mean(treatment_failure)
    ) %>%
    filter(n_patients >= 3)  # 至少3名患者的组合
  
  print("时序+联合模式组合特征与临床结局:")
  print(timing_combo_outcomes)
  
  # 可视化组合特征与耐药风险的关系
  p17 <- ggplot(timing_combo_outcomes, aes(x = reorder(timing_combination, -resistance_rate), y = resistance_rate)) +
    geom_bar(stat = "identity", aes(fill = resistance_rate)) +
    scale_fill_gradient(low = "lightgreen", high = "darkred") +
    labs(title = "抗生素时序+联合模式与耐药风险", 
         x = "时序-联合模式组合", 
         y = "耐药菌检出率") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p17)
  
  # Logistic回归分析组合特征与耐药风险
  combo_resistance_model <- glm(
    resistance_detected ~ timing_combination + age + SOFA_score,
    family = binomial, data = patients
  )
  
  summary(combo_resistance_model)
  
  # 10.3 序贯模式与联合方式的耐药影响
  if(nrow(sequential_patterns) > 0) {
    # 计算每种药物对的耐药风险
    pair_resistance <- data.frame()
    
    unique_pairs <- unique(paste(sequential_patterns$from_class, "→", sequential_patterns$to_class))
    
    for(pair in unique_pairs) {
      # 提取药物名称
      split_pair <- strsplit(pair, " → ")[[1]]
      from_drug <- split_pair[1]
      to_drug <- split_pair[2]
      
      # 找出使用该药物对的患者
      patients_with_pair <- unique(sequential_patterns$patient_id[
        sequential_patterns$from_class == from_drug & sequential_patterns$to_class == to_drug
      ])
      
      # 计算这些患者的耐药率
      if(length(patients_with_pair) > 0) {
        resistance_rate <- mean(patients$resistance_detected[patients$patient_id %in% patients_with_pair])
        death_rate <- mean(patients$death_28d[patients$patient_id %in% patients_with_pair])
        
        pair_resistance <- rbind(pair_resistance, data.frame(
          pair = pair,
          n_patients = length(patients_with_pair),
          resistance_rate = resistance_rate,
          death_rate = death_rate
        ))
      }
    }
    
    # 筛选至少有3名患者的药物对
    common_pairs <- pair_resistance[pair_resistance$n_patients >= 3, ]
    
    if(nrow(common_pairs) > 0) {
      print("常见药物序贯对与临床结局:")
      print(common_pairs)
      
      # 可视化常见药物对的耐药风险
      p18 <- ggplot(common_pairs, aes(x = reorder(pair, -resistance_rate), y = resistance_rate)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
                  position = position_dodge(width = 0.9), vjust = -0.5) +
        labs(title = "常见抗生素序贯对与耐药风险", 
             x = "抗生素序贯对", 
             y = "耐药菌检出率") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(p18)
    }
  }
}

# ================================
# 11. 综合分析与可视化
# ================================

# 11.1 创建联合模型预测耐药风险
# 整合多种特征预测耐药风险
pred_vars_resistance <- c("antibiotic_count", "shannon_index", "special_antibiotic", 
                          "combination_mode", "age", "SOFA_score", "infection_site")

# 如果时间特征和剂量特征存在，添加到预测变量
if("total_exposure_hours" %in% names(patients)) {
  pred_vars_resistance <- c(pred_vars_resistance, "total_exposure_hours", "change_frequency", "max_concurrent")
}

if("total_ddd" %in% names(patients)) {
  pred_vars_resistance <- c(pred_vars_resistance, "total_ddd", "ddd_per_antibiotic")
}

# 创建公式
formula_resistance <- as.formula(paste("resistance_detected ~", paste(pred_vars_resistance, collapse = " + ")))

# 构建模型
full_resistance_model <- glm(formula_resistance, family = binomial, data = patients)
summary(full_resistance_model)

# 变量重要性分析
full_model_or <- exp(coef(full_resistance_model))
full_model_ci <- exp(confint(full_resistance_model))

full_model_importance <- data.frame(
  Variable = names(full_model_or),
  OR = full_model_or,
  Lower_CI = full_model_ci[, 1],
  Upper_CI = full_model_ci[, 2],
  P_value = summary(full_resistance_model)$coefficients[, 4]
)

# 筛选显著变量(P<0.05)
significant_vars <- full_model_importance[full_model_importance$P_value < 0.05, ]
print("耐药风险预测的显著变量:")
print(significant_vars)

# 11.2 创建综合的风险因素可视化
p19 <- ggplot(full_model_importance[-1, ], aes(x = reorder(Variable, OR), y = OR)) +
  geom_point(aes(color = P_value < 0.05), size = 3) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, color = P_value < 0.05), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("gray", "blue"), name = "统计学显著") +
  coord_flip() +
  labs(title = "耐药风险综合模型 - 变量重要性", 
       x = "变量", 
       y = "优势比(OR) 及 95%置信区间") +
  theme_minimal()

print(p19)

# 11.3 耐药风险随时间变化的趋势
# 假设微生物学数据中有检测时间
# 计算每天的耐药菌检出率
days <- 1:20  # 假设最长住院日为20天
resistance_by_day <- sapply(days, function(day) {
  # 计算到第day天为止检出耐药菌的患者比例
  detections_by_day <- sapply(1:n_patients, function(pid) {
    patient_micro <- microbiology[microbiology$patient_id == pid, ]
    if (nrow(patient_micro) > 0) {
      # 检查是否在第day天或之前检出耐药菌
      any(patient_micro$culture_time <= day*24 & patient_micro$resistance != "None")
    } else {
      FALSE
    }
  })
  mean(detections_by_day)
})

# 创建时间趋势数据框
time_trend <- data.frame(
  Day = days,
  Resistance_rate = resistance_by_day
)

# 绘制耐药率随时间变化的趋势
p20 <- ggplot(time_trend, aes(x = Day, y = Resistance_rate)) +
  geom_line(color = "darkred", size = 1) +
  geom_point() +
  labs(title = "耐药菌检出率随住院日变化趋势", 
       x = "住院日", 
       y = "累计耐药菌检出率") +
  theme_minimal()

print(p20)

# 11.4 创建综合风险评分系统
# 基于模型系数构建评分系统
coef_values <- coef(full_resistance_model)
coef_values <- coef_values[!is.na(coef_values)]  # 移除NA值

# 计算每个变量的评分权重(标准化系数)
max_coef <- max(abs(coef_values[-1]))  # 排除截距
score_weights <- round(10 * coef_values / max_coef)  # 标准化为-10到10的范围

# 打印评分系统
score_system <- data.frame(
  Variable = names(score_weights),
  Weight = score_weights,
  Original_coef = coef_values
)
print("耐药风险评分系统:")
print(score_system)

# 计算每个患者的风险评分
calculate_risk_score <- function(patient, weights) {
  score <- weights["(Intercept)"]  # 基础分数(截距)
  
  # 添加各变量的得分
  for (var in names(weights)[-1]) {  # 排除截距
    if (var %in% names(patient)) {
      if (is.numeric(patient[[var]])) {
        # 数值变量
        score <- score + patient[[var]] * weights[var]
      } else if (is.factor(patient[[var]]) || is.character(patient[[var]])) {
        # 因子或字符变量
        var_level <- paste0(var, patient[[var]])
        if (var_level %in% names(weights)) {
          score <- score + weights[var_level]
        }
      } else if (is.logical(patient[[var]])) {
        # 逻辑变量
        if (patient[[var]]) {
          score <- score + weights[var]
        }
      }
    }
  }
  
  return(score)
}

# 给每个患者计算风险评分
patients$risk_score <- sapply(1:nrow(patients), function(i) {
  calculate_risk_score(patients[i, ], score_weights)
})

# 评估风险评分的预测能力
risk_score_roc <- roc(patients$resistance_detected, patients$risk_score)
risk_score_auc <- auc(risk_score_roc)

p21 <- ggroc(risk_score_roc) +
  labs(title = paste("耐药风险评分ROC曲线, AUC =", round(risk_score_auc, 3)),
       x = "1-特异性", y = "敏感性") +
  theme_minimal()

print(p21)

# 计算不同评分区间的实际耐药风险
patients$score_group <- cut(patients$risk_score, 
                            breaks = quantile(patients$risk_score, probs = seq(0, 1, by = 0.25)),
                            labels = c("低风险", "中低风险", "中高风险", "高风险"),
                            include.lowest = TRUE)

risk_by_score <- patients %>%
  group_by(score_group) %>%
  summarize(
    n_patients = n(),
    resistance_rate = mean(resistance_detected)
  )

print("风险评分组与实际耐药率:")
print(risk_by_score)

p22 <- ggplot(risk_by_score, aes(x = score_group, y = resistance_rate, fill = score_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "风险评分与实际耐药率", 
       x = "风险评分组", 
       y = "耐药菌检出率") +
  theme_minimal() +
  theme(legend.position = "none")

print(p22)

# 11.5 不同亚组风险评分的预测能力比较
risk_score_by_subgroup <- list()

for (sg in names(subgroups)) {
  sg_var <- subgroups[[sg]][1]
  
  for (group in unique(patients[[sg_var]])) {
    sg_data <- patients[patients[[sg_var]] == group, ]
    sg_name <- paste0(sg, "_", group)
    
    # 计算亚组中风险评分的预测能力
    sg_risk_roc <- roc(sg_data$resistance_detected, sg_data$risk_score)
    sg_risk_auc <- auc(sg_risk_roc)
    
    risk_score_by_subgroup[[sg_name]] <- list(
      AUC = sg_risk_auc,
      ROC = sg_risk_roc
    )
  }
}

# 创建亚组风险评分性能比较表
risk_score_performance <- data.frame(
  Group = c("全人群"),
  AUC = c(round(risk_score_auc, 3)),
  stringsAsFactors = FALSE
)

for (sg_name in names(risk_score_by_subgroup)) {
  risk_score_performance <- rbind(risk_score_performance, data.frame(
    Group = sg_name,
    AUC = round(risk_score_by_subgroup[[sg_name]]$AUC, 3),
    stringsAsFactors = FALSE
  ))
}

print("不同亚组风险评分预测能力比较:")
print(risk_score_performance)

# ================================
# 12. 研究结论与临床建议
# ================================

# 12.1 创建关键研究发现总结
# 关键分析结果汇总
key_findings <- data.frame(
  Analysis = character(),
  Result = character(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# 抗生素种类数阈值
if(exists("threshold_all") && !is.na(threshold_all)) {
  key_findings <- rbind(key_findings, data.frame(
    Analysis = "抗生素种类数阈值点(全人群)",
    Result = paste0(round(threshold_all, 2), " (95% CI: ", 
                    round(boot_ci$basic[4], 2), "-", round(boot_ci$basic[5], 2), ")"),
    Significance = ifelse(threshold_all > 0, "超过此阈值耐药风险显著增加", "未发现显著阈值效应"),
    stringsAsFactors = FALSE
  ))
}

# 添加各亚组阈值点
for(i in 1:nrow(threshold_summary)) {
  key_findings <- rbind(key_findings, data.frame(
    Analysis = paste0("抗生素种类数阈值点(", threshold_summary$Group[i], ")"),
    Result = as.character(threshold_summary$Threshold[i]),
    Significance = paste0("阈值前后风险比变化: ", threshold_summary$OR_ratio[i], "倍"),
    stringsAsFactors = FALSE
  ))
}

# 时序特征与耐药风险
if(exists("time_results_all")) {
  # 首次使用时间
  first_use_or <- time_results_all$OR[time_results_all$Variable == "first_use_time"]
  first_use_p <- time_results_all$P_value[time_results_all$Variable == "first_use_time"]
  if(!is.na(first_use_or) && !is.na(first_use_p)) {
    key_findings <- rbind(key_findings, data.frame(
      Analysis = "首次抗生素使用时间与耐药风险",
      Result = paste0("OR = ", round(first_use_or, 2), " (每延迟1小时)"),
      Significance = ifelse(first_use_p < 0.05, "统计学显著", "无统计学显著性"),
      stringsAsFactors = FALSE
    ))
  }
  
  # 总暴露时间
  exposure_or <- time_results_all$OR[time_results_all$Variable == "total_exposure_hours"]
  exposure_p <- time_results_all$P_value[time_results_all$Variable == "total_exposure_hours"]
  if(!is.na(exposure_or) && !is.na(exposure_p)) {
    key_findings <- rbind(key_findings, data.frame(
      Analysis = "抗生素总暴露时间与耐药风险",
      Result = paste0("OR = ", round(exposure_or, 2), " (每增加1小时)"),
      Significance = ifelse(exposure_p < 0.05, "统计学显著", "无统计学显著性"),
      stringsAsFactors = FALSE
    ))
  }
}

# 剂量特征与耐药风险
if(exists("dose_results_all")) {
  # 累计DDD
  total_ddd_or <- dose_results_all$OR[dose_results_all$Variable == "total_ddd"]
  total_ddd_p <- dose_results_all$P_value[dose_results_all$Variable == "total_ddd"]
  if(!is.na(total_ddd_or) && !is.na(total_ddd_p)) {
    key_findings <- rbind(key_findings, data.frame(
      Analysis = "累计DDD与耐药风险",
      Result = paste0("OR = ", round(total_ddd_or, 2), " (每增加1个DDD)"),
      Significance = ifelse(total_ddd_p < 0.05, "统计学显著", "无统计学显著性"),
      stringsAsFactors = FALSE
    ))
  }
}

# 联合用药模式与耐药风险
if(exists("combination_vs_outcome_all")) {
  # 序贯vs单一
  seq_vs_single_rate_diff <- combination_vs_outcome_all$resistance_rate[combination_vs_outcome_all$combination_mode == "Sequential"] -
    combination_vs_outcome_all$resistance_rate[combination_vs_outcome_all$combination_mode == "Single"]
  
  # 并行vs序贯
  parallel_vs_seq_rate_diff <- combination_vs_outcome_all$resistance_rate[combination_vs_outcome_all$combination_mode == "Parallel"] -
    combination_vs_outcome_all$resistance_rate[combination_vs_outcome_all$combination_mode == "Sequential"]
  
  key_findings <- rbind(key_findings, data.frame(
    Analysis = "联合用药模式与耐药风险",
    Result = paste0("并行vs序贯: 耐药率差异 ", round(parallel_vs_seq_rate_diff * 100, 1), "%"),
    Significance = ifelse(abs(parallel_vs_seq_rate_diff) > 0.1, "临床显著差异", "差异不显著"),
    stringsAsFactors = FALSE
  ))
}

# 序贯方向与耐药风险
if(exists("direction_stats_all")) {
  # 升级vs降级
  esc_vs_deesc_rate_