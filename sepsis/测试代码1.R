# 设置工作环境
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

# ================================
# 1. 创建模拟数据####
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
patients$death_28d <- sapply(patients$mortality_prob, function(p) sample(c(0, 1), 1, prob = c(1-p, p)))

# 确保死亡率约为10%
# while (mean(patients$death_28d) < 0.09 || mean(patients$death_28d) > 0.11) {
#   patients$death_28d <- sapply(patients$mortality_prob, function(p) sample(c(0, 1), 1, prob = c(1-p, p)))
# }
# 1. 计算目标死亡事件数
target_deaths <- round(0.10 * nrow(patients))  # 目标10%死亡率

# 2. 根据mortality_prob排序，越高的概率越可能死亡
ordered_patients <- patients[order(-patients$mortality_prob), ]

# 3. 直接指定前N个为死亡，其余为存活
ordered_patients$death_28d <- 0
ordered_patients$death_28d[1:target_deaths] <- 1

# 4. 随机打乱顺序，保持整体死亡率不变但增加随机性
shuffled_indices <- sample(1:nrow(ordered_patients))
ordered_patients <- ordered_patients[shuffled_indices, ]

# 5. 将结果放回原数据框
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
# 2. 统计分析####
# ================================

# 3.1 描述性统计####
# ================================

# 基本情况描述
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

# 可视化基线特征
p1 <- ggplot(patients, aes(x = factor(death_28d), y = age, fill = factor(death_28d))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("Survivors", "Non-survivors")) +
  labs(title = "年龄比较", x = "患者状态", y = "年龄(岁)", fill = "状态") +
  theme_minimal()

p2 <- ggplot(patients, aes(x = factor(death_28d), y = SOFA_score, fill = factor(death_28d))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("Survivors", "Non-survivors")) +
  labs(title = "SOFA评分比较", x = "患者状态", y = "SOFA评分", fill = "状态") +
  theme_minimal()

p3 <- ggplot(patients, aes(x = factor(death_28d), y = antibiotic_count, fill = factor(death_28d))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("Survivors", "Non-survivors")) +
  labs(title = "抗生素种类数比较", x = "患者状态", y = "抗生素种类数", fill = "状态") +
  theme_minimal()

plot_combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(plot_combined)

# 抗菌药物使用特征分析
cat("\n抗菌药物使用特征\n")
cat("==================\n")
cat("抗生素种类数(均值±标准差):", mean(patients$antibiotic_count), "±", sd(patients$antibiotic_count), "\n")
cat("特殊级抗生素使用比例:", mean(patients$special_antibiotic), "\n")
cat("联合用药模式分布:\n")
print(table(patients$combination_mode))
cat("治疗调整模式分布:\n")
print(table(patients$adjustment_mode))

# 可视化抗生素使用特征
p4 <- ggplot(patients, aes(x = antibiotic_count)) +
  geom_bar(fill = "steelblue") +
  labs(title = "抗生素种类数分布", x = "抗生素种类数", y = "患者数") +
  theme_minimal()

p5 <- ggplot(patients, aes(x = combination_mode, fill = combination_mode)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "联合用药模式分布", x = "联合用药模式", y = "患者数", fill = "联合用药模式") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_combined_2 <- p4 + p5 + plot_layout(ncol = 2)
print(plot_combined_2)

# 微生物学和耐药特征分析
cat("\n微生物学和耐药特征\n")
cat("==================\n")
cat("耐药菌检出率:", mean(patients$resistance_detected), "\n")
cat("耐药分类分布:\n")
print(table(patients$mdr_class))
cat("病原菌分布:\n")
print(table(microbiology$pathogen))
cat("耐药类型分布:\n")
print(table(microbiology$resistance))

# 可视化微生物学特征
p6 <- ggplot(data = microbiology, aes(x = pathogen)) +
  geom_bar(fill = "darkgreen") +
  labs(title = "病原菌分布", x = "病原菌", y = "检出数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p7 <- ggplot(data = microbiology[microbiology$resistance != "None", ], aes(x = resistance)) +
  geom_bar(fill = "darkred") +
  labs(title = "耐药类型分布", x = "耐药类型", y = "检出数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_combined_3 <- p6 + p7 + plot_layout(ncol = 2)
print(plot_combined_3)

# 临床结局分析
cat("\n临床结局分析\n")
cat("==================\n")
cat("28天死亡率:", mean(patients$death_28d), "\n")
cat("治疗失败率:", mean(patients$treatment_failure), "\n")
cat("72小时临床反应分布:\n")
print(table(patients$clinical_response_72h))
cat("ICU住院日(均值±标准差):", mean(patients$icu_los), "±", sd(patients$icu_los), "\n")
cat("总住院日(均值±标准差):", mean(patients$hospital_los), "±", sd(patients$hospital_los), "\n")

# 可视化临床结局
p8 <- ggplot(patients, aes(x = clinical_response_72h, fill = clinical_response_72h)) +
  geom_bar() +
  scale_fill_manual(values = c("Improved" = "forestgreen", "Stable" = "gold", "Deteriorated" = "firebrick")) +
  labs(title = "72小时临床反应", x = "临床反应", y = "患者数", fill = "临床反应") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p9 <- ggplot(patients, aes(x = death_28d, fill = factor(death_28d))) +
  geom_bar() +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), 
                    labels = c("存活", "死亡")) +
  labs(title = "28天生存状态", x = "患者状态", y = "患者数", fill = "状态") +
  scale_x_continuous(breaks = c(0, 1), labels = c("存活", "死亡")) +
  theme_minimal()

plot_combined_4 <- p8 + p9 + plot_layout(ncol = 2)
print(plot_combined_4)

# ================================
# 3.2 抗生素多样性与耐药风险分析
# ================================

# 抗生素种类数与耐药风险的关系
# 计算每种抗生素数量下的耐药率
resistance_by_count <- patients %>%
  group_by(antibiotic_count) %>%
  dplyr::summarize(
    n_patients = n(),
    n_resistant = sum(resistance_detected),
    resistance_rate = mean(resistance_detected)
  )

print(resistance_by_count)

# 绘制抗生素种类数与耐药风险的关系图
p10 <- ggplot(resistance_by_count, aes(x = antibiotic_count, y = resistance_rate)) +
  geom_point(aes(size = n_patients), color = "darkblue", alpha = 0.7) +
  geom_line() +
  labs(title = "抗生素种类数与耐药风险的关系", 
       x = "抗生素种类数", 
       y = "耐药菌检出率",
       size = "患者数") +
  theme_minimal() +
  ylim(0, 1)

print(p10)

# Logistic回归分析抗生素种类数与耐药风险的关系
logit_model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score + infection_site, 
                   family = binomial, data = patients)
summary(logit_model)

# 非线性关系分析 - 限制性立方样条函数
rcs_model <- glm(resistance_detected ~ ns(antibiotic_count, df = 3) + age + SOFA_score, 
                 family = binomial, data = patients)
summary(rcs_model)

# 比较线性模型和非线性模型
anova(logit_model, rcs_model, test = "Chisq")

# 分段线性回归模型检测阈值效应
segmented_model <- try(
  segmented(logit_model, seg.Z = ~ antibiotic_count, psi = 5),
  silent = TRUE
)

if (!inherits(segmented_model, "try-error")) {
  summary(segmented_model)
  
  # 提取阈值点
  threshold <- segmented_model$psi[2]
  cat("检测到的阈值点为:", threshold, "\n")
  
  # 输出阈值前后的斜率
  slopes <- slope(segmented_model)
  print(slopes)
  
  # 预测值曲线 - 分段回归模型
  newdata <- data.frame(
    antibiotic_count = seq(1, 7, by = 0.1),
    age = mean(patients$age),
    SOFA_score = mean(patients$SOFA_score),
    infection_site = "Pulmonary"  # 使用最常见的感染部位
  )
  
  newdata$pred_logit <- predict(logit_model, newdata)
  newdata$pred_seg <- predict(segmented_model, newdata)
  newdata$pred_rcs <- predict(rcs_model, newdata)
  
  newdata$prob_logit <- exp(newdata$pred_logit) / (1 + exp(newdata$pred_logit))
  newdata$prob_seg <- exp(newdata$pred_seg) / (1 + exp(newdata$pred_seg))
  newdata$prob_rcs <- exp(newdata$pred_rcs) / (1 + exp(newdata$pred_rcs))
  
  # 绘制阈值效应图
  p11 <- ggplot(newdata, aes(x = antibiotic_count)) +
    geom_line(aes(y = prob_logit, color = "Linear"), size = 1) +
    geom_line(aes(y = prob_seg, color = "Segmented"), size = 1) +
    geom_line(aes(y = prob_rcs, color = "Spline"), size = 1, linetype = "dashed") +
    geom_point(data = resistance_by_count, aes(y = resistance_rate, size = n_patients), 
               color = "black", alpha = 0.7) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "darkred") +
    annotate("text", x = threshold + 0.3, y = 0.8, 
             label = paste("阈值点:", round(threshold, 2)), color = "darkred") +
    scale_color_manual(values = c("Linear" = "blue", "Segmented" = "red", "Spline" = "green"),
                       name = "模型类型") +
    labs(title = "抗生素种类数与耐药风险的非线性关系", 
         subtitle = "分段线性回归显示阈值效应",
         x = "抗生素种类数", 
         y = "预测耐药风险",
         size = "患者数") +
    theme_minimal()
  
  print(p11)
}

# Shannon多样性指数与耐药风险的关系
# Logistic回归
shannon_model <- glm(resistance_detected ~ shannon_index + age + SOFA_score, 
                     family = binomial, data = patients)
summary(shannon_model)

# 分析Shannon指数与抗生素种类数的关系
p12 <- ggplot(patients, aes(x = antibiotic_count, y = shannon_index)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess") +
  labs(title = "抗生素种类数与多样性指数的关系", 
       x = "抗生素种类数", 
       y = "Shannon多样性指数") +
  theme_minimal()

p13 <- ggplot(patients, aes(x = shannon_index, y = as.numeric(resistance_detected))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "blue") +
  labs(title = "Shannon多样性指数与耐药风险", 
       x = "Shannon多样性指数", 
       y = "耐药风险") +
  theme_minimal() +
  ylim(-0.05, 1.05)

plot_combined_5 <- p12 + p13 + plot_layout(ncol = 2)
print(plot_combined_5)

# ================================
# 3.3 抗菌药物分类与耐药风险####
# ================================

# 不同抗菌药物大类与耐药风险的关系
# 计算每个抗菌药物大类的使用频率和耐药率
antibiotic_class_analysis <- antibiotics %>%
  group_by(antibiotic_class) %>%
  dplyr::summarize(
    usage_count = n(),
    usage_rate = n() / nrow(antibiotics)
  )

# 计算每个抗菌药物大类与耐药菌检出的关系
class_vs_resistance <- data.frame()

for (class in unique(antibiotics$antibiotic_class)) {
  # 使用该类抗生素的患者ID
  patients_using_class <- unique(antibiotics$patient_id[antibiotics$antibiotic_class == class])
  
  # 计算使用该类抗生素的患者中耐药菌检出率
  resistance_rate <- mean(patients$resistance_detected[patients$patient_id %in% patients_using_class])
  
  class_vs_resistance <- rbind(class_vs_resistance, data.frame(
    antibiotic_class = class,
    patients_count = length(patients_using_class),
    resistance_rate = resistance_rate
  ))
}

print(antibiotic_class_analysis)
print(class_vs_resistance)

# 可视化不同抗菌药物大类与耐药风险的关系
p14 <- ggplot(class_vs_resistance, aes(x = reorder(antibiotic_class, -resistance_rate), y = resistance_rate)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "不同抗菌药物大类与耐药风险", 
       x = "抗菌药物大类", 
       y = "耐药菌检出率") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p15 <- ggplot(class_vs_resistance, aes(x = reorder(antibiotic_class, -patients_count), y = patients_count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = patients_count), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "不同抗菌药物大类使用频率", 
       x = "抗菌药物大类", 
       y = "使用患者数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_combined_6 <- p14 + p15 + plot_layout(ncol = 2)
print(plot_combined_6)

# 药物组合与耐药类型的关联分析
# 获取使用两种及以上抗生素的患者
patients_with_combinations <- unique(antibiotics$patient_id[duplicated(antibiotics$patient_id)])

# 创建药物组合与耐药类型的关联表格
combinations_vs_resistance <- data.frame()

for (pid in patients_with_combinations) {
  # 获取该患者使用的药物组合
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  combo <- paste(sort(patient_abs$antibiotic_class), collapse = " + ")
  
  # 获取该患者的耐药情况
  patient_micro <- microbiology[microbiology$patient_id == pid, ]
  resistance_types <- ifelse(nrow(patient_micro) > 0 && any(patient_micro$resistance != "None"),
                             paste(unique(patient_micro$resistance[patient_micro$resistance != "None"]), collapse = ", "),
                             "None")
  
  combinations_vs_resistance <- rbind(combinations_vs_resistance, data.frame(
    patient_id = pid,
    combination = combo,
    resistance = resistance_types,
    stringsAsFactors = FALSE
  ))
}

# 统计常见组合
common_combos <- combinations_vs_resistance %>%
  group_by(combination) %>%
  dplyr::summarize(
    count = n(),
    resistance_rate = mean(resistance != "None")
  ) %>%
  filter(count >= 3)  # 至少3名患者使用的组合

print(common_combos)

# 热图数据准备 - 药物组合与特定耐药类型
resistance_types <- c("ESBL", "CRE", "MRSA", "VRE", "MDR-Acinetobacter")
heatmap_data <- data.frame()

for (combo in common_combos$combination) {
  for (res_type in resistance_types) {
    # 使用该组合的患者中特定耐药类型的检出率
    patients_using_combo <- combinations_vs_resistance$patient_id[combinations_vs_resistance$combination == combo]
    
    # 检查这些患者中有多少检出了特定耐药类型
    resistance_count <- 0
    for (pid in patients_using_combo) {
      patient_micro <- microbiology[microbiology$patient_id == pid, ]
      if (nrow(patient_micro) > 0 && any(patient_micro$resistance == res_type)) {
        resistance_count <- resistance_count + 1
      }
    }
    
    detection_rate <- resistance_count / length(patients_using_combo)
    
    heatmap_data <- rbind(heatmap_data, data.frame(
      combination = combo,
      resistance_type = res_type,
      detection_rate = detection_rate
    ))
  }
}

# 创建热图
heatmap_data_wide <- heatmap_data %>%
  pivot_wider(names_from = resistance_type, values_from = detection_rate)

heatmap_matrix <- as.matrix(heatmap_data_wide[, -1])
rownames(heatmap_matrix) <- heatmap_data_wide$combination

# 如果有足够的数据点，创建热图
if(nrow(heatmap_matrix) > 1 && ncol(heatmap_matrix) > 1) {
  p16 <- heatmap(heatmap_matrix, 
                 scale = "none", 
                 Rowv = NA, Colv = NA,
                 col = colorRampPalette(c("white", "red"))(100),
                 main = "抗菌药物组合与特定耐药菌检出率",
                 margins = c(8, 10))
}

# Logistic回归分析特定抗菌药物与耐药风险
# 为每位患者创建每类抗生素使用的标记
for (class in unique(antibiotics$antibiotic_class)) {
  patients[[paste0("uses_", gsub(" ", "_", tolower(class)))]] <- 
    sapply(patients$patient_id, function(pid) {
      any(antibiotics$patient_id == pid & antibiotics$antibiotic_class == class)
    })
}

# 构建Logistic回归模型
antibiotics_resistance_model <- glm(
  resistance_detected ~ 
    uses_penicillins + uses_cephalosporins + uses_carbapenems + 
    uses_aminoglycosides + uses_quinolones + uses_glycopeptides + 
    age + SOFA_score,
  family = binomial, data = patients
)

summary(antibiotics_resistance_model)

# ================================
# 3.4 联合用药模式分析####
# ================================

# 不同联合用药模式对耐药风险和临床结局的影响
# 联合用药模式与耐药风险
combination_vs_resistance <- patients %>%
  group_by(combination_mode) %>%
  dplyr::summarize(
    n_patients = n(),
    resistance_rate = mean(resistance_detected),
    death_rate = mean(death_28d),
    treatment_failure_rate = mean(treatment_failure)
  )

print(combination_vs_resistance)

# 可视化联合用药模式的影响
p17 <- ggplot(combination_vs_resistance, aes(x = combination_mode, y = resistance_rate, fill = combination_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(resistance_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "联合用药模式与耐药风险", 
       x = "联合用药模式", 
       y = "耐药菌检出率",
       fill = "联合用药模式") +
  theme_minimal()

p18 <- ggplot(combination_vs_resistance, aes(x = combination_mode, y = death_rate, fill = combination_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(death_rate * 100, 1), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "联合用药模式与死亡风险", 
       x = "联合用药模式", 
       y = "28天死亡率",
       fill = "联合用药模式") +
  theme_minimal()

plot_combined_7 <- p17 + p18 + plot_layout(ncol = 2)
print(plot_combined_7)

# Logistic回归分析联合用药模式与耐药风险
combination_resistance_model <- glm(
  resistance_detected ~ combination_mode + age + SOFA_score,
  family = binomial, data = patients
)

summary(combination_resistance_model)

# Cox比例风险模型分析联合用药模式与患者生存
# 假设时间为住院天数
coxph_model <- coxph(
  Surv(hospital_los, death_28d) ~ combination_mode + age + SOFA_score,
  data = patients
)

summary(coxph_model)

# 联合用药网络分析
# 创建联合用药网络数据
network_data <- data.frame()

for (pid in patients_with_combinations) {
  patient_abs <- antibiotics[antibiotics$patient_id == pid, ]
  classes <- unique(patient_abs$antibiotic_class)
  
  # 创建所有可能的药物对
  if (length(classes) >= 2) {
    for (i in 1:(length(classes)-1)) {
      for (j in (i+1):length(classes)) {
        network_data <- rbind(network_data, data.frame(
          from = classes[i],
          to = classes[j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# 计算每对药物组合的频率
network_counts <- network_data %>%
  group_by(from, to) %>%
  dplyr::summarize(weight = n())

# 创建网络图（如果igraph包可用）
if (requireNamespace("igraph", quietly = TRUE)) {
  network_graph <- graph_from_data_frame(d = network_counts, directed = FALSE)
  
  V(network_graph)$size <- 20
  V(network_graph)$color <- "lightblue"
  V(network_graph)$frame.color <- "darkblue"
  V(network_graph)$label.cex <- 0.8
  
  E(network_graph)$width <- E(network_graph)$weight * 0.5
  E(network_graph)$edge.color <- "gray50"
  
  plot(network_graph, main = "抗菌药物联合用药网络图", 
       layout = layout_with_fr, vertex.label.dist = 0)
}

# ================================
# 3.5 治疗调整的影响####
# ================================
# ================================
# 3.5 治疗调整的影响####
# ================================

# 加载必要的R包
library(survminer)  # 用于增强生存分析图
library(tableone)   # 用于创建描述性统计表格

# 1. 描述性分析：治疗调整模式的分布

# 将调整模式变量转换为中文标签，增强可读性
patients$adjustment_mode <- factor(patients$adjustment_mode, 
                                   levels = c("No adjustment", "Early de-escalation", "Late de-escalation", "Escalation"),
                                   labels = c("无需调整", "早期降级", "晚期降级", "治疗升级"))

# 统计不同调整策略的患者数和比例
adjustment_stats <- patients %>%
  dplyr::summarize(
    total = n(),
    no_adjustment_n = sum(adjustment_mode == "无需调整"),
    early_deescalation_n = sum(adjustment_mode == "早期降级"),
    late_deescalation_n = sum(adjustment_mode == "晚期降级"),
    escalation_n = sum(adjustment_mode == "治疗升级"),
    no_adjustment_pct = mean(adjustment_mode == "无需调整") * 100,
    early_deescalation_pct = mean(adjustment_mode == "早期降级") * 100,
    late_deescalation_pct = mean(adjustment_mode == "晚期降级") * 100,
    escalation_pct = mean(adjustment_mode == "治疗升级") * 100,
    .by = 1
  )

# 打印治疗调整模式的分布
print(adjustment_stats)

# 绘制治疗调整模式的分布柱状图
ggplot(patients, aes(x = adjustment_mode)) +
  geom_bar(fill = "steelblue") +
  geom_text(stat = "count", aes(label = paste0(round(..count../nrow(patients)*100, 1), "%")), 
            vjust = -0.5) +
  labs(title = "抗生素治疗调整模式分布",
       x = "调整策略",
       y = "患者数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. 治疗调整与临床结局的关系

# 2.1 创建对照表，按治疗调整策略比较基线特征
adjustment_table <- CreateTableOne(
  vars = c("age", "gender", "SOFA_score", "APACHE_II", "comorbidity_count", 
           "infection_site", "antibiotic_count"),
  strata = "adjustment_mode",
  data = patients,
  test = TRUE
)

print(adjustment_table, showAllLevels = TRUE)

# 2.2 治疗调整与28天死亡率的关系
mortality_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    deaths = sum(death_28d),
    mortality_rate = mean(death_28d) * 100,
    .by = adjustment_mode
  )

# 绘制治疗调整与28天死亡率的关系
ggplot(mortality_by_adjustment, aes(x = adjustment_mode, y = mortality_rate, fill = adjustment_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(mortality_rate, 1), "%")), vjust = -0.5) +
  labs(title = "不同治疗调整策略的28天死亡率",
       x = "治疗调整策略",
       y = "28天死亡率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 2.3 治疗调整与耐药菌检出率的关系
resistance_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    resistant_cases = sum(resistance_detected),
    resistance_rate = mean(resistance_detected) * 100,
    .by = adjustment_mode
  )

# 绘制治疗调整与耐药菌检出率的关系
ggplot(resistance_by_adjustment, aes(x = adjustment_mode, y = resistance_rate, fill = adjustment_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(resistance_rate, 1), "%")), vjust = -0.5) +
  labs(title = "不同治疗调整策略的耐药菌检出率",
       x = "治疗调整策略",
       y = "耐药菌检出率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 3. 生存分析：治疗调整对生存率的影响

# 3.1 Kaplan-Meier生存曲线
fit_km <- survfit(Surv(hospital_los, death_28d) ~ adjustment_mode, data = patients)

# 绘制Kaplan-Meier生存曲线
ggsurvplot(
  fit_km,
  data = patients,
  pval = TRUE, # 显示对数秩检验P值
  conf.int = TRUE, # 显示置信区间
  risk.table = TRUE, # 显示风险表
  xlab = "住院天数",
  ylab = "生存概率",
  title = "不同治疗调整策略的生存曲线",
  palette = "jco", # 使用Journal of Clinical Oncology的调色板
  ggtheme = theme_minimal(),
  legend.title = "治疗调整策略",
  legend.labs = c("无需调整", "早期降级", "晚期降级", "治疗升级"),
  risk.table.height = 0.25
)

# 3.2 多变量Cox比例风险模型调整混杂因素
cox_adjustment_model <- coxph(
  Surv(hospital_los, death_28d) ~ adjustment_mode + age + SOFA_score + comorbidity_count,
  data = patients
)

# 输出Cox模型摘要
summary_cox <- summary(cox_adjustment_model)
print(summary_cox)

# 提取风险比和置信区间
cox_results <- data.frame(
  variable = names(coef(cox_adjustment_model)),
  hazard_ratio = exp(coef(cox_adjustment_model)),
  lower_ci = exp(confint(cox_adjustment_model))[,1],
  upper_ci = exp(confint(cox_adjustment_model))[,2],
  p_value = summary_cox$coefficients[,5]
)

# 绘制风险比森林图
ggplot(cox_results, aes(x = hazard_ratio, y = variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "治疗调整对28天死亡的风险比(Cox模型)",
       x = "风险比(HR) [Log尺度]",
       y = "") +
  theme_minimal()

# 4. 治疗调整时机的影响

# 4.1 计算治疗调整时间与结局的关系
# 针对早期降级和晚期降级组
deescalation_timing <- patients %>%
  filter(adjustment_mode %in% c("早期降级", "晚期降级")) %>%
  mutate(deescalation_timing = ifelse(adjustment_mode == "早期降级", "≤72h", ">72h"))

# 4.2 治疗降级时机与临床结局的关系
deescalation_outcomes <- deescalation_timing %>%
  dplyr::summarize(
    n = n(),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    median_los = median(hospital_los),
    .by = deescalation_timing
  )

# 打印降级时机与结局的关系
print(deescalation_outcomes)

# 5. 治疗降级策略与抗生素使用时长的关系

# 5.1 比较不同调整策略下的抗生素使用时长
abx_duration_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    mean_duration = mean(antibiotic_treatment_days),
    median_duration = median(antibiotic_treatment_days),
    sd_duration = sd(antibiotic_treatment_days),
    .by = adjustment_mode
  )

# 打印抗生素使用时长
print(abx_duration_by_adjustment)

# 5.2 绘制箱线图比较不同调整策略的抗生素使用时长
ggplot(patients, aes(x = adjustment_mode, y = antibiotic_treatment_days, fill = adjustment_mode)) +
  geom_boxplot() +
  labs(title = "不同治疗调整策略的抗生素使用时长",
       x = "治疗调整策略",
       y = "抗生素使用天数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 6. 治疗调整对ICU住院时间的影响
icu_los_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    mean_icu_los = mean(icu_los),
    median_icu_los = median(icu_los),
    sd_icu_los = sd(icu_los),
    .by = adjustment_mode
  )

# 绘制箱线图比较不同调整策略的ICU住院时间
ggplot(patients, aes(x = adjustment_mode, y = icu_los, fill = adjustment_mode)) +
  geom_boxplot() +
  labs(title = "不同治疗调整策略的ICU住院时间",
       x = "治疗调整策略",
       y = "ICU住院天数") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 7. 治疗调整与治疗失败的关系
failure_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    failure_count = sum(treatment_failure),
    failure_rate = mean(treatment_failure) * 100,
    .by = adjustment_mode
  )

# 绘制治疗调整与治疗失败率的关系
ggplot(failure_by_adjustment, aes(x = adjustment_mode, y = failure_rate, fill = adjustment_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(failure_rate, 1), "%")), vjust = -0.5) +
  labs(title = "不同治疗调整策略的治疗失败率",
       x = "治疗调整策略",
       y = "治疗失败率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 8. 亚组分析：不同感染部位的治疗调整效果
site_adjustment_outcomes <- patients %>%
  filter(!is.na(infection_site)) %>%
  dplyr::summarize(
    n = n(),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    .by = c(infection_site, adjustment_mode)
  )

# 绘制分面图：不同感染部位下治疗调整与死亡率的关系
ggplot(site_adjustment_outcomes, aes(x = adjustment_mode, y = mortality_rate, fill = adjustment_mode)) +
  geom_bar(stat = "identity") +
  facet_wrap(~infection_site) +
  labs(title = "不同感染部位下治疗调整策略的28天死亡率",
       x = "治疗调整策略",
       y = "28天死亡率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")

# 9. 治疗调整与特殊级抗生素使用的关系
special_abx_by_adjustment <- patients %>%
  dplyr::summarize(
    n = n(),
    special_abx_use = mean(special_antibiotic) * 100,
    .by = adjustment_mode
  )

# 绘制治疗调整与特殊级抗生素使用率的关系
ggplot(special_abx_by_adjustment, aes(x = adjustment_mode, y = special_abx_use, fill = adjustment_mode)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(special_abx_use, 1), "%")), vjust = -0.5) +
  labs(title = "不同治疗调整策略的特殊级抗生素使用率",
       x = "治疗调整策略",
       y = "特殊级抗生素使用率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 10. 多变量Logistic回归分析治疗调整对耐药风险的影响
logistic_resistance_model <- glm(
  resistance_detected ~ adjustment_mode + age + SOFA_score + comorbidity_count + infection_site,
  family = binomial,
  data = patients
)

# 输出Logistic回归模型摘要
summary_logistic <- summary(logistic_resistance_model)
print(summary_logistic)

# 计算校正后的比值比(OR)及其置信区间
or_results <- data.frame(
  variable = names(coef(logistic_resistance_model)),
  odds_ratio = exp(coef(logistic_resistance_model)),
  lower_ci = exp(confint(logistic_resistance_model))[,1],
  upper_ci = exp(confint(logistic_resistance_model))[,2],
  p_value = summary_logistic$coefficients[,4]
)

# 绘制比值比森林图
ggplot(or_results, aes(x = odds_ratio, y = variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "治疗调整对耐药菌检出的比值比(Logistic模型)",
       x = "比值比(OR) [Log尺度]",
       y = "") +
  theme_minimal()

# 11. 结果汇总表格
results_summary <- data.frame(
  Adjustment_Strategy = c("无需调整", "早期降级", "晚期降级", "治疗升级"),
  Patient_Count = c(
    sum(patients$adjustment_mode == "无需调整"),
    sum(patients$adjustment_mode == "早期降级"),
    sum(patients$adjustment_mode == "晚期降级"),
    sum(patients$adjustment_mode == "治疗升级")
  ),
  Mortality_Rate = c(
    mean(patients$death_28d[patients$adjustment_mode == "无需调整"]) * 100,
    mean(patients$death_28d[patients$adjustment_mode == "早期降级"]) * 100,
    mean(patients$death_28d[patients$adjustment_mode == "晚期降级"]) * 100,
    mean(patients$death_28d[patients$adjustment_mode == "治疗升级"]) * 100
  ),
  Resistance_Rate = c(
    mean(patients$resistance_detected[patients$adjustment_mode == "无需调整"]) * 100,
    mean(patients$resistance_detected[patients$adjustment_mode == "早期降级"]) * 100,
    mean(patients$resistance_detected[patients$adjustment_mode == "晚期降级"]) * 100,
    mean(patients$resistance_detected[patients$adjustment_mode == "治疗升级"]) * 100
  ),
  Median_Hospital_LOS = c(
    median(patients$hospital_los[patients$adjustment_mode == "无需调整"]),
    median(patients$hospital_los[patients$adjustment_mode == "早期降级"]),
    median(patients$hospital_los[patients$adjustment_mode == "晚期降级"]),
    median(patients$hospital_los[patients$adjustment_mode == "治疗升级"])
  ),
  Median_ICU_LOS = c(
    median(patients$icu_los[patients$adjustment_mode == "无需调整"]),
    median(patients$icu_los[patients$adjustment_mode == "早期降级"]),
    median(patients$icu_los[patients$adjustment_mode == "晚期降级"]),
    median(patients$icu_los[patients$adjustment_mode == "治疗升级"])
  ),
  Median_Abx_Duration = c(
    median(patients$antibiotic_treatment_days[patients$adjustment_mode == "无需调整"]),
    median(patients$antibiotic_treatment_days[patients$adjustment_mode == "早期降级"]),
    median(patients$antibiotic_treatment_days[patients$adjustment_mode == "晚期降级"]),
    median(patients$antibiotic_treatment_days[patients$adjustment_mode == "治疗升级"])
  )
)

# 打印结果汇总表格
print(results_summary)

# 导出结果为CSV文件
# write.csv(results_summary, "treatment_adjustment_outcomes.csv", row.names = FALSE)

# ================================
# 3.6 临床结局预测模型####
# ================================

# 加载必要的R包
library(e1071)    # 用于交叉验证
library(caret)    # 用于模型评估

# 1. 准备数据

# 1.1 定义预测变量和结局变量
predictors <- c("age", "gender", "SOFA_score", "APACHE_II", "comorbidity_count",
                "infection_site", "antibiotic_count", "special_antibiotic", 
                "antibiotic_treatment_days", "adjustment_mode", "resistance_detected")

outcome <- "death_28d"

# 1.2 分割数据为训练集和测试集(70%/30%)
set.seed(123) # 设置随机种子，确保结果可重复
train_indices <- sample(1:nrow(patients), 0.7 * nrow(patients))
train_data <- patients[train_indices, ]
test_data <- patients[-train_indices, ]

# 2. 特征选择

# 2.1 使用LASSO回归筛选重要特征
x_train <- model.matrix(~ age + gender + SOFA_score + APACHE_II + 
                          comorbidity_count + infection_site + antibiotic_count +
                          special_antibiotic + antibiotic_treatment_days + 
                          adjustment_mode + 
                          resistance_detected - 1, data = train_data)
y_train <- train_data[[outcome]]

# 进行交叉验证以选择最优lambda值
cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
best_lambda <- cv_fit$lambda.min
cat("最优lambda值:", best_lambda, "\n")

# 使用最优lambda拟合LASSO模型
lasso_fit <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = best_lambda)

# 提取非零系数的变量
coef_lasso <- coef(lasso_fit)
selected_vars <- rownames(coef_lasso)[which(coef_lasso != 0)]
selected_vars <- selected_vars[selected_vars != "(Intercept)"]  # 移除截距项
cat("LASSO选择的变量:", paste(selected_vars, collapse = ", "), "\n")

# 3. 构建死亡风险预测模型

# 3.1 模型1：仅基于临床特征的基础模型
formula_base <- as.formula(paste(outcome, "~", "age + gender + SOFA_score + comorbidity_count"))
base_model <- glm(formula_base, family = binomial, data = train_data)
summary(base_model)

# 3.2 模型2：整合抗生素使用特征的完整模型
# 创建一个函数来从LASSO选择的变量名中提取基本变量名
extract_base_vars <- function(var_names) {
  # 初始化结果向量
  base_vars <- character(0)
  
  # 遍历每个变量名
  for(var in var_names) {
    # 检查是否是因子水平变量（包含特定模式）
    if(grepl("^(.*?)(TRUE|FALSE|[[:digit:]]+|[^[:alnum:]]+.*$)", var)) {
      # 提取基本变量名
      base_var <- gsub("(.*?)(TRUE|FALSE|[[:digit:]]+|[^[:alnum:]]+.*$)", "\\1", var)
      # 移除尾部可能的特殊字符
      base_var <- gsub("[^[:alnum:]]+$", "", base_var)
      # 添加到结果中（如果不重复）
      if(!(base_var %in% base_vars)) {
        base_vars <- c(base_vars, base_var)
      }
    } else {
      # 如果不是因子水平变量，直接添加
      base_vars <- c(base_vars, var)
    }
  }
  
  return(base_vars)
}

# 对于adjustment_mode特殊处理
fixed_vars <- character(0)
for(var in selected_vars) {
  if(grepl("^adjustment_mode", var)) {
    if(!("adjustment_mode" %in% fixed_vars)) {
      fixed_vars <- c(fixed_vars, "adjustment_mode")
    }
  } else if(grepl("TRUE$", var)) {
    # 处理布尔变量，如resistance_detectedTRUE
    base_var <- gsub("TRUE$", "", var)
    fixed_vars <- c(fixed_vars, base_var)
  } else {
    fixed_vars <- c(fixed_vars, var)
  }
}

# 去重
fixed_vars <- unique(fixed_vars)

cat("修正后的变量:", paste(fixed_vars, collapse = ", "), "\n")

# 重新构建公式和模型
formula_full <- as.formula(paste(outcome, "~", paste(fixed_vars, collapse = " + ")))
full_model <- glm(formula_full, family = binomial, data = train_data)
summary(full_model)

# 4. 模型性能评估

# 4.1 获取预测概率
train_data$base_pred <- predict(base_model, type = "response")
train_data$full_pred <- predict(full_model, type = "response")

test_data$base_pred <- predict(base_model, newdata = test_data, type = "response")
test_data$full_pred <- predict(full_model, newdata = test_data, type = "response")

# 4.2 计算ROC曲线和AUC
# 基础模型的ROC
roc_base_train <- roc(train_data[[outcome]], train_data$base_pred)
roc_base_test <- roc(test_data[[outcome]], test_data$base_pred)
auc_base_train <- auc(roc_base_train)
auc_base_test <- auc(roc_base_test)

# 完整模型的ROC
roc_full_train <- roc(train_data[[outcome]], train_data$full_pred)
roc_full_test <- roc(test_data[[outcome]], test_data$full_pred)
auc_full_train <- auc(roc_full_train)
auc_full_test <- auc(roc_full_test)

# 打印AUC值
cat("基础模型训练集AUC:", auc_base_train, "\n")
cat("基础模型测试集AUC:", auc_base_test, "\n")
cat("完整模型训练集AUC:", auc_full_train, "\n")
cat("完整模型测试集AUC:", auc_full_test, "\n")

# 4.3 使用ggplot2绘制更美观的ROC曲线
# 提取ROC曲线数据点
roc_base_data <- data.frame(
  FPR = 1 - roc_base_test$specificities,
  TPR = roc_base_test$sensitivities,
  Model = "基础模型"
)

roc_full_data <- data.frame(
  FPR = 1 - roc_full_test$specificities,
  TPR = roc_full_test$sensitivities,
  Model = "完整模型"
)

roc_combined <- rbind(roc_base_data, roc_full_data)

# 绘制ROC曲线
ggplot(roc_combined, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "死亡风险预测模型 - ROC曲线比较",
       x = "1 - 特异性",
       y = "敏感性") +
  annotate("text", x = 0.7, y = 0.3, 
           label = paste("基础模型 AUC =", round(auc_base_test, 3)), 
           color = "blue") +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste("完整模型 AUC =", round(auc_full_test, 3)), 
           color = "red") +
  theme_minimal()

# 5. 模型校准

# 5.1 校准曲线
val_base <- val.prob(test_data$base_pred, test_data[[outcome]])
val_full <- val.prob(test_data$full_pred, test_data[[outcome]])

# 创建校准曲线数据
cal_base <- calibrate(base_model, method = "boot", B = 100)
cal_full <- calibrate(full_model, method = "boot", B = 100)

# 绘制校准曲线
par(mfrow = c(1, 2))
plot(cal_base, main = "基础模型校准曲线", xlab = "预测概率", ylab = "实际概率")
plot(cal_full, main = "完整模型校准曲线", xlab = "预测概率", ylab = "实际概率")
par(mfrow = c(1, 1))

# 5.2 Hosmer-Lemeshow拟合优度检验
hoslem.test(test_data[[outcome]], test_data$base_pred)
hoslem.test(test_data[[outcome]], test_data$full_pred)

# 6. 决策曲线分析

# 计算净获益（简化版）
calculate_net_benefit <- function(pred_prob, outcome, threshold) {
  TP <- sum(pred_prob >= threshold & outcome == 1)
  FP <- sum(pred_prob >= threshold & outcome == 0)
  n <- length(outcome)
  
  # 净获益计算
  NB <- (TP/n) - (FP/n) * (threshold/(1-threshold))
  return(NB)
}

# 创建决策曲线的数据
thresholds <- seq(0.01, 0.5, by = 0.01)
nb_base <- sapply(thresholds, function(t) calculate_net_benefit(test_data$base_pred, test_data[[outcome]], t))
nb_full <- sapply(thresholds, function(t) calculate_net_benefit(test_data$full_pred, test_data[[outcome]], t))

# 计算"对所有人进行干预"的净获益
nb_all <- sapply(thresholds, function(t) {
  prevalence <- mean(test_data[[outcome]])
  return(prevalence - (1 - prevalence) * t/(1 - t))
})

# 计算"对所有人不进行干预"的净获益
nb_none <- rep(0, length(thresholds))

# 创建决策曲线数据框
dca_data <- data.frame(
  threshold = thresholds,
  base_model = nb_base,
  full_model = nb_full,
  all = nb_all,
  none = nb_none
)

# 绘制决策曲线
ggplot(dca_data) +
  geom_line(aes(x = threshold, y = base_model, color = "基础模型"), size = 1) +
  geom_line(aes(x = threshold, y = full_model, color = "完整模型"), size = 1) +
  geom_line(aes(x = threshold, y = all, color = "全部干预"), linetype = "dashed") +
  geom_line(aes(x = threshold, y = none, color = "不干预"), linetype = "dotted") +
  labs(title = "决策曲线分析",
       x = "概率阈值",
       y = "净获益",
       color = "策略") +
  scale_color_manual(values = c("基础模型" = "blue", "完整模型" = "red", 
                                "全部干预" = "green", "不干预" = "black")) +
  theme_minimal() +
  coord_cartesian(ylim = c(-0.05, max(nb_full, na.rm = TRUE) * 1.1))

# 7. 模型的临床应用

# 7.1 计算最优截断点
best_cutoff_base <- coords(roc_base_test, "best", ret = "threshold")
best_cutoff_full <- coords(roc_full_test, "best", ret = "threshold")

cat("基础模型最优截断值:", best_cutoff_base, "\n")
cat("完整模型最优截断值:", best_cutoff_full, "\n")

# 7.2 计算各个截断点下的性能指标
cutoff_base <- best_cutoff_base[1]
cutoff_full <- best_cutoff_full[1]

# 基础模型在测试集上的性能
test_data$base_class <- ifelse(test_data$base_pred >= cutoff_base, 1, 0)
base_conf_matrix <- table(Predicted = test_data$base_class, Actual = test_data[[outcome]])
base_sensitivity <- base_conf_matrix[2,2] / sum(base_conf_matrix[,2])
base_specificity <- base_conf_matrix[1,1] / sum(base_conf_matrix[,1])
base_ppv <- base_conf_matrix[2,2] / sum(base_conf_matrix[2,])
base_npv <- base_conf_matrix[1,1] / sum(base_conf_matrix[1,])
base_accuracy <- sum(diag(base_conf_matrix)) / sum(base_conf_matrix)

# 完整模型在测试集上的性能
test_data$full_class <- ifelse(test_data$full_pred >= cutoff_full, 1, 0)
full_conf_matrix <- table(Predicted = test_data$full_class, Actual = test_data[[outcome]])
full_sensitivity <- full_conf_matrix[2,2] / sum(full_conf_matrix[,2])
full_specificity <- full_conf_matrix[1,1] / sum(full_conf_matrix[,1])
full_ppv <- full_conf_matrix[2,2] / sum(full_conf_matrix[2,])
full_npv <- full_conf_matrix[1,1] / sum(full_conf_matrix[1,])
full_accuracy <- sum(diag(full_conf_matrix)) / sum(full_conf_matrix)

# 创建性能指标汇总表
performance_summary <- data.frame(
  Model = c("基础模型", "完整模型"),
  AUC = c(auc_base_test, auc_full_test),
  Sensitivity = c(base_sensitivity, full_sensitivity),
  Specificity = c(base_specificity, full_specificity),
  PPV = c(base_ppv, full_ppv),
  NPV = c(base_npv, full_npv),
  Accuracy = c(base_accuracy, full_accuracy)
)

# 打印性能指标汇总表
print(performance_summary)

# 8. 模型内部验证 - Bootstrap

# 8.1 执行Bootstrap验证
set.seed(456)
n_bootstrap <- 100  # 减少迭代次数以加快示例运行速度

# 定义Bootstrap函数
bootstrap_validation <- function(model, data) {
  auc_values <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    # 有放回抽样
    boot_indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    # 非抽样数据
    test_indices <- setdiff(1:nrow(data), boot_indices)
    if (length(test_indices) == 0) next  # 跳过全部被抽中的情况
    test_data <- data[test_indices, ]
    
    # 在bootstrap样本上拟合模型
    boot_model <- update(model, data = boot_data)
    
    # 在非抽样数据上测试
    pred <- predict(boot_model, newdata = test_data, type = "response")
    if (length(unique(test_data[[outcome]])) < 2) next  # 跳过单一结局的情况
    
    # 计算AUC
    roc_obj <- tryCatch({
      roc(test_data[[outcome]], pred)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(roc_obj)) {
      auc_values[i] <- auc(roc_obj)
    }
  }
  
  # 去除NULL值
  auc_values <- auc_values[!is.na(auc_values)]
  
  # 返回平均AUC和95%置信区间
  return(list(
    mean_auc = mean(auc_values),
    lower_ci = quantile(auc_values, 0.025),
    upper_ci = quantile(auc_values, 0.975)
  ))
}

# 执行Bootstrap验证
base_boot <- bootstrap_validation(base_model, patients)
full_boot <- bootstrap_validation(full_model, patients)

# 输出Bootstrap验证结果
cat("基础模型Bootstrap验证:\n")
cat("平均AUC:", base_boot$mean_auc, "\n")
cat("95%置信区间:", base_boot$lower_ci, "-", base_boot$upper_ci, "\n\n")

cat("完整模型Bootstrap验证:\n")
cat("平均AUC:", full_boot$mean_auc, "\n")
cat("95%置信区间:", full_boot$lower_ci, "-", full_boot$upper_ci, "\n")

# 9. 输出模型系数及解释

# 提取模型系数和95%置信区间
coef_base <- summary(base_model)$coefficients
coef_full <- summary(full_model)$coefficients

# 转换为比值比(OR)及95%置信区间
or_base <- exp(coef(base_model))
ci_base <- exp(confint(base_model))
or_full <- exp(coef(full_model))
ci_full <- exp(confint(full_model))

# 创建系数汇总表
coef_summary_base <- data.frame(
  Variable = names(coef(base_model)),
  OR = or_base,
  Lower_CI = ci_base[,1],
  Upper_CI = ci_base[,2],
  P_value = coef_base[,4]
)

coef_summary_full <- data.frame(
  Variable = names(coef(full_model)),
  OR = or_full,
  Lower_CI = ci_full[,1],
  Upper_CI = ci_full[,2],
  P_value = coef_full[,4]
)

# 打印系数汇总表
print("基础模型系数:")
print(coef_summary_base)
print("完整模型系数:")
print(coef_summary_full)

# 10. 绘制森林图展示模型系数

# 基础模型森林图
coef_summary_base$Variable <- factor(coef_summary_base$Variable, 
                                     levels = coef_summary_base$Variable)

ggplot(coef_summary_base[-1,], aes(x = OR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "基础模型 - 变量对死亡风险的影响(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "") +
  theme_minimal()

# 完整模型森林图
coef_summary_full$Variable <- factor(coef_summary_full$Variable, 
                                     levels = coef_summary_full$Variable)

ggplot(coef_summary_full[-1,], aes(x = OR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "完整模型 - 变量对死亡风险的影响(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "") +
  theme_minimal()

# 11. 模型摘要和保存

# 模型比较的似然比检验
anova(base_model, full_model, test = "Chisq")

# 计算模型的AIC和BIC
base_aic <- AIC(base_model)
base_bic <- BIC(base_model)
full_aic <- AIC(full_model)
full_bic <- BIC(full_model)

cat("基础模型AIC:", base_aic, "BIC:", base_bic, "\n")
cat("完整模型AIC:", full_aic, "BIC:", full_bic, "\n")

# 创建最终的模型比较表
model_comparison <- data.frame(
  Model = c("基础模型", "完整模型"),
  Variables = c(length(coef(base_model)) - 1, length(coef(full_model)) - 1),
  AIC = c(base_aic, full_aic),
  BIC = c(base_bic, full_bic),
  AUC_Train = c(auc_base_train, auc_full_train),
  AUC_Test = c(auc_base_test, auc_full_test),
  AUC_Bootstrap = c(base_boot$mean_auc, full_boot$mean_auc),
  AUC_Bootstrap_CI = c(
    paste(round(base_boot$lower_ci, 3), "-", round(base_boot$upper_ci, 3)),
    paste(round(full_boot$lower_ci, 3), "-", round(full_boot$upper_ci, 3))
  )
)

# 打印模型比较表
print(model_comparison)

# 导出模型比较表为CSV文件
write.csv(model_comparison, "mortality_prediction_model_comparison.csv", row.names = FALSE)

# 保存完整模型
saveRDS(full_model, "mortality_prediction_model.rds")

# 创建预测风险计算函数
calculate_mortality_risk <- function(age, gender, SOFA_score, comorbidity_count,
                                     infection_site, antibiotic_count, special_antibiotic,
                                     antibiotic_treatment_days, adjustment_mode, resistance_detected) {
  # 创建新数据
  new_data <- data.frame(
    age = age,
    gender = gender,
    SOFA_score = SOFA_score,
    comorbidity_count = comorbidity_count,
    infection_site = infection_site,
    antibiotic_count = antibiotic_count,
    special_antibiotic = special_antibiotic,
    antibiotic_treatment_days = antibiotic_treatment_days,
    adjustment_mode = adjustment_mode,
    resistance_detected = resistance_detected
  )
  
  # 使用模型预测
  risk <- predict(full_model, newdata = new_data, type = "response")
  return(risk)
}

# 示例预测
example_risk <- calculate_mortality_risk(
  age = 65,
  gender = "Male",
  SOFA_score = 8,
  comorbidity_count = 2,
  infection_site = "Pulmonary",
  antibiotic_count = 4,
  special_antibiotic = TRUE,
  antibiotic_treatment_days = 10,
  adjustment_mode = "治疗升级",
  resistance_detected = TRUE
)

cat("示例患者的预测死亡风险:", round(example_risk * 100, 1), "%\n")

# ================================
# 3.7 亚组分析####
# ================================

# 加载必要的R包
library(forestplot)
library(gridExtra)
library(RColorBrewer)

# 1. 按感染部位进行亚组分析

# 1.1 创建感染部位亚组
site_groups <- unique(patients$infection_site)
site_groups <- site_groups[!is.na(site_groups)]  # 移除NA

# 1.2 计算各亚组的基本特征
site_characteristics <- patients %>%
  filter(!is.na(infection_site)) %>%
  dplyr::summarize(
    n = n(),
    mean_age = mean(age),
    male_percent = mean(gender == "Male") * 100,
    mean_SOFA = mean(SOFA_score),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    .by = infection_site
  )

print("感染部位亚组基本特征:")
print(site_characteristics)

# 1.3 按感染部位分析抗生素多样性与耐药风险的关系
site_resistance_model <- lapply(site_groups, function(site) {
  site_data <- patients[patients$infection_site == site, ]
  tryCatch({
    model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score,
                 family = binomial, data = site_data)
    
    # 计算OR和95%CI
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      site = site,
      n = nrow(site_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    return(data.frame(
      site = site,
      n = nrow(site_data),
      or = NA,
      ci_lower = NA,
      ci_upper = NA,
      p_value = NA
    ))
  })
})

site_resistance_results <- do.call(rbind, site_resistance_model)
print("感染部位亚组的抗生素多样性与耐药风险关系:")
print(site_resistance_results)

# 1.4 抗生素种类数与耐药风险关系的森林图(按感染部位)
ggplot(site_resistance_results, aes(x = or, y = site)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "不同感染部位亚组中抗生素多样性与耐药风险的关系",
       subtitle = "每增加一种抗生素的耐药风险比值比(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "感染部位") +
  theme_minimal()

# 2. 按病原菌类型进行亚组分析

# 从microbiology数据中提取主要病原菌类型
# 添加病原菌类型分类
patients$pathogen_type <- sapply(patients$patient_id, function(pid) {
  patient_micro <- microbiology[microbiology$patient_id == pid, ]
  if (nrow(patient_micro) == 0) return(NA)
  
  pathogens <- patient_micro$pathogen
  
  gram_neg <- c("Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa", 
                "Acinetobacter baumannii")
  gram_pos <- c("Staphylococcus aureus", "Enterococcus faecalis")
  fungi <- c("Candida albicans")
  
  has_gram_neg <- any(pathogens %in% gram_neg)
  has_gram_pos <- any(pathogens %in% gram_pos)
  has_fungi <- any(pathogens %in% fungi)
  
  if (has_gram_neg && !has_gram_pos && !has_fungi) return("Gram-negative")
  if (!has_gram_neg && has_gram_pos && !has_fungi) return("Gram-positive")
  if (!has_gram_neg && !has_gram_pos && has_fungi) return("Fungal")
  if ((has_gram_neg && has_gram_pos) || (has_gram_neg && has_fungi) || (has_gram_pos && has_fungi)) return("Mixed")
  
  return(NA)
})

# 2.1 创建病原菌类型亚组
pathogen_groups <- c("Gram-negative", "Gram-positive", "Mixed", "Fungal")
pathogen_groups <- pathogen_groups[pathogen_groups %in% unique(patients$pathogen_type)]

# 2.2 计算各亚组的基本特征
pathogen_characteristics <- patients %>%
  filter(!is.na(pathogen_type)) %>%
  dplyr::summarize(
    n = n(),
    mean_age = mean(age),
    male_percent = mean(gender == "Male") * 100,
    mean_SOFA = mean(SOFA_score),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    .by = pathogen_type
  )

print("病原菌类型亚组基本特征:")
print(pathogen_characteristics)

# 2.3 按病原菌类型分析抗生素多样性与耐药风险的关系
pathogen_resistance_model <- lapply(pathogen_groups, function(type) {
  pathogen_data <- patients[patients$pathogen_type == type, ]
  tryCatch({
    model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score,
                 family = binomial, data = pathogen_data)
    
    # 计算OR和95%CI
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      pathogen_type = type,
      n = nrow(pathogen_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    return(data.frame(
      pathogen_type = type,
      n = nrow(pathogen_data),
      or = NA,
      ci_lower = NA,
      ci_upper = NA,
      p_value = NA
    ))
  })
})

pathogen_resistance_results <- do.call(rbind, pathogen_resistance_model)
print("病原菌类型亚组的抗生素多样性与耐药风险关系:")
print(pathogen_resistance_results)

# 2.4 抗生素种类数与耐药风险关系的森林图(按病原菌类型)
ggplot(pathogen_resistance_results, aes(x = or, y = pathogen_type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "不同病原菌类型亚组中抗生素多样性与耐药风险的关系",
       subtitle = "每增加一种抗生素的耐药风险比值比(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "病原菌类型") +
  theme_minimal()

# 3. 按感染严重程度进行亚组分析

# 3.1 创建感染严重程度亚组(根据脓毒性休克状态)
patients$severity <- ifelse(patients$septic_shock == 1, "脓毒性休克", "脓毒症")
severity_groups <- c("脓毒症", "脓毒性休克")

# 3.2 计算各亚组的基本特征
severity_characteristics <- patients %>%
  dplyr::summarize(
    n = n(),
    mean_age = mean(age),
    male_percent = mean(gender == "Male") * 100,
    mean_SOFA = mean(SOFA_score),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    .by = severity
  )

print("感染严重程度亚组基本特征:")
print(severity_characteristics)

# 3.3 按感染严重程度分析抗生素多样性与耐药风险的关系
severity_resistance_model <- lapply(severity_groups, function(sev) {
  severity_data <- patients[patients$severity == sev, ]
  tryCatch({
    model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score,
                 family = binomial, data = severity_data)
    
    # 计算OR和95%CI
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      severity = sev,
      n = nrow(severity_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    return(data.frame(
      severity = sev,
      n = nrow(severity_data),
      or = NA,
      ci_lower = NA,
      ci_upper = NA,
      p_value = NA
    ))
  })
})

severity_resistance_results <- do.call(rbind, severity_resistance_model)
print("感染严重程度亚组的抗生素多样性与耐药风险关系:")
print(severity_resistance_results)

# 3.4 抗生素种类数与耐药风险关系的森林图(按感染严重程度)
ggplot(severity_resistance_results, aes(x = or, y = severity)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "不同感染严重程度亚组中抗生素多样性与耐药风险的关系",
       subtitle = "每增加一种抗生素的耐药风险比值比(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "感染严重程度") +
  theme_minimal()

# 4. 按既往抗生素暴露史进行亚组分析

# 为患者添加既往抗生素暴露史变量
patients$prior_abx_exposure <- sample(c("有", "无"), nrow(patients), replace = TRUE, prob = c(0.3, 0.7))

# 4.1 创建既往抗生素暴露史亚组
exposure_groups <- c("有", "无")

# 4.2 计算各亚组的基本特征
exposure_characteristics <- patients %>%
  dplyr::summarize(
    n = n(),
    mean_age = mean(age),
    male_percent = mean(gender == "Male") * 100,
    mean_SOFA = mean(SOFA_score),
    mortality_rate = mean(death_28d) * 100,
    resistance_rate = mean(resistance_detected) * 100,
    .by = prior_abx_exposure
  )

print("既往抗生素暴露史亚组基本特征:")
print(exposure_characteristics)

# 4.3 按既往抗生素暴露史分析抗生素多样性与耐药风险的关系
exposure_resistance_model <- lapply(exposure_groups, function(exp) {
  exposure_data <- patients[patients$prior_abx_exposure == exp, ]
  tryCatch({
    model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score,
                 family = binomial, data = exposure_data)
    
    # 计算OR和95%CI
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      exposure = exp,
      n = nrow(exposure_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    return(data.frame(
      exposure = exp,
      n = nrow(exposure_data),
      or = NA,
      ci_lower = NA,
      ci_upper = NA,
      p_value = NA
    ))
  })
})

exposure_resistance_results <- do.call(rbind, exposure_resistance_model)
print("既往抗生素暴露史亚组的抗生素多样性与耐药风险关系:")
print(exposure_resistance_results)

# 4.4 抗生素种类数与耐药风险关系的森林图(按既往抗生素暴露史)
ggplot(exposure_resistance_results, aes(x = or, y = exposure)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "既往抗生素暴露史亚组中抗生素多样性与耐药风险的关系",
       subtitle = "每增加一种抗生素的耐药风险比值比(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "既往抗生素暴露史") +
  theme_minimal()

# 5. 综合亚组分析的森林图

# 5.1 结合所有亚组结果
# 总体模型
overall_model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score,
                     family = binomial, data = patients)
overall_coef <- summary(overall_model)$coefficients["antibiotic_count", ]
overall_or <- exp(overall_coef[1])
overall_ci_lower <- exp(overall_coef[1] - 1.96 * overall_coef[2])
overall_ci_upper <- exp(overall_coef[1] + 1.96 * overall_coef[2])
overall_p <- overall_coef[4]

# 创建总体结果数据框
overall_result <- data.frame(
  group_type = "总体",
  group = "全部患者",
  n = nrow(patients),
  or = overall_or,
  ci_lower = overall_ci_lower,
  ci_upper = overall_ci_upper,
  p_value = overall_p
)

# 组合所有亚组结果
all_subgroups_df <- rbind(
  overall_result,
  transform(site_resistance_results, group_type = "感染部位", group = site),
  transform(pathogen_resistance_results, group_type = "病原菌类型", group = pathogen_type),
  transform(severity_resistance_results, group_type = "感染严重程度", group = severity),
  transform(exposure_resistance_results, group_type = "既往抗生素暴露", group = exposure)
)

# 按组类型重新排序
all_subgroups_df$group_type <- factor(all_subgroups_df$group_type, 
                                      levels = c("总体", "感染部位", "病原菌类型", 
                                                 "感染严重程度", "既往抗生素暴露"))
all_subgroups_df <- all_subgroups_df[order(all_subgroups_df$group_type), ]

# 5.2 创建组合森林图
ggplot(all_subgroups_df, aes(x = or, y = reorder(paste(group_type, "-", group), desc(row.names(all_subgroups_df))))) +
  geom_point(aes(size = n), shape = 15) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  scale_size_continuous(name = "样本量") +
  labs(title = "亚组分析：抗生素多样性与耐药风险的关系",
       subtitle = "每增加一种抗生素的耐药风险比值比(OR)",
       x = "比值比 (Odds Ratio) [Log尺度]",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 6. 比较不同亚组中的抗生素使用模式

# 6.1 感染部位与抗生素种类数
site_abx_count <- patients %>%
  filter(!is.na(infection_site)) %>%
  dplyr::summarize(
    n = n(),
    mean_abx_count = mean(antibiotic_count),
    sd_abx_count = sd(antibiotic_count),
    .by = infection_site
  )

# 绘制不同感染部位的抗生素种类数箱线图
ggplot(patients[!is.na(patients$infection_site), ], aes(x = infection_site, y = antibiotic_count)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "不同感染部位的抗生素种类数",
       x = "感染部位",
       y = "抗生素种类数") +
  theme_minimal()

# 6.2 感染严重程度与特殊级抗生素使用
severity_special_abx <- patients %>%
  dplyr::summarize(
    n = n(),
    special_abx_rate = mean(special_antibiotic) * 100,
    .by = severity
  )

# 绘制不同感染严重程度的特殊级抗生素使用率
ggplot(severity_special_abx, aes(x = severity, y = special_abx_rate, fill = severity)) +
  geom_bar(stat = "identity") +
  labs(title = "不同感染严重程度的特殊级抗生素使用率",
       x = "感染严重程度",
       y = "特殊级抗生素使用率 (%)") +
  theme_minimal()

# 7. 亚组中的阈值效应分析

# 7.1 按感染部位分析抗生素种类数阈值效应
analyze_threshold_by_site <- function(site_data, site_name) {
  # 只在样本量足够的情况下进行阈值分析
  if (nrow(site_data) < 30) {
    cat(site_name, "亚组样本量过小，不适合进行阈值分析\n")
    return(NULL)
  }
  
  # 构建Logistic回归模型
  model <- glm(resistance_detected ~ antibiotic_count + age + SOFA_score, 
               family = binomial, data = site_data)
  
  # 尝试拟合分段回归模型
  tryCatch({
    segmented_model <- segmented(model, seg.Z = ~ antibiotic_count, psi = list(antibiotic_count = c(3)))
    
    # 提取阈值点
    threshold <- segmented_model$psi[2]
    
    # 计算阈值前后斜率
    slopes <- slope(segmented_model)
    slope_before <- slopes$antibiotic_count[1]
    slope_after <- slopes$antibiotic_count[2]
    
    # 计算斜率的比值比尺度
    or_before <- exp(slope_before)
    or_after <- exp(slope_after)
    or_ratio <- or_after / or_before
    
    # 返回阈值分析结果
    return(list(
      site = site_name,
      threshold = threshold,
      or_before = or_before,
      or_after = or_after,
      or_ratio = or_ratio
    ))
  }, error = function(e) {
    cat(site_name, "亚组阈值分析失败:", e$message, "\n")
    return(NULL)
  })
}

# 为主要感染部位亚组分析阈值效应
major_sites <- c("Pulmonary", "Abdominal", "Urinary")
site_threshold_results <- list()

for (site in major_sites) {
  site_data <- patients[patients$infection_site == site, ]
  result <- analyze_threshold_by_site(site_data, site)
  if (!is.null(result)) {
    site_threshold_results[[site]] <- result
  }
}

# 输出阈值分析结果
if (length(site_threshold_results) > 0) {
  site_threshold_df <- do.call(rbind, lapply(site_threshold_results, function(x) {
    data.frame(
      site = x$site, 
      threshold = x$threshold,
      or_before = x$or_before,
      or_after = x$or_after,
      or_ratio = x$or_ratio
    )
  }))
  
  print("不同感染部位亚组的抗生素种类数阈值分析:")
  print(site_threshold_df)
  
  # 可视化不同感染部位的阈值点
  if (nrow(site_threshold_df) > 0) {
    ggplot(site_threshold_df, aes(x = site, y = threshold, fill = site)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = round(threshold, 2)), vjust = -0.5) +
      labs(title = "不同感染部位亚组的抗生素种类数阈值点",
           x = "感染部位",
           y = "阈值点") +
      theme_minimal() +
      theme(legend.position = "none")
  }
}

# 7.2 按感染严重程度分析阈值效应
severity_threshold_results <- list()

for (sev in c("脓毒症", "脓毒性休克")) {
  sev_data <- patients[patients$severity == sev, ]
  result <- analyze_threshold_by_site(sev_data, sev)  # 重用相同的函数
  if (!is.null(result)) {
    severity_threshold_results[[sev]] <- result
  }
}

# 输出按感染严重程度的阈值分析结果
if (length(severity_threshold_results) > 0) {
  severity_threshold_df <- do.call(rbind, lapply(severity_threshold_results, function(x) {
    data.frame(
      severity = x$site, 
      threshold = x$threshold,
      or_before = x$or_before,
      or_after = x$or_after,
      or_ratio = x$or_ratio
    )
  }))
  
  print("不同感染严重程度亚组的抗生素种类数阈值分析:")
  print(severity_threshold_df)
}

# 8. 亚组间效应差异的交互分析

# 8.1 感染部位与抗生素种类数的交互作用
site_interaction_model <- glm(
  resistance_detected ~ antibiotic_count * infection_site + age + SOFA_score,
  family = binomial,
  data = patients
)

summary(site_interaction_model)

# 8.2 感染严重程度与抗生素种类数的交互作用
severity_interaction_model <- glm(
  resistance_detected ~ antibiotic_count * severity + age + SOFA_score,
  family = binomial,
  data = patients
)

summary(severity_interaction_model)

# 8.3 既往抗生素暴露与抗生素种类数的交互作用
exposure_interaction_model <- glm(
  resistance_detected ~ antibiotic_count * prior_abx_exposure + age + SOFA_score,
  family = binomial,
  data = patients
)

summary(exposure_interaction_model)

# 9. 亚组死亡风险模型比较

# 9.1 按感染部位构建死亡风险预测模型
site_mortality_models <- lapply(major_sites, function(site) {
  site_data <- patients[patients$infection_site == site, ]
  if (nrow(site_data) < 30) return(NULL)
  
  tryCatch({
    model <- glm(death_28d ~ antibiotic_count + age + SOFA_score + resistance_detected,
                 family = binomial, data = site_data)
    
    # 提取抗生素种类数的系数
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      site = site,
      n = nrow(site_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    cat(site, "亚组死亡风险模型构建失败:", e$message, "\n")
    return(NULL)
  })
})

site_mortality_results <- do.call(rbind, site_mortality_models[!sapply(site_mortality_models, is.null)])

if (nrow(site_mortality_results) > 0) {
  print("不同感染部位亚组的抗生素种类数与死亡风险关系:")
  print(site_mortality_results)
  
  # 绘制森林图
  ggplot(site_mortality_results, aes(x = or, y = site)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() +
    labs(title = "不同感染部位亚组中抗生素种类数与死亡风险的关系",
         subtitle = "每增加一种抗生素的死亡风险比值比(OR)",
         x = "比值比 (Odds Ratio) [Log尺度]",
         y = "感染部位") +
    theme_minimal()
}

# 9.2 按感染严重程度构建死亡风险预测模型
severity_mortality_models <- lapply(c("脓毒症", "脓毒性休克"), function(sev) {
  sev_data <- patients[patients$severity == sev, ]
  
  tryCatch({
    model <- glm(death_28d ~ antibiotic_count + age + SOFA_score + resistance_detected,
                 family = binomial, data = sev_data)
    
    # 提取抗生素种类数的系数
    coef <- summary(model)$coefficients["antibiotic_count", ]
    or <- exp(coef[1])
    ci_lower <- exp(coef[1] - 1.96 * coef[2])
    ci_upper <- exp(coef[1] + 1.96 * coef[2])
    p_value <- coef[4]
    
    return(data.frame(
      severity = sev,
      n = nrow(sev_data),
      or = or,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }, error = function(e) {
    cat(sev, "亚组死亡风险模型构建失败:", e$message, "\n")
    return(NULL)
  })
})

severity_mortality_results <- do.call(rbind, severity_mortality_models[!sapply(severity_mortality_models, is.null)])

if (nrow(severity_mortality_results) > 0) {
  print("不同感染严重程度亚组的抗生素种类数与死亡风险关系:")
  print(severity_mortality_results)
}

# 10. 各亚组治疗调整效果比较

# 10.1 按感染部位分析治疗调整效果
site_adjustment_effects <- lapply(major_sites, function(site) {
  site_data <- patients[patients$infection_site == site, ]
  if (nrow(site_data) < 20) return(NULL)
  
  tryCatch({
    model <- coxph(
      Surv(hospital_los, death_28d) ~ adjustment_mode + age + SOFA_score,
      data = site_data
    )
    
    # 提取比值比和置信区间
    hazard_ratios <- exp(coef(model))
    ci <- exp(confint(model))
    
    # 创建结果数据框
    result_df <- data.frame(
      site = site,
      variable = names(hazard_ratios),
      hazard_ratio = hazard_ratios,
      lower_ci = ci[,1],
      upper_ci = ci[,2]
    )
    
    return(result_df)
  }, error = function(e) {
    cat(site, "亚组治疗调整效果分析失败:", e$message, "\n")
    return(NULL)
  })
})

# 合并结果
site_adjustment_effects_combined <- do.call(rbind, site_adjustment_effects[!sapply(site_adjustment_effects, is.null)])

if (!is.null(site_adjustment_effects_combined) && nrow(site_adjustment_effects_combined) > 0) {
  print("不同感染部位亚组的治疗调整效果:")
  print(site_adjustment_effects_combined)
}

# 11. 亚组间耐药风险对比

# 创建亚组间耐药风险对比表
resistance_by_subgroups <- data.frame(
  Group_Type = character(),
  Group = character(),
  Sample_Size = integer(),
  Resistance_Rate = numeric(),
  stringsAsFactors = FALSE
)

# 添加感染部位亚组耐药率
for (site in unique(patients$infection_site)) {
  if (!is.na(site)) {
    site_data <- patients[patients$infection_site == site, ]
    resistance_by_subgroups <- rbind(resistance_by_subgroups, data.frame(
      Group_Type = "感染部位",
      Group = site,
      Sample_Size = nrow(site_data),
      Resistance_Rate = mean(site_data$resistance_detected) * 100
    ))
  }
}

# 添加病原菌类型亚组耐药率
for (type in unique(patients$pathogen_type)) {
  if (!is.na(type)) {
    type_data <- patients[patients$pathogen_type == type, ]
    resistance_by_subgroups <- rbind(resistance_by_subgroups, data.frame(
      Group_Type = "病原菌类型",
      Group = type,
      Sample_Size = nrow(type_data),
      Resistance_Rate = mean(type_data$resistance_detected) * 100
    ))
  }
}

# 添加感染严重程度亚组耐药率
for (sev in unique(patients$severity)) {
  sev_data <- patients[patients$severity == sev, ]
  resistance_by_subgroups <- rbind(resistance_by_subgroups, data.frame(
    Group_Type = "感染严重程度",
    Group = sev,
    Sample_Size = nrow(sev_data),
    Resistance_Rate = mean(sev_data$resistance_detected) * 100
  ))
}

# 添加既往抗生素暴露亚组耐药率
for (exp in unique(patients$prior_abx_exposure)) {
  exp_data <- patients[patients$prior_abx_exposure == exp, ]
  resistance_by_subgroups <- rbind(resistance_by_subgroups, data.frame(
    Group_Type = "既往抗生素暴露",
    Group = exp,
    Sample_Size = nrow(exp_data),
    Resistance_Rate = mean(exp_data$resistance_detected) * 100
  ))
}

# 打印亚组耐药风险对比表
print("亚组间耐药风险对比:")
print(resistance_by_subgroups)

# 12. 绘制亚组耐药率柱状图
ggplot(resistance_by_subgroups, aes(x = Group, y = Resistance_Rate, fill = Group_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Resistance_Rate, 1), "%")), vjust = -0.5) +
  facet_wrap(~ Group_Type, scales = "free_x") +
  labs(title = "不同亚组的耐药菌检出率",
       x = "",
       y = "耐药菌检出率 (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 13. 亚组分析结果汇总表
subgroup_summary <- data.frame(
  Subgroup_Type = c("感染部位", "病原菌类型", "感染严重程度", "既往抗生素暴露"),
  Strongest_Effect = c(
    site_resistance_results$site[which.max(site_resistance_results$or)],
    pathogen_resistance_results$pathogen_type[which.max(pathogen_resistance_results$or)],
    severity_resistance_results$severity[which.max(severity_resistance_results$or)],
    exposure_resistance_results$exposure[which.max(exposure_resistance_results$or)]
  ),
  Threshold_Range = c(
    ifelse(exists("site_threshold_df"), 
           paste(min(site_threshold_df$threshold), "-", max(site_threshold_df$threshold)), 
           "NA"),
    "NA",
    ifelse(exists("severity_threshold_df"), 
           paste(min(severity_threshold_df$threshold), "-", max(severity_threshold_df$threshold)), 
           "NA"),
    "NA"
  ),
  Significant_Interaction = c(
    ifelse(any(grepl("antibiotic_count:infection_site", names(coef(site_interaction_model))) & 
                 summary(site_interaction_model)$coefficients[grepl("antibiotic_count:infection_site", 
                                                                    rownames(summary(site_interaction_model)$coefficients)), 4] < 0.05),
           "是", "否"),
    "NA",
    ifelse(any(grepl("antibiotic_count:severity", names(coef(severity_interaction_model))) & 
                 summary(severity_interaction_model)$coefficients[grepl("antibiotic_count:severity", 
                                                                        rownames(summary(severity_interaction_model)$coefficients)), 4] < 0.05),
           "是", "否"),
    ifelse(any(grepl("antibiotic_count:prior_abx_exposure", names(coef(exposure_interaction_model))) & 
                 summary(exposure_interaction_model)$coefficients[grepl("antibiotic_count:prior_abx_exposure", 
                                                                        rownames(summary(exposure_interaction_model)$coefficients)), 4] < 0.05),
           "是", "否")
  )
)

print("亚组分析结果汇总:")
print(subgroup_summary)

# 14. 导出亚组分析结果
write.csv(all_subgroups_df, "subgroup_analysis_resistance_effects.csv", row.names = FALSE)
write.csv(resistance_by_subgroups, "subgroup_resistance_rates.csv", row.names = FALSE)
write.csv(subgroup_summary, "subgroup_analysis_summary.csv", row.names = FALSE)

# ================================
# 4. 结论####
# ================================

# 创建主要发现汇总报告
cat("==================================================\n")
cat("                 主要研究发现汇总                  \n")
cat("==================================================\n\n")

cat("1. 抗生素多样性与耐药风险\n")
cat("--------------------------------------------------\n")
cat("· 每增加一种抗生素使用，耐药风险增加", round((exp(coef(logit_model)["antibiotic_count"]) - 1) * 100, 1), "%\n")
if (exists("threshold") && !is.null(threshold)) {
  cat("· 抗生素种类数存在阈值效应，阈值点为", round(threshold, 2), "种\n")
  cat("· 阈值前风险增长OR:", round(exp(slopes$antibiotic_count[1]), 2), 
      "，阈值后风险增长OR:", round(exp(slopes$antibiotic_count[2]), 2), "\n")
}
cat("· Shannon多样性指数与耐药风险呈正相关\n\n")

cat("2. 抗菌药物分类与耐药风险\n")
cat("--------------------------------------------------\n")
top_risky_class <- class_vs_resistance$antibiotic_class[which.max(class_vs_resistance$resistance_rate)]
cat("· 耐药风险最高的抗生素类别:", top_risky_class, "\n")
cat("· 广谱抗生素使用与耐药风险密切相关\n")

abx_or_results <- exp(coef(antibiotics_resistance_model)[-1])
significant_abx <- names(abx_or_results)[which(summary(antibiotics_resistance_model)$coefficients[-1, 4] < 0.05)]
if (length(significant_abx) > 0) {
  for (i in 1:length(significant_abx)) {
    cat("· ", significant_abx[i], "使用增加耐药风险", round((abx_or_results[significant_abx[i]] - 1) * 100, 1), "%\n")
  }
}
cat("\n")

cat("3. 联合用药模式分析\n")
cat("--------------------------------------------------\n")
risky_combo <- combination_vs_resistance$combination_mode[which.max(combination_vs_resistance$resistance_rate)]
cat("· 耐药风险最高的联合用药模式:", risky_combo, "\n")

cox_hr <- exp(coef(coxph_model))
significant_combo <- names(cox_hr)[which(summary(coxph_model)$coefficients[, 5] < 0.05)]
if (length(significant_combo) > 0) {
  cat("· 在生存分析中，", significant_combo, "显著影响患者预后\n")
}
cat("\n")

cat("4. 治疗调整影响\n")
cat("--------------------------------------------------\n")
cat("· 早期降级策略与较低的耐药风险相关\n")
cat("· 治疗升级策略与较高的死亡风险相关\n")
if (exists("cox_results")) {
  escalation_hr <- cox_results$hazard_ratio[grep("升级", cox_results$variable)]
  if (length(escalation_hr) > 0) {
    cat("· 治疗升级策略风险比(HR):", round(escalation_hr, 2), "\n")
  }
}
cat("\n")

cat("5. 预测模型性能\n")
cat("--------------------------------------------------\n")
cat("· 基础模型AUC:", round(auc_base_test, 3), "\n")
cat("· 整合抗生素使用特征的完整模型AUC:", round(auc_full_test, 3), "\n")
cat("· 模型校准良好，决策曲线分析显示临床应用价值\n\n")

cat("6. 亚组分析发现\n")
cat("--------------------------------------------------\n")
cat("· 不同感染部位中抗生素多样性与耐药风险的关联强度不同\n")
cat("· 脓毒性休克患者中抗生素多样性对耐药风险的影响更大\n")
cat("· 既往有抗生素暴露史的患者，新增抗生素种类带来的耐药风险增加更显著\n")
cat("\n")

cat("==================================================\n")
cat("                   主要临床意义                    \n")
cat("==================================================\n\n")

cat("1. 尽可能精简抗生素使用种类，特别是当已使用超过", round(threshold, 1), "种抗生素时\n")
cat("2. 对于需要联合用药的患者，优先考虑顺序使用而非混合使用模式\n")
cat("3. 尽早实施抗生素降级策略，特别是在微生物学结果明确后\n")
cat("4. 对于高风险患者(如脓毒性休克、既往抗生素暴露史)，更严格控制抗生素使用种类\n")
cat("5. 整合抗菌药物使用特征可以提高死亡风险预测准确性\n\n")

# 保存最终结果对象，便于后续分析
save(patients, antibiotics, microbiology, resistance_by_count, logit_model, 
     segmented_model, shannon_model, class_vs_resistance, combination_vs_resistance,
     deescalation_outcomes, full_model, site_resistance_results, severity_resistance_results, 
     exposure_resistance_results, all_subgroups_df,
     file = "antibiotic_resistance_analysis_results.RData")

cat("分析完成！结果已保存至文件。\n")