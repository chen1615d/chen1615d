# 1.数据生成####
# 加载必要的包
rm(list = ls())
library(tidyverse)
library(lubridate)
library(pROC)
library(rms)
library(segmented)
library(ggplot2)
library(mgcv)
library(tableone)
library(glmnet)
library(survival)
library(survminer)

# 设置随机种子确保可重复性
set.seed(2025)

# 生成模拟患者基本信息数据
generate_patient_data <- function(n = 200) {
  # 患者ID
  patient_id <- paste0("P", sprintf("%03d", 1:n))
  
  # 基本特征
  age <- round(rnorm(n, 65, 12))
  gender <- sample(c("男", "女"), n, replace = TRUE, prob = c(0.6, 0.4))
  weight <- round(rnorm(n, 65, 10))
  height <- round(rnorm(n, 165, 8))
  bmi <- round(weight / ((height/100)^2), 1)
  
  # 基础疾病 (0=无, 1=有)
  hypertension <- rbinom(n, 1, 0.45)
  diabetes <- rbinom(n, 1, 0.25)
  heart_disease <- rbinom(n, 1, 0.20)
  renal_disease <- rbinom(n, 1, 0.15)
  liver_disease <- rbinom(n, 1, 0.10)
  copd <- rbinom(n, 1, 0.15)
  cancer <- rbinom(n, 1, 0.08)
  
  # 住院信息
  department <- sample(c("急诊医学科", "呼吸内科", "普通外科", "神经内科", "心内科"), 
                       n, replace = TRUE)
  los_total <- round(rlnorm(n, log(14), 0.5))
  los_icu <- round(rlnorm(n, log(7), 0.6))
  mechanical_ventilation <- rbinom(n, 1, 0.6)
  vasopressor_use <- rbinom(n, 1, 0.5)
  
  # 生理参数
  sofa_score <- round(runif(n, 2, 15))
  apache_ii <- round(rnorm(n, 18, 5))
  temperature_admission <- round(rnorm(n, 38.2, 1.1), 1)
  map_admission <- round(rnorm(n, 75, 15))
  hr_admission <- round(rnorm(n, 105, 20))
  rr_admission <- round(rnorm(n, 24, 6))
  
  # 实验室检查
  wbc_d1 <- round(rlnorm(n, log(12), 0.5), 1)
  neutrophil_pct_d1 <- round(runif(n, 65, 95), 1)
  lymphocyte_pct_d1 <- round(100 - neutrophil_pct_d1 - runif(n, 2, 8), 1)
  crp_d1 <- round(rlnorm(n, log(120), 0.8))
  procalcitonin_d1 <- round(rlnorm(n, log(5), 1.2), 2)
  
  # 感染部位
  infection_lung <- rbinom(n, 1, 0.45)
  infection_abdomen <- rbinom(n, 1, 0.25)
  infection_urinary <- rbinom(n, 1, 0.15)
  infection_bloodstream <- rbinom(n, 1, 0.25)
  infection_skin <- rbinom(n, 1, 0.08)
  infection_catheter <- rbinom(n, 1, 0.12)
  infection_multiple <- ifelse(rowSums(cbind(infection_lung, infection_abdomen, 
                                             infection_urinary, infection_bloodstream, 
                                             infection_skin, infection_catheter)) > 1, 1, 0)
  
  # 创建数据框
  patient_data <- data.frame(
    patient_id = patient_id,
    age = age,
    gender = gender,
    weight = weight,
    height = height,
    bmi = bmi,
    hypertension = hypertension,
    diabetes = diabetes,
    heart_disease = heart_disease,
    renal_disease = renal_disease,
    liver_disease = liver_disease,
    copd = copd,
    cancer = cancer,
    department = department,
    los_total = los_total,
    los_icu = los_icu,
    mechanical_ventilation = mechanical_ventilation,
    vasopressor_use = vasopressor_use,
    sofa_score = sofa_score,
    apache_ii = apache_ii,
    temperature_admission = temperature_admission,
    map_admission = map_admission,
    hr_admission = hr_admission,
    rr_admission = rr_admission,
    wbc_d1 = wbc_d1,
    neutrophil_pct_d1 = neutrophil_pct_d1,
    lymphocyte_pct_d1 = lymphocyte_pct_d1,
    crp_d1 = crp_d1,
    procalcitonin_d1 = procalcitonin_d1,
    infection_lung = infection_lung,
    infection_abdomen = infection_abdomen,
    infection_urinary = infection_urinary,
    infection_bloodstream = infection_bloodstream,
    infection_skin = infection_skin,
    infection_catheter = infection_catheter,
    infection_multiple = infection_multiple
  )
  
  return(patient_data)
}

# 生成抗生素使用数据
generate_antibiotics_data <- function(patient_data) {
  # 定义抗生素列表及其分类
  antibiotics <- list(
    "青霉素类" = c("哌拉西林他唑巴坦", "阿莫西林克拉维酸钾", "氨苄西林舒巴坦"),
    "头孢菌素类" = c("头孢曲松", "头孢吡肟", "头孢他啶", "头孢噻肟"),
    "碳青霉烯类" = c("美罗培南", "亚胺培南", "厄他培南"),
    "氨基糖苷类" = c("阿米卡星", "庆大霉素", "妥布霉素"),
    "喹诺酮类" = c("左氧氟沙星", "环丙沙星", "莫西沙星"),
    "糖肽类" = c("万古霉素", "替考拉宁"),
    "其他" = c("利奈唑胺", "替加环素", "复方磺胺甲恶唑", "甲硝唑")
  )
  
  # 展平抗生素列表
  all_antibiotics <- unlist(antibiotics)
  
  # 获取患者ID
  patient_ids <- patient_data$patient_id
  
  # 为每个患者随机分配1-8种抗生素
  antibiotics_data <- list()
  
  for (id in patient_ids) {
    # 获取该患者的SOFA评分作为风险参考
    sofa_score <- patient_data$sofa_score[patient_data$patient_id == id]
    
    # 根据SOFA评分调整抗生素使用数量概率
    prob_vector <- c(0.05, 0.15, 0.20, 0.25, 0.15, 0.10, 0.06, 0.04)
    if (sofa_score > 10) {
      # 严重患者更可能使用多种抗生素
      prob_vector <- c(0.02, 0.08, 0.15, 0.25, 0.20, 0.15, 0.10, 0.05)
    }
    
    # 决定该患者使用的抗生素数量
    num_antibiotics <- sample(1:8, 1, prob = prob_vector)
    
    # 选择抗生素
    selected_antibiotics <- sample(all_antibiotics, num_antibiotics)
    
    # 为每种抗生素生成使用信息
    for (abx in selected_antibiotics) {
      # 确定抗生素类别
      abx_class <- names(antibiotics)[sapply(antibiotics, function(x) abx %in% x)]
      
      # 决定是否为特殊级抗生素
      is_special <- abx %in% c("美罗培南", "亚胺培南", "万古霉素", "利奈唑胺", "替加环素")
      
      # 使用天数: 特殊级抗生素一般使用时间短一些
      days_of_use <- if(is_special) round(runif(1, 3, 10)) else round(runif(1, 5, 14))
      
      # 开始时间: 入院后0-3天内开始
      start_day <- round(runif(1, 0, 3), 1)
      
      # 结束时间
      end_day <- start_day + days_of_use
      
      # 给药途径
      route <- sample(c("静脉", "口服"), 1, prob = c(0.8, 0.2))
      if (abx %in% c("万古霉素", "美罗培南", "亚胺培南", "利奈唑胺")) {
        route <- "静脉"  # 这些药物通常静脉给药
      }
      
      # 日剂量 (DDD)
      ddd <- round(runif(1, 0.8, 1.5), 2)
      
      # 将该条记录添加到列表
      antibiotics_data[[length(antibiotics_data) + 1]] <- data.frame(
        patient_id = id,
        antibiotic_name = abx,
        antibiotic_class = abx_class,
        is_special = as.integer(is_special),
        start_day = start_day,
        end_day = end_day,
        days_of_use = days_of_use,
        route = route,
        ddd = ddd
      )
    }
  }
  
  # 合并所有记录
  antibiotics_df <- do.call(rbind, antibiotics_data)
  
  return(antibiotics_df)
}

# 生成耐药监测数据
generate_resistance_data <- function(patient_data, antibiotics_data) {
  # 汇总每个患者的抗生素使用情况
  patient_abx_summary <- antibiotics_data %>%
    group_by(patient_id) %>%
    dplyr::summarize(
      num_antibiotics = n_distinct(antibiotic_name),
      num_special = sum(is_special),
      total_days = sum(days_of_use),
      abx_classes = n_distinct(antibiotic_class)
    )
  
  # 合并患者基本信息
  patient_combined <- patient_data %>%
    left_join(patient_abx_summary, by = "patient_id")
  
  # 针对每个患者生成耐药监测结果
  resistance_data <- list()
  
  for (i in 1:nrow(patient_combined)) {
    patient <- patient_combined[i, ]
    
    # 决定是否检出耐药菌
    # 假设抗生素种类数、特殊级抗生素使用和SOFA评分会影响耐药风险
    
    # 计算耐药风险基础概率 (注意：这里我们设计一个阈值效应）
    # 当抗生素数量为1-6时，耐药概率随数量增加而增加
    # 当抗生素数量>6时，增长显著减缓（阈值效应）
    base_prob <- case_when(
      patient$num_antibiotics == 1 ~ 0.05,
      patient$num_antibiotics == 2 ~ 0.10,
      patient$num_antibiotics == 3 ~ 0.20,
      patient$num_antibiotics == 4 ~ 0.30,
      patient$num_antibiotics == 5 ~ 0.45,
      patient$num_antibiotics == 6 ~ 0.55, # 阈值点
      patient$num_antibiotics == 7 ~ 0.60,
      patient$num_antibiotics >= 8 ~ 0.65
    )
    
    # 调整因素
    adj_factor <- 1.0
    if (patient$sofa_score > 10) adj_factor <- adj_factor * 1.3
    if (patient$num_special >= 2) adj_factor <- adj_factor * 1.2
    if (patient$mechanical_ventilation == 1) adj_factor <- adj_factor * 1.2
    if (patient$infection_multiple == 1) adj_factor <- adj_factor * 1.3
    
    # 最终耐药风险
    resistance_prob <- min(base_prob * adj_factor, 0.95)
    
    # 决定是否有耐药菌
    has_resistance <- rbinom(1, 1, resistance_prob)
    
    if (has_resistance == 1) {
      # 决定耐药菌类型
      pathogen_types <- c(
        "大肠埃希菌-ESBL", "肺炎克雷伯菌-KPC", "铜绿假单胞菌-MDR", 
        "鲍曼不动杆菌-CR", "金黄色葡萄球菌-MRSA"
      )
      
      # 根据感染部位调整病原菌概率
      pathogen_probs <- c(0.25, 0.20, 0.20, 0.20, 0.15)
      
      if (patient$infection_lung == 1) {
        pathogen_probs <- c(0.10, 0.25, 0.25, 0.25, 0.15)
      } else if (patient$infection_abdomen == 1) {
        pathogen_probs <- c(0.40, 0.30, 0.10, 0.10, 0.10)
      } else if (patient$infection_urinary == 1) {
        pathogen_probs <- c(0.50, 0.20, 0.15, 0.05, 0.10)
      }
      
      pathogen <- sample(pathogen_types, 1, prob = pathogen_probs)
      
      # 决定耐药级别
      resistance_levels <- c("MDR", "XDR", "PDR")
      resistance_level <- sample(resistance_levels, 1, 
                                 prob = c(0.7, 0.25, 0.05))
      
      # 决定检出时间
      detection_day <- round(runif(1, 2, min(10, patient$los_total)), 1)
      
      # 将该记录添加到列表
      resistance_data[[length(resistance_data) + 1]] <- data.frame(
        patient_id = patient$patient_id,
        has_resistance = has_resistance,
        pathogen = pathogen,
        resistance_level = resistance_level,
        detection_day = detection_day
      )
    } else {
      # 无耐药菌的记录
      resistance_data[[length(resistance_data) + 1]] <- data.frame(
        patient_id = patient$patient_id,
        has_resistance = has_resistance,
        pathogen = NA,
        resistance_level = NA,
        detection_day = NA
      )
    }
  }
  
  # 合并所有记录
  resistance_df <- do.call(rbind, resistance_data)
  
  return(resistance_df)
}

# 生成临床结局数据 - 修复后的版本
generate_outcomes <- function(patient_data, antibiotics_data, resistance_data) {
  # 合并所有数据
  antibiotics_summary <- antibiotics_data %>%
    group_by(patient_id) %>%
    dplyr::summarize(
      num_antibiotics = n_distinct(antibiotic_name),
      num_special = sum(is_special),
      total_abx_days = sum(days_of_use),
      abx_classes = n_distinct(antibiotic_class),
      first_abx_time = min(start_day)
    )
  
  combined_data <- patient_data %>%
    left_join(antibiotics_summary, by = "patient_id") %>%
    left_join(resistance_data, by = "patient_id")
  
  # 模拟死亡结局
  # 计算基础死亡风险
  combined_data$base_mortality_risk <- with(combined_data, {
    # 年龄风险
    age_risk <- (age - 60) / 30
    
    # SOFA评分风险 
    sofa_risk <- sofa_score / 10
    
    # 抗生素相关风险（设计阈值效应）
    abx_risk <- case_when(
      num_antibiotics <= 2 ~ -0.1,
      num_antibiotics == 3 ~ 0,
      num_antibiotics == 4 ~ 0.1,
      num_antibiotics == 5 ~ 0.2,
      num_antibiotics == 6 ~ 0.3,
      num_antibiotics >= 7 ~ 0.4
    )
    
    # 首次给药时间风险（越晚风险越高）
    time_risk <- first_abx_time * 0.1
    
    # 耐药风险
    resistance_risk <- ifelse(has_resistance == 1, 0.3, 0)
    
    # 整合风险
    integrated_risk <- 0.1 + age_risk + sofa_risk + abx_risk + time_risk + resistance_risk
    
    # 其他调整因素 - 修复部分：使用向量化操作替代if语句
    integrated_risk <- integrated_risk + 0.2 * (mechanical_ventilation == 1)
    integrated_risk <- integrated_risk + 0.2 * (vasopressor_use == 1)
    integrated_risk <- integrated_risk + 0.15 * (infection_multiple == 1)
    
    # 确保风险在合理范围内
    pmin(pmax(integrated_risk, 0.05), 0.95)
  })
  
  # 生成28天死亡结局
  set.seed(42) # 确保可重复性
  combined_data$mortality_28d <- rbinom(nrow(combined_data), 1, combined_data$base_mortality_risk)
  
  # 生成生存时间
  combined_data$survival_time <- with(combined_data, {
    if_else(mortality_28d == 1,
            round(rexp(length(mortality_28d), 1/15) * 28), # 死亡患者的生存时间
            28) # 存活患者设为28天
  })
  
  # 删除中间变量
  combined_data$base_mortality_risk <- NULL
  
  return(combined_data)
}

# 生成模拟数据
set.seed(2025)
patient_data <- generate_patient_data(200)
antibiotics_data <- generate_antibiotics_data(patient_data)
resistance_data <- generate_resistance_data(patient_data, antibiotics_data)
final_data <- generate_outcomes(patient_data, antibiotics_data, resistance_data)

# 抗生素多样性指数计算
calculate_adi <- function(antibiotics_data) {
  # 按患者计算抗生素种类和使用天数
  abx_diversity <- antibiotics_data %>%
    group_by(patient_id, antibiotic_name) %>%
    dplyr::summarize(days = sum(days_of_use), .groups = "drop") %>%
    group_by(patient_id) %>%
    mutate(
      total_days = sum(days),
      proportion = days / total_days
    ) %>%
    dplyr::summarize(
      ADI = -sum(proportion * log(proportion)),
      .groups = "drop"
    )
  
  return(abx_diversity)
}

# 计算抗生素多样性指数
adi_data <- calculate_adi(antibiotics_data)

# 将ADI合并到最终数据中
final_data <- final_data %>%
  left_join(adi_data, by = "patient_id")

# 保存数据以便后续分析
saveRDS(final_data, "sepsis_final_data.rds")
saveRDS(antibiotics_data, "sepsis_antibiotics_data.rds")

# 2.描述性分析####
# 加载保存的数据
final_data <- readRDS("sepsis_final_data.rds")
antibiotics_data <- readRDS("sepsis_antibiotics_data.rds")

# 创建抗生素分组变量
final_data$abx_group <- cut(final_data$num_antibiotics, 
                            breaks = c(0, 1, 2, 5, Inf),
                            labels = c("单一抗生素", "2种抗生素", "3-5种抗生素", ">5种抗生素"),
                            right = FALSE)

# 创建表格1：基线特征按抗生素使用分组比较
vars_to_include <- c("age", "gender", "bmi", "hypertension", "diabetes", "heart_disease", 
                     "sofa_score", "apache_ii", "mechanical_ventilation", "vasopressor_use",
                     "los_icu", "infection_lung", "infection_abdomen", "infection_bloodstream",
                     "infection_multiple", "has_resistance", "mortality_28d")

categorical_vars <- c("gender", "hypertension", "diabetes", "heart_disease", 
                      "mechanical_ventilation", "vasopressor_use", "infection_lung",
                      "infection_abdomen", "infection_bloodstream", "infection_multiple",
                      "has_resistance", "mortality_28d")

table1 <- CreateTableOne(vars = vars_to_include, 
                         strata = "abx_group", 
                         data = final_data, 
                         factorVars = categorical_vars)

# 打印表格1
table1_print <- print(table1, showAllLevels = TRUE, printToggle = FALSE, nonnormal = c("sofa_score", "los_icu"))
table1_print

# 分析抗生素使用模式
abx_patterns <- antibiotics_data %>%
  group_by(antibiotic_class) %>%
  dplyr::summarize(
    n_patients = n_distinct(patient_id),
    avg_duration = mean(days_of_use, na.rm = TRUE),
    pct_special = mean(is_special) * 100
  ) %>%
  arrange(desc(n_patients))

print(abx_patterns)

# 分析最常见的抗生素组合（前5种）
top_combinations <- antibiotics_data %>%
  group_by(patient_id) %>%
  dplyr::summarize(combination = paste(sort(unique(antibiotic_name)), collapse = " + ")) %>%
  count(combination) %>%
  arrange(desc(n)) %>%
  head(5)

print(top_combinations)

# 可视化抗生素种类数分布
ggplot(final_data, aes(x = num_antibiotics)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "抗生素种类数分布", x = "抗生素种类数", y = "患者数") +
  theme_minimal()

# 可视化抗生素种类数与耐药率的关系
resistance_by_abx <- final_data %>%
  group_by(num_antibiotics) %>%
  dplyr::summarize(
    n_patients = n(),
    n_resistant = sum(has_resistance),
    resistance_rate = mean(has_resistance) * 100
  )

ggplot(resistance_by_abx, aes(x = num_antibiotics, y = resistance_rate)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "抗生素种类数与耐药率关系", 
       x = "抗生素种类数", 
       y = "耐药率 (%)") +
  theme_minimal()

# 可视化抗生素种类数与28天死亡率的关系
mortality_by_abx <- final_data %>%
  group_by(num_antibiotics) %>%
  dplyr::summarize(
    n_patients = n(),
    n_deaths = sum(mortality_28d),
    mortality_rate = mean(mortality_28d) * 100
  )

ggplot(mortality_by_abx, aes(x = num_antibiotics, y = mortality_rate)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "抗生素种类数与28天死亡率关系", 
       x = "抗生素种类数", 
       y = "28天死亡率 (%)") +
  theme_minimal()

# 3. 抗生素多样性与耐药风险阈值效应分析####
# 阈值效应分析

# 1. 分段线性回归寻找阈值
# 首先构建基础模型
resistance_model <- glm(has_resistance ~ num_antibiotics + age + sofa_score + 
                          mechanical_ventilation + infection_multiple, 
                        family = binomial, data = final_data)

# 使用分段回归寻找阈值点
segmented_model <- segmented(resistance_model, seg.Z = ~num_antibiotics, psi = c(5))
summary(segmented_model)

# 提取断点(阈值)
breakpoint <- summary(segmented_model)$psi[2]

# 生成预测数据用于绘图
new_data <- expand.grid(
  num_antibiotics = seq(1, 8, by = 0.1),
  age = mean(final_data$age),
  sofa_score = mean(final_data$sofa_score),
  mechanical_ventilation = 0,
  infection_multiple = 0
)

# 预测值
pred <- predict(segmented_model, newdata = new_data, type = "response")
new_data$predicted_risk <- pred

# 绘制阈值效应图
ggplot(new_data, aes(x = num_antibiotics, y = predicted_risk)) +
  geom_line(size = 1.2, color = "blue") +
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red") +
  annotate("text", x = breakpoint + 0.5, y = 0.7, 
           label = paste("阈值 =", round(breakpoint, 1)), color = "red") +
  labs(title = "抗生素多样性与耐药风险的阈值效应", 
       x = "抗生素种类数", 
       y = "预测耐药风险") +
  theme_minimal()

# 2. 限制性立方样条函数分析
library(rms)
dd <- datadist(final_data)
options(datadist = "dd")

# 使用RCS样条函数拟合非线性关系
spline_model <- lrm(has_resistance ~ rcs(num_antibiotics, 4) + age + 
                      sofa_score + mechanical_ventilation + infection_multiple, 
                    data = final_data)

# 输出模型摘要
print(spline_model)

# 检验非线性关系显著性
anova(spline_model)

# 绘制调整后的剂量-反应曲线
plot_spline <- Predict(spline_model, num_antibiotics, 
                       fun = plogis,  # 转换为概率
                       conf.int = 0.95)

ggplot(as.data.frame(plot_spline), aes(x = num_antibiotics, y = yhat)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(title = "抗生素多样性与耐药风险的非线性关系", 
       x = "抗生素种类数", 
       y = "预测耐药风险") +
  theme_minimal()


# 3. 广义加性模型 - 修复版本
gam_model <- gam(has_resistance ~ s(num_antibiotics, k=5) + age + 
                   sofa_score + mechanical_ventilation + infection_multiple, 
                 family = binomial, data = final_data)
# 输出模型摘要
summary(gam_model)

# 绘制平滑曲线
plot(gam_model, select = 1, shade = TRUE, xlab = "抗生素种类数", 
     ylab = "s(抗生素种类数)", main = "广义加性模型: 抗生素数量的平滑效应")

# 计算不同抗生素种类数下的耐药风险
abx_groups <- data.frame(num_antibiotics = 1:8) %>%
  mutate(
    group_label = paste(num_antibiotics, "种抗生素"),
    n_patients = sapply(num_antibiotics, function(x) sum(final_data$num_antibiotics == x)),
    n_resistant = sapply(num_antibiotics, function(x) sum(final_data$has_resistance[final_data$num_antibiotics == x])),
    resistance_rate = n_resistant / n_patients * 100
  )

# 可视化阈值效应
ggplot(abx_groups, aes(x = num_antibiotics, y = resistance_rate)) +
  geom_point(aes(size = n_patients), color = "blue") +
  geom_line() +
  geom_vline(xintercept = round(breakpoint), linetype = "dashed", color = "red") +
  labs(title = "抗生素种类数与耐药率关系", 
       subtitle = paste("阈值点:", round(breakpoint, 1), "种抗生素"),
       x = "抗生素种类数", 
       y = "耐药率 (%)",
       size = "患者数") +
  theme_minimal()

# 计算阈值前后的斜率变化
slope_before <- diff(abx_groups$resistance_rate[1:round(breakpoint)]) / 
  diff(abx_groups$num_antibiotics[1:round(breakpoint)])
slope_after <- diff(abx_groups$resistance_rate[round(breakpoint):8]) / 
  diff(abx_groups$num_antibiotics[round(breakpoint):8])

cat("阈值前平均斜率:", mean(slope_before, na.rm = TRUE), "\n")
cat("阈值后平均斜率:", mean(slope_after, na.rm = TRUE), "\n")
cat("斜率变化率:", (mean(slope_after, na.rm = TRUE) / mean(slope_before, na.rm = TRUE) - 1) * 100, "%\n")

# 4. 多因素分析与预测模型构建####
# 多因素分析：抗生素使用与临床结局

# 1. 耐药风险多因素Logistic回归
resistance_multiv <- glm(has_resistance ~ num_antibiotics + age + sofa_score + 
                           mechanical_ventilation + vasopressor_use + infection_multiple +
                           num_special + first_abx_time, 
                         family = binomial, data = final_data)

# 输出结果
summary(resistance_multiv)
odds_ratios <- exp(coef(resistance_multiv))
ci <- exp(confint(resistance_multiv))
resistance_results <- data.frame(
  OR = odds_ratios,
  Lower_CI = ci[, 1],
  Upper_CI = ci[, 2],
  P_value = summary(resistance_multiv)$coefficients[, 4]
)

print(resistance_results)

# 2. 死亡风险多因素Logistic回归
mortality_multiv <- glm(mortality_28d ~ num_antibiotics + age + sofa_score + 
                          mechanical_ventilation + vasopressor_use + infection_multiple +
                          num_special + first_abx_time + has_resistance, 
                        family = binomial, data = final_data)

# 输出结果
summary(mortality_multiv)
odds_ratios <- exp(coef(mortality_multiv))
ci <- exp(confint(mortality_multiv))
mortality_results <- data.frame(
  OR = odds_ratios,
  Lower_CI = ci[, 1],
  Upper_CI = ci[, 2],
  P_value = summary(mortality_multiv)$coefficients[, 4]
)

print(mortality_results)

# 使用分类变量进行分析：抗生素使用分组
final_data$abx_category <- cut(final_data$num_antibiotics, 
                               breaks = c(0, 1, 2, 5, Inf), 
                               labels = c("1种", "2种", "3-5种", ">5种"))

mortality_cat <- glm(mortality_28d ~ abx_category + age + sofa_score + 
                       mechanical_ventilation + vasopressor_use + infection_multiple +
                       num_special + first_abx_time + has_resistance, 
                     family = binomial, data = final_data)

# 输出结果
summary(mortality_cat)
odds_ratios <- exp(coef(mortality_cat))
ci <- exp(confint(mortality_cat))
mortality_cat_results <- data.frame(
  OR = odds_ratios,
  Lower_CI = ci[, 1],
  Upper_CI = ci[, 2],
  P_value = summary(mortality_cat)$coefficients[, 4]
)

print(mortality_cat_results)

# 3. LASSO回归进行变量筛选
# 准备数据
x_vars <- model.matrix(~ num_antibiotics + age + gender + bmi + hypertension + 
                         diabetes + heart_disease + renal_disease + liver_disease + 
                         copd + cancer + sofa_score + apache_ii + 
                         mechanical_ventilation + vasopressor_use + 
                         infection_lung + infection_abdomen + infection_urinary + 
                         infection_bloodstream + infection_multiple + 
                         num_special + first_abx_time + has_resistance + ADI - 1, 
                       data = final_data)

y_mortality <- final_data$mortality_28d

# 执行LASSO回归
set.seed(42)
lasso_model <- cv.glmnet(x_vars, y_mortality, family = "binomial", alpha = 1)

# 查看最佳lambda值
best_lambda <- lasso_model$lambda.1se
cat("最佳lambda值:", best_lambda, "\n")

# 查看选定的变量
selected_vars <- which(coef(lasso_model, s = best_lambda) != 0)
selected_names <- colnames(x_vars)[selected_vars]
cat("LASSO选择的变量:", paste(selected_names, collapse = ", "), "\n")

# 基于LASSO结果构建死亡预测模型
formula_str <- paste("mortality_28d ~", paste(selected_names, collapse = " + "))
final_model <- glm(as.formula(formula_str), family = binomial, data = final_data)

# 输出最终模型结果
summary(final_model)

# 4. 构建基础模型与抗生素模型比较
# 基础模型：仅包含临床特征
base_model <- glm(mortality_28d ~ age + sofa_score + mechanical_ventilation + 
                    vasopressor_use + infection_multiple, 
                  family = binomial, data = final_data)

# 抗生素模型：增加抗生素使用特征
abx_model <- glm(mortality_28d ~ age + sofa_score + mechanical_ventilation + 
                   vasopressor_use + infection_multiple + num_antibiotics + 
                   num_special + first_abx_time + has_resistance, 
                 family = binomial, data = final_data)

# 比较模型性能
base_pred <- predict(base_model, type = "response")
abx_pred <- predict(abx_model, type = "response")

# ROC曲线比较
roc_base <- roc(final_data$mortality_28d, base_pred)
roc_abx <- roc(final_data$mortality_28d, abx_pred)

# 绘制ROC曲线
plot(roc_base, col = "blue", main = "ROC曲线比较")
lines(roc_abx, col = "red")
legend("bottomright", 
       legend = c(paste("基础模型 AUC =", round(auc(roc_base), 3)),
                  paste("抗生素模型 AUC =", round(auc(roc_abx), 3))),
       col = c("blue", "red"), 
       lwd = 2)

# 比较两个模型的AUC
roc_test <- roc.test(roc_base, roc_abx)
cat("AUC比较P值:", roc_test$p.value, "\n")

# 5. 基于RMS包创建列线图
library(rms)
dd <- datadist(final_data)
options(datadist = "dd")

# 转换到RMS格式
final_model_rms <- lrm(as.formula(formula_str), data = final_data)

# 创建列线图
nomogram <- nomogram(final_model_rms, 
                     fun = plogis,  # 转换为概率
                     lp = FALSE,    # 不显示线性预测值
                     funlabel = "28天死亡概率")

plot(nomogram)

# 6. Bootstrap内部验证
set.seed(42)
validate_stats <- validate(final_model_rms, method = "boot", B = 100)
print(validate_stats)

# 计算校准斜率和截距
cal <- calibrate(final_model_rms, method = "boot", B = 100)
plot(cal, main = "校准曲线")
# 5. 耐药风险分层系统与决策支持工具####
# 基于模型构建耐药风险分层系统

# 使用最终模型预测耐药风险
final_data$predicted_resistance <- predict(resistance_multiv, type = "response")

# 定义风险分层
final_data$resistance_risk_level <- cut(final_data$predicted_resistance,
                                        breaks = c(0, 0.1, 0.3, 1),
                                        labels = c("低风险", "中风险", "高风险"))

# 分析不同风险分层的特征
risk_summary <- final_data %>%
  group_by(resistance_risk_level) %>%
  dplyr::summarize(
    n_patients = n(),
    avg_antibiotics = mean(num_antibiotics),
    avg_special = mean(num_special),
    avg_sofa = mean(sofa_score),
    actual_resistance_rate = mean(has_resistance) * 100,
    mortality_rate = mean(mortality_28d) * 100
  )

print(risk_summary)

# 可视化不同风险分层的耐药率和死亡率
risk_long <- risk_summary %>%
  select(resistance_risk_level, actual_resistance_rate, mortality_rate) %>%
  pivot_longer(cols = c(actual_resistance_rate, mortality_rate),
               names_to = "outcome",
               values_to = "rate")

ggplot(risk_long, aes(x = resistance_risk_level, y = rate, fill = outcome)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "耐药风险分层与临床结局", 
       x = "耐药风险分层", 
       y = "比率 (%)") +
  scale_fill_discrete(name = "结局", 
                      labels = c("实际耐药率", "28天死亡率")) +
  theme_minimal()

# 构建抗生素使用决策支持工具

# 基于阈值效应的抗生素数量建议
threshold_value <- round(breakpoint)

# 为不同风险患者提供抗生素使用建议
risk_recommendations <- data.frame(
  risk_level = c("低风险", "中风险", "高风险"),
  max_recommended_abx = c(2, threshold_value, threshold_value),
  strategy = c(
    "优先单一窄谱抗生素，治疗5-7天，尽早降阶梯",
    paste0("合理联合用药，建议≤", threshold_value, "种，定期评估继续/停用指征"),
    paste0("精准联合治疗，注意控制在", threshold_value, "种以内，考虑循环使用策略")
  )
)

print(risk_recommendations)

# 可视化阈值指导的干预策略
patient_distribution <- final_data %>%
  count(resistance_risk_level, num_antibiotics) %>%
  group_by(resistance_risk_level) %>%
  mutate(pct = n / sum(n) * 100)

# 添加阈值线与建议区域
ggplot(patient_distribution, aes(x = num_antibiotics, y = pct, fill = resistance_risk_level)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_vline(xintercept = threshold_value, linetype = "dashed", color = "red") +
  annotate("text", x = threshold_value + 0.3, y = max(patient_distribution$pct), 
           label = paste("阈值 =", threshold_value), color = "red") +
  annotate("rect", xmin = -Inf, xmax = threshold_value, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "green") +
  annotate("rect", xmin = threshold_value, xmax = Inf, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "red") +
  labs(title = "基于耐药风险分层的抗生素使用分布", 
       subtitle = paste("阈值效应指导：建议控制在", threshold_value, "种抗生素以内"),
       x = "抗生素种类数", 
       y = "患者比例 (%)") +
  facet_wrap(~ resistance_risk_level) +
  theme_minimal()

# 决策流程图（文本描述）
cat("==== 脓毒血症抗生素使用决策流程 ====\n\n")
cat("1. 评估患者耐药风险\n")
cat("   - 低风险: 无特殊基础疾病、SOFA<5、无多部位感染\n")
cat("   - 中风险: 有基础疾病、SOFA 5-10、单一部位感染\n")
cat("   - 高风险: 多种基础疾病、SOFA>10、多部位感染或既往耐药菌感染\n\n")
cat("2. 初始抗生素选择\n")
cat("   - 低风险: 单一经验性抗生素\n")
cat("   - 中风险: 2-3种抗生素联合\n")
cat("   - 高风险: 3-", threshold_value, "种抗生素精准联合\n\n")
cat("3. 监测与调整\n")
cat("   - 获得药敏结果后及时调整\n")
cat("   - 密切关注总抗生素种类数，避免超过", threshold_value, "种\n")
cat("   - 如必须使用>", threshold_value, "种抗生素，加强耐药监测\n\n")
cat("4. 降阶梯策略\n")
cat("   - 低风险: 48-72小时评估降阶梯\n")
cat("   - 中风险: 96小时评估降阶梯\n")
cat("   - 高风险: 持续评估，条件允许时优先减少抗生素种类\n\n")
cat("5. 疗程优化\n")
cat("   - 低风险: 5-7天为主\n")
cat("   - 中风险: 7-10天\n")
cat("   - 高风险: 10-14天，个体化评估\n\n")
