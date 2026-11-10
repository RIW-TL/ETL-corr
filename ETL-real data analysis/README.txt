NHANES n0=100、全部 source：两图精简复现包

运行：在 RStudio 中打开 MAIN.R，点击 Source 即可；也可运行 source("MAIN.R")。
仅使用基础 R，无须安装额外包或下载数据。生成文件在 results 文件夹：
  Figure6_fullsource.pdf
  Figure7_components.pdf
  summary_metrics.csv（由100次误差重新计算的均值、SD、SE）

MAIN.R 是唯一程序；data 文件夹仅有3份必要输入：
  correlation_pairs.csv：1431对变量的两个 source 相关系数差异；
  variance_coordinates.csv：54个变量的信息及两个 source/target 方差比；
  split_errors.csv：4种方法在100次划分中的3项误差，共400行。

本包由已保存分析结果重现两张图，不从原始观测重新拟合 ETL。
均值和标准误由每次划分的误差重新计算，误差棒为±1.96 SE。
第二张图保持13.2×4.2英寸、三位小数；右图仅画共同方差估计一次。
既有A/M选择使用相同测试样本矩，属于exploratory/oracle结果。
完整原始数据及拟合程序仍在此前交付的完整压缩包中。

数据背景：NHANES 2017–March 2020；54个变量；target为60岁及以上。
每次target训练100人，测试2242人；两个source分别为20–39岁1999人、
40–59岁2218人，均使用全部source。Figure6对应第一次训练划分。
数据来源：https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020
