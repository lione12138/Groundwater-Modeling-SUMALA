# MODFLOW 6 Groundwater Flow Modeling - SUMALA Aquifer System
# MODFLOW 6 地下水流动模型 - SUMALA 含水层系统

---

## English Version

### Project Overview

This project presents a comprehensive groundwater flow modeling study of the SUMALA aquifer system using MODFLOW 6 (v6.6.2). The model simulates both steady-state and transient groundwater flow conditions, incorporating multiple lithological units and boundary conditions.

### Project Structure

```
SUMALA_Groundwater_Model/
│
├── README.md                          # Project documentation (bilingual)
├── requirements.txt                   # Python dependencies
├── .gitignore                         # Git ignore rules
│
├── 01_Steady_Calibrate/               # 📊 Steady-State Calibration
│   ├── Code/                          # Python calibration scripts
│   │   ├── 01_lhs_calibration.py      # LHS-based batch calibration
│   │   ├── 02_lhs_calibration_gpt.py  # Calibration using .gpt parsing
│   │   ├── 03_single_run_shp.py       # Single-run testing (shapefile)
│   │   ├── 04_single_run_gpt.py       # Single-run testing (.gpt file)
│   │   ├── 05_sensitivity_analysis.py # OAT sensitivity analysis
│   │   ├── 06_method_comparison.py    # Compare different methods
│   │   ├── 07_plot_lithology.py       # Lithology distribution plot
│   │   └── 08_plot_results.py         # Results visualization
│   ├── Model/                         # Model input files (.gpt, .dis, etc.)
│   ├── Lithology/                     # Lithology shapefiles
│   └── Out/                           # Calibration output results
│       ├── lhs_runs/                  # LHS calibration results
│       ├── run_flopy_shp/             # FloPy shapefile runs
│       ├── run_flopy_gpt/             # FloPy .gpt file runs
│       ├── run_modelmuse/             # ModelMuse GUI runs
│       └── Sensitivity_Analysis/      # Sensitivity analysis results
│
├── 02_Transient_Calibrate/            # ⏱️ Transient Calibration
│   ├── Code/                          # Python transient calibration scripts
│   │   ├── 01_lhs_calibration.py      # LHS-based transient calibration
│   │   ├── 02_single_run.py           # Single-run transient testing
│   │   ├── 03_sensitivity_analysis.py # Transient sensitivity analysis
│   │   ├── 04_bayesian_optimization.py# Bayesian optimization
│   │   ├── 05_correlation_analysis.py # Parameter correlation analysis
│   │   └── 06_plot_results.py         # Interpolation comparison plots
│   ├── Model/                         # Transient model files
│   ├── Correlation/                   # Correlation analysis results
│   └── Output/                        # Simulation outputs
│       ├── lhs_runs/                  # LHS calibration results
│       ├── run_test/                  # Single test runs
│       ├── run_modelmuse/             # ModelMuse GUI runs
│       └── Sensitivity_Analysis/      # Sensitivity analysis results
│
├── 03_Dewatering/                     # 🏗️ Dewatering Design
│   ├── Minimum_Rate/                  # 💧 Minimum pumping rate (steady-state)
│   │   └── [MODFLOW model files]      # Steady-state dewatering model
│   └── Minimum_Time/                  # ⏰ Minimum pumping time (transient)
│       ├── Plot_Dewatering_Results.py # Dewatering visualization
│       └── [MODFLOW model files]      # Transient dewatering model
│
├── Config/                            # ⚙️ Configuration & Reference Data
│   ├── Lithology_Parameters.xlsx      # Lithology hydraulic parameters
│   ├── borehole_logs.xlsx             # Borehole log data
│   ├── Initial_Head.txt               # Initial head configuration
│   ├── Final_Head.txt                 # Final head configuration
│   └── Observation_Points.txt         # Observation point coordinates
│
├── Instructions/                      # 📚 Course Materials & Observations
│   ├── Code/                          # Template visualization scripts
│   ├── *.csv                          # Calibration observation data
│   └── *.pdf                          # Course documentation
│
├── Report/                            # 📝 Project Report & Documentation
└── Reference/                         # 📖 Reference Materials
```

### Model Specifications

| Component | Description |
|-----------|-------------|
| **Model Type** | MODFLOW 6 (v6.6.2) |
| **Grid** | Structured grid (125 columns × 125 rows × 2 layers) |
| **Cell Size** | 25m × 25m |
| **Domain** | 3.125 km × 3.125 km |
| **Layers** | Layer 1: Unconfined (sandy/gravel); Layer 2: Confined (limestone/sandstone) |

### Boundary Conditions

- **CHD (Constant Head)**: Fixed heads at domain boundaries
- **GHB (General Head Boundary)**: Northern and southern boundaries
- **RIV (River)**: Central river feature
- **WEL (Well)**: Pumping well (P1) at -1500 m³/d
- **RCH (Recharge)**: Spatially variable recharge rates

### Calibration Methodology

#### Steady-State Calibration
- **Method**: Latin Hypercube Sampling (LHS)
- **Parameters**: 10 parameters (6 Kx values + initial head + RIV/GHB conductance + recharge factor)
- **Samples**: 200-500 combinations
- **Target**: Minimize RMSE, maximize R² for 12 observation points

#### Transient Calibration
- **Additional Parameters**: Specific storage (Ss) and Specific yield (Sy) for each lithology
- **Observation Points**: P1 and Pz12 drawdown data
- **Methods**: LHS, Bayesian Optimization

### Key Results

| Condition | RMSE (m) | R² |
|-----------|----------|-----|
| Natural (SS) | <0.15 | >0.98 |
| Pumping (SS) | <0.15 | >0.98 |
| Transient | TBD | TBD |

### Requirements

- Python 3.10+
- FloPy (≥3.4.0)
- NumPy, Pandas, SciPy
- Matplotlib, Seaborn
- Fiona, Shapely (for shapefile processing)
- MODFLOW 6 executable

### Installation

```bash
pip install flopy numpy pandas scipy matplotlib seaborn fiona shapely
```

### Usage

1. **Steady-State Calibration**:
   ```bash
   cd 01_Steady_Calibrate/Code
   python 01_lhs_calibration.py
   ```

2. **Transient Calibration**:
   ```bash
   cd 02_Transient_Calibrate/Code
   python 01_lhs_calibration.py
   ```

3. **Dewatering Analysis**:
   ```bash
   cd 03_Dewatering/Minimum_Time
   python Plot_Dewatering_Results.py
   ```

4. **Visualization**:
   ```bash
   python 08_plot_results.py
   ```

### Authors

- Graduate Student, Ghent University
- Course: Ground Water Modelling (2024-2025)

### License

This project is for educational purposes as part of the Ground Water Modelling course at Ghent University.

---

## 中文版本

### 项目概述

本项目使用 MODFLOW 6 (v6.6.2) 对 SUMALA 含水层系统进行全面的地下水流动模拟研究。模型模拟了稳态和瞬态地下水流动条件，包含多个岩性单元和边界条件。

### 项目结构

```
SUMALA_Groundwater_Model/
│
├── README.md                          # 项目说明（双语）
├── requirements.txt                   # Python 依赖
├── .gitignore                         # Git 忽略规则
│
├── 01_Steady_Calibrate/               # 📊 稳态校准
│   ├── Code/                          # Python 校准脚本
│   │   ├── 01_lhs_calibration.py      # 基于 LHS 的批量校准
│   │   ├── 02_lhs_calibration_gpt.py  # 使用 .gpt 解析的校准
│   │   ├── 03_single_run_shp.py       # 单次测试（shapefile）
│   │   ├── 04_single_run_gpt.py       # 单次测试（.gpt 文件）
│   │   ├── 05_sensitivity_analysis.py # OAT 敏感性分析
│   │   ├── 06_method_comparison.py    # 不同方法对比
│   │   ├── 07_plot_lithology.py       # 岩性分布图
│   │   └── 08_plot_results.py         # 结果可视化
│   ├── Model/                         # 模型输入文件（.gpt, .dis 等）
│   ├── Lithology/                     # 岩性 shapefile 文件
│   └── Out/                           # 校准输出结果
│       ├── lhs_runs/                  # LHS 校准结果
│       ├── run_flopy_shp/             # FloPy shapefile 运行
│       ├── run_flopy_gpt/             # FloPy .gpt 文件运行
│       ├── run_modelmuse/             # ModelMuse GUI 运行
│       └── Sensitivity_Analysis/      # 敏感性分析结果
│
├── 02_Transient_Calibrate/            # ⏱️ 瞬态校准
│   ├── Code/                          # Python 瞬态校准脚本
│   │   ├── 01_lhs_calibration.py      # 基于 LHS 的瞬态校准
│   │   ├── 02_single_run.py           # 瞬态单次测试
│   │   ├── 03_sensitivity_analysis.py # 瞬态敏感性分析
│   │   ├── 04_bayesian_optimization.py# 贝叶斯优化
│   │   ├── 05_correlation_analysis.py # 参数相关性分析
│   │   └── 06_plot_results.py         # 插值对比图
│   ├── Model/                         # 瞬态模型文件
│   ├── Correlation/                   # 相关性分析结果
│   └── Output/                        # 模拟输出
│       ├── lhs_runs/                  # LHS 校准结果
│       ├── run_test/                  # 单次测试运行
│       ├── run_modelmuse/             # ModelMuse GUI 运行
│       └── Sensitivity_Analysis/      # 敏感性分析结果
│
├── 03_Dewatering/                     # 🏗️ 降水设计
│   ├── Minimum_Rate/                  # 💧 最小抽水率（稳态）
│   │   └── [MODFLOW 模型文件]         # 稳态降水模型
│   └── Minimum_Time/                  # ⏰ 最小抽水时间（瞬态）
│       ├── Plot_Dewatering_Results.py # 降水可视化
│       └── [MODFLOW 模型文件]         # 瞬态降水模型
│
├── Config/                            # ⚙️ 配置与参考数据
│   ├── Lithology_Parameters.xlsx      # 岩性水力参数
│   ├── borehole_logs.xlsx             # 钻孔记录数据
│   ├── Initial_Head.txt               # 初始水头配置
│   ├── Final_Head.txt                 # 最终水头配置
│   └── Observation_Points.txt         # 观测点坐标
│
├── Instructions/                      # 📚 课程材料与观测数据
│   ├── Code/                          # 模板可视化脚本
│   ├── *.csv                          # 校准观测数据
│   └── *.pdf                          # 课程文档
│
├── Report/                            # 📝 项目报告与文档
└── Reference/                         # 📖 参考资料
```

### 模型规格

| 组件 | 描述 |
|------|------|
| **模型类型** | MODFLOW 6 (v6.6.2) |
| **网格** | 结构化网格 (125列 × 125行 × 2层) |
| **单元大小** | 25m × 25m |
| **模型域** | 3.125 km × 3.125 km |
| **分层** | 第1层：非承压层（砂/砾石）；第2层：承压层（石灰岩/砂岩） |

### 边界条件

- **CHD（定水头边界）**：模型域边界固定水头
- **GHB（一般水头边界）**：北部和南部边界
- **RIV（河流边界）**：中部河流特征
- **WEL（井）**：抽水井 (P1)，流量 -1500 m³/d
- **RCH（入渗补给）**：空间变化的入渗率

### 校准方法

#### 稳态校准
- **方法**：拉丁超立方采样 (LHS)
- **参数**：10 个参数（6 个渗透系数 Kx + 初始水头 + RIV/GHB 导水度 + 入渗因子）
- **采样数**：200-500 组组合
- **目标**：最小化 RMSE，最大化 12 个观测点的 R²

#### 瞬态校准
- **附加参数**：各岩性的比储存系数 (Ss) 和给水度 (Sy)
- **观测点**：P1 和 Pz12 的降深数据
- **方法**：LHS、贝叶斯优化

### 主要结果

| 条件 | RMSE (m) | R² |
|------|----------|-----|
| 天然状态（稳态） | <0.15 | >0.98 |
| 抽水状态（稳态） | <0.15 | >0.98 |
| 瞬态 | 待定 | 待定 |

### 环境要求

- Python 3.10+
- FloPy (≥3.4.0)
- NumPy, Pandas, SciPy
- Matplotlib, Seaborn
- Fiona, Shapely（用于 shapefile 处理）
- MODFLOW 6 可执行文件

### 安装

```bash
pip install flopy numpy pandas scipy matplotlib seaborn fiona shapely
```

### 使用方法

1. **稳态校准**：
   ```bash
   cd 01_Steady_Calibrate/Code
   python 01_lhs_calibration.py
   ```

2. **瞬态校准**：
   ```bash
   cd 02_Transient_Calibrate/Code
   python 01_lhs_calibration.py
   ```

3. **降水分析**：
   ```bash
   cd 03_Dewatering/Minimum_Time
   python Plot_Dewatering_Results.py
   ```

4. **可视化**：
   ```bash
   python 08_plot_results.py
   ```

### 作者

- 根特大学研究生
- 课程：地下水模拟 (2024-2025)

### 许可证

本项目仅用于根特大学地下水模拟课程的教学目的。

---

## GitHub Repository

🔗 [https://github.com/lione12138/Groundwate_Modeling_SUMALA](https://github.com/lione12138/Groundwate_Modeling_SUMALA)
