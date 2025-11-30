# Lithology Data Source Comparison for MODFLOW Model

## Background
During the groundwater modeling process, we compared two methods for extracting lithology distribution data from ModMuse GUI to use with FloPy for model calibration.

## Methods Comparison

### Method 1: Shapefile Export
- **Data Source**: KX.shp file exported from ModMuse
- **Processing**: Used Flopy's `GridIntersect` tool to identify grid cells covered by each lithology zone
- **Results**:
  - Layer 1: 429 cells differ from GUI (~2.6% of total grid)
  - Layer 2: 261 cells differ from GUI (~1.6% of total grid)

### Method 2: GPT Project File Parsing
- **Data Source**: Directly parsed ScreenObject polygon definitions from ModMuse's .gpt project file
- **Processing**: Same `GridIntersect` method for cell assignment
- **Results**:
  - Layer 1: 886 cells differ from GUI (~5.4% of total grid)
  - Layer 2: 602 cells differ from GUI (~3.7% of total grid)

## Analysis of Differences

Theoretically, the .gpt file should be the most accurate as it is ModMuse's original data source. However, the Shapefile method actually produces results closer to the GUI. This is due to:

1. **Cell Assignment Algorithm Differences**
   - ModMuse GUI: Likely uses cell-center point method or proprietary algorithms
   - Flopy GridIntersect: Uses area-intersection method
   - Boundary cells are assigned differently between the two approaches

2. **Shapefile as Pre-processed Output**
   - The exported Shapefile may already contain the results of GUI's internal cell-assignment processing
   - Therefore it more accurately reflects the actual lithology distribution used by the GUI

## Conclusion

**Recommendation**: Use the Shapefile method when reproducing ModMuse GUI model results.

Although GPT parsing is theoretically more direct, the Shapefile method demonstrates better consistency with GUI in practice. This study adopted the Shapefile method for parameter calibration.

**Final Model Performance**:
| Condition | RMSE (m) | R² |
|-----------|----------|-----|
| Steady-state | 0.0300 | 0.9458 |
| Pumping | 0.2121 | 0.8421 |

---

# MODFLOW 模型岩性数据源对比分析

## 背景
在地下水建模过程中，我们对比了两种从 ModMuse GUI 提取岩性分布数据的方法，以便配合 FloPy 进行模型校准。

## 方法对比

### 方法一：Shapefile 导出
- **数据源**：从 ModMuse 导出的 KX.shp 文件
- **处理方式**：使用 Flopy 的 `GridIntersect` 工具识别各岩性区覆盖的网格单元
- **结果**：
  - Layer 1：与 GUI 相差 429 个单元格（约占总网格 2.6%）
  - Layer 2：与 GUI 相差 261 个单元格（约占总网格 1.6%）

### 方法二：GPT 项目文件解析
- **数据源**：直接从 ModMuse 的 .gpt 项目文件解析 ScreenObject 多边形定义
- **处理方式**：同样使用 `GridIntersect` 进行单元格分配
- **结果**：
  - Layer 1：与 GUI 相差 886 个单元格（约占总网格 5.4%）
  - Layer 2：与 GUI 相差 602 个单元格（约占总网格 3.7%）

## 差异原因分析

理论上，.gpt 文件是 ModMuse 的原始数据源，应该最准确。然而实际结果显示 Shapefile 方法反而更接近 GUI。原因如下：

1. **单元格分配算法差异**
   - ModMuse GUI：可能使用单元格中心点判定法或其他专有算法
   - Flopy GridIntersect：采用面积交集法
   - 两者对边界单元格的判定存在差异

2. **Shapefile 是预处理结果**
   - ModMuse 导出的 Shapefile 可能已包含 GUI 内部单元格化处理的结果
   - 因此更准确地反映了 GUI 实际使用的岩性分布

## 结论

**建议**：对于需要复现 ModMuse GUI 模型结果的场景，推荐使用 Shapefile 方法。

虽然从理论角度 .gpt 解析更直接，但 Shapefile 方法在实践中与 GUI 的一致性更好。本研究最终采用 Shapefile 方法进行参数校准。

**最终模型性能**：
| 工况 | RMSE (m) | R² |
|------|----------|-----|
| 稳态 | 0.0300 | 0.9458 |
| 抽水 | 0.2121 | 0.8421 |
