# 最终论文图的15张小鼠原图与分析代码

本项目对应作者确认的最终小鼠拼图：**2、4、8 h 三行，五个实验组五列，共15张完整原图；每张图只分析最终拼图显示的3只小鼠，共45个ROI结果。**

五个组按图中从左到右排列：mRNA、DOPC-mRNA、DOPC-Chole-mRNA、66%DOPC33%DOPG-mRNA、66%DOPC33%DOPG-Chole-mRNA。作者确认实验为5组、每组3只小鼠，同一批动物在2、4、8 h接受成像，每次FFz注射后5分钟成像。

## 项目内容

| 文件 | 用途 |
|---|---|
| `analyze_mice.py` | 针对这15张图的分析程序 |
| `data/images/` | 15张与最终拼图逐格核对的完整TIFF，原像素和色条保留 |
| `data/manifest.json` | 每张图的文件校验值、色条上下限、3只鼠的ROI坐标和图中位置 |
| `data/image_index.csv` | 时间、组别、原文件名和项目内图片路径 |
| `requirements.txt` | 固定版本的运行依赖 |
| `results/` | 本项目实际运行生成的45个ROI数值及核查结果 |
| `tests/` | 图像计算规则的测试 |

发布范围固定为以上15张TIFF。代码逐项读取清单，只计算清单指定的3个ROI。输出数值直接来自这些原图；本项目不会读取历史作图数值进行匹配，也不宣称旧版图表数值已被复现。

## 运行

使用Python 3.12，在项目目录运行：

```sh
python -m pip install -r requirements.txt
python analyze_mice.py
```

结果默认写入 `outputs/`。运行不修改原图。进行科学规则测试：

```sh
python -m unittest discover -s tests -v
```

## 图像定量方法

每只鼠使用清单内固定的矩形ROI，纵向范围为 `y=100:780`；横向坐标沿用与该鼠位置对应的历史逐鼠边界，坐标上限不包含在ROI内。只分析最终图中显示的鼠，不在运行时重新筛选。

程序在同一张图中检测色条，对色条每行取RGB中位数并每隔一行抽样。ROI内候选像素须满足 `max(R,G,B)-min(R,G,B)>25` 且 `max(R,G,B)>55`。将候选像素与色条RGB做欧氏距离最近匹配，距离不大于55的像素纳入；按该图色条上下限线性赋值，然后对纳入像素求和。未做背景扣除、ROI面积归一化或跨图归一化。

该指标称为 **image-derived integrated bioluminescence signal (a.u.)**，即图像推算的积分生物发光信号。原图是带色条的RGB导出图，像素求和不是仪器原生、经物理校准的总光子通量。零表示ROI内没有通过上述规则的显示信号像素。

原图中的左右位置用于定位ROI，不能自动视为跨时间的动物ID。此程序执行图像定量，不执行组间检验或重复测量统计。

## 版本对应

原图均按最终截图的鼠体姿势、足尾及发光位置和颜色核对。4H的025、027使用嵌套目录中的重导出版本；8H的040使用带“(2)”文件名的更新版；8H的043使用与最终处理裁图一致的更新版；8H最后一组使用050的较早导出版本及图中第1、2、3只鼠。实际文件校验值和ROI以清单为准。

## Manuscript method

Bioluminescence was quantified from the 15 full exported pseudocolor images corresponding to the final figure. Three predefined mouse-specific rectangular regions of interest were analyzed per image. Chromatic pixels were identified by RGB channel differences and matched to the nearest RGB color on the color bar of the same image. Accepted pixels were assigned values by linear mapping between the displayed color-bar limits, and these values were summed within each ROI. The result is reported as image-derived integrated bioluminescence signal in arbitrary units. No background subtraction, ROI-area normalization or cross-image normalization was applied. A zero indicates that no displayed pixels passed the inclusion criteria. The supplied manifest specifies the image hashes, scale limits and ROI coordinates for all 45 measurements.

## 上传

已于2026-10-02检查这15张原始TIFF：每张为1.888–1.999 MiB，总计29.009 MiB（30,418,146字节）。GitHub网页允许单文件不超过25 MiB、每次最多100个文件，因此15张原图可直接上传，不需要缩图、降低质量或改格式。依据：[GitHub官方上传限制](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)。

解压后把项目内容放到GitHub仓库根目录，保持文件夹结构。ZIP用于整体交付，仓库内按各个文件上传。此文件夹是本地准备的上传材料，没有自动上传或发布功能。
