# IDEA 中文入门

[English project overview](../README.md)

IDEA 将“直接测量哪些 SNP”和“预测某个性状时使用哪些 SNP”分开处理。
IDEA-design 面向基因型填充进行稀疏面板选择；填充后，IDEA-GP 按性状选择预测标记。

## 先验证安装和小示例

在仓库根目录执行：

```bash
python -m pip install -e ".[test]"
bash scripts/run_toy_idea_design.sh
pytest -q
```

小示例输出为 `results/toy/toy_selected_snps.txt`：

```text
1_200
1_300
1_500
```

示例的六个候选 SNP、权重和校准汇总均为人工构造，只用于检查程序运行。
它们不含 UK Biobank 或企业猪个体记录，也不是论文实际使用的校准参数。

## 使用自己的数据

| 目的 | 入口 | 说明 |
| --- | --- | --- |
| 从预处理特征生成面板 | `idea-design run --select-only` | 通过 `--target-snp-total` 指定跨染色体的总标记数 |
| 运行小示例所用的选择器 | `idea-design select` | 按物种和群体目录读取输入；`--num-bins` 指位置分箱数 |
| 运行填充与评估 | `idea-design run --config configs/idea.example.yaml` | 先替换参考基因型、基准数据、JAR 和工具路径 |
| 按性状进行预测 | `idea-gp run` | 提供 PLINK 二进制数据前缀及含表头的 `FID IID Pheno` 表型文件 |

详细参数和输入格式见 [IDEA-design](idea_design.md)、[IDEA-GP](idea_gp.md) 和
[安装及外部依赖](installation.md)。网页工具
[alphaindex.zju.edu.cn/IDEA](https://alphaindex.zju.edu.cn/IDEA)
支持 IDEA-design 面板生成。

## 仓库与论文复现材料

仓库目前提供核心代码、说明文档、配置模板、小示例和测试。
论文实际校准输入、固定权重、最终 SNP 清单及结果分析脚本仍需与最终稿匹配后归档，
其状态列于[复现记录](manuscript_reproducibility.md)。
UK Biobank 和企业猪个体数据通过各自数据管理方授权获取，不在此仓库分发；
公共参考数据入口见[数据获取说明](data_availability.md)。

报告程序问题时，请使用 [GitHub Issues](https://github.com/nariew111/IDEA/issues)，
附上命令、代码提交版本、工具版本和脱敏后的报错信息。
