# 使用说明

程序与testdata放在不含空格的路径下。testdata文件夹与hgt_demo.py同级。

## 测试

已有Python≥3.10即可运行：

```bash
bash testdata/run_example.sh
```

预期为2个片段对、2行主表，并显示PASS。结果在results_example/05_summary/。重复测试时指定新目录：

```bash
bash testdata/run_example.sh results_example_02
```

testdata是模拟数据，使用预计算的比对、ANI和ORF结果。

## 分析自己的数据

```bash
conda env create -f environment.yml
conda activate hgt-demo

bash hgt_pipeline.sh \
  -i MAGs/ -o results/ \
  --taxonomy taxonomy.tsv \
  -t 8
```

必须提供MAG核苷酸文件和taxonomy。默认使用database/SARG3.2.1-S.fasta，含12,746条蛋白序列。--arg-db可指定其他库；--hgt-only只做DNA片段分析。已有BLASTn、ANI、ORF及ARG比对结果可通过对应参数复用，详见英文README与输入格式说明。

## 数据库缓存

首次运行建立ARG数据库索引，后续自动复用。FASTA内容或比对工具改变时建立对应索引。默认缓存为~/.cache/hgt_demo/arg_databases，可用--arg-db-cache指定目录；--rebuild-arg-db强制重建。各次输出的04_arg_annotation/database.json记录使用的索引。
