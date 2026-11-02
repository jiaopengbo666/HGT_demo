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
  --arg-db /path/to/SARG.fasta \
  -t 8
```

必须提供MAG核苷酸文件和taxonomy。ARG蛋白库由用户提供。已有BLASTn、ANI、ORF及ARG比对结果可通过对应参数复用，详见英文README与输入格式说明。
