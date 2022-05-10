# MaxEntPy ANNOVAR

基于[MaxEntPy](https://github.com/kepbod/maxentpy)和ANNVOAR分析结果，进行MaxEnt分析
   
## 参数说明

- -i/--annovar: 输入文件，ANNOVAR输出结果
- -g/--genome: 基因组fasta
- -r/--refgene: refGene文件
- -c/--gene_detail: ANNOVAR结果中基因详情坐在列，一般包括ANNOVAR GeneDetail、AAChange等包含转录本信息列
- -o/--output: 输出文件

## 安装运行
```shell
zhuy@ubuntu:/projects/example$ pip install maxentpy_annovar
zhuy@ubuntu:/projects/example$ maxentpy-annovar.py -i test.split.annovar.txt -o test.maxent.txt -g /path/to/hg19.fa -r /path/to/refgene.txt -c Detail
```

### 结果说明
1. Chr: 染色体
2. Start: 变异起始位置
3. End: 变异终止位置
4. Ref：Ref
5. Alt: Alt
6. Maxent_type:
6. Maxent_pred:
7. Maxent_score_ref:
8. Maxent_score_alt:
9. Maxent_score_var:
10. Maxent_foldchange: