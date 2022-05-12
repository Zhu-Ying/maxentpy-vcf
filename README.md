# MaxEntPy ANNOVAR

基于[MaxEntPy](https://github.com/kepbod/maxentpy)和SNV VCF，进行MaxEnt分析
   
## 参数说明

- -i/--vcf: 输入文件，VCF文件
- -g/--genome: 基因组fasta
- -r/--refgene: refGene文件
- -o/--output: 输出文件

## 安装运行
```shell
zhuy@ubuntu:/projects/example$ pip install maxentpy_vcf
zhuy@ubuntu:/projects/example$ maxentpy-vcf.py -i test.vcf -o test.maxent.txt -g /path/to/hg19.fa -r /path/to/refgene.txt
```

### 结果说明
1. Chr: 染色体
2. Start: 变异起始位置
3. End: 变异终止位置
4. Ref: REF
5. Alt: ALT
6. Maxent_type: MaxEnt Splicing 类型（donor/acceptor）
6. Maxent_pred: MaxEnt 预测预测类型（T/D/U）
7. Maxent_score_ref: MaxEnt REF 得分
8. Maxent_score_alt: MaxEnt ALT 得分
9. Maxent_score_var: MaxEnt 综合得分
10. Maxent_foldchange: MaxEnt ALT/REF 得分比值