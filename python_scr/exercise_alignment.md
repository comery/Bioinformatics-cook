
### exercise 1

转换序列坐标

根据agp文件还原原始坐标：

```
superscaffold_chr1.pat	1	1647570	1	W	JAAOMF010000039.1	1	1647570	+
superscaffold_chr1.pat	1647571	1648570	2	N	NA	1	1000	+
superscaffold_chr1.pat	1648571	257535152	3	W	JAAOMF010000001.1	1	255886582	+
superscaffold_chr1.mat	1	147137239	1	W	JAAOME010000002.1	1	147137239	-
superscaffold_chr1.mat	147137240	147138239	2	N	NA	1	1000	+
superscaffold_chr1.mat	147138240	255412214	3	W	JAAOME010000006.1	1	108273975	-
```

'W' WGS scaffold

'N' 代表gap

已知变异位点的坐标，是根据superscaffold_chr1.pat 或者superscaffold_chr1.mat 生成的，现在需要转换成原始scaffold上的坐标，比如 superscaffold_chr1.pat	1对应JAAOMF010000039.1 的1

要求：
需要用两种策略实现，第一种对内存使用没有限制，第二种需要内存使用不能超过1G。

温馨提示：
注意 序列的方向，“-” 表示反向互补序列。

数据：

[test.makeup.agp](./test.makeup.agp)

[test.sv.filterN.bed](./test.sv.filterN.bed)

### exercise 2

有两个基因组版本的比对信息，如下所示：

```
1	122826	525068	402342	122826	122727	99.72	517912	529938	23.72	23.16	chr1-1	ptg000292l
127502	188418	426627	365841	60917	60787	99.77	517912	529938	11.76	11.47	chr1-1	ptg000292l
153892	224008	355482	285366	70117	70117	99.95	517912	529938	13.54	13.23	chr1-1	ptg000292l
```

需要做的是，将chr1-1（ref）上面的对应序列（可以被querry）比对上的部分用querry对应的序列替换掉。

这里需要考虑多个情况：
1. 比如相邻的两个block有overlap的情况；
2. 如何判断一个比对是真正的比对位置，而非是重复序列导致的非共线性比对；
3. ...

所需的数据：
在NCBI上找一个物种的两个不同的组装版本然后进行如下操作得到类似上面的比对信息：

```shell
ref=pecies_v1.fa
asm=species_v2.fa
p="TEST"
mummer='/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/All_kinds_align/mummer4/bin'
$mummer/nucmer --maxmatch -l 100 -c 500 -p $p $ref $asm
$mummer/delta-filter -m -i 90 -l 100 $p.delta > $p.delta.filt
$mummer/dnadiff -d $p.delta -p $p.diff
```
需要根据TEST.diff.1coords 生成比对信息进行处理。


### exercies 3

接练习2，有两个基因组版本的比对信息，如下所示：

```
1	122826	525068	402342	122826	122727	99.72	517912	529938	23.72	23.16	chr1-1	ptg000292l
127502	188418	426627	365841	60917	60787	99.77	517912	529938	11.76	11.47	chr1-1	ptg000292l
153892	224008	355482	285366	70117	70117	99.95	517912	529938	13.54	13.23	chr1-1	ptg000292l
```

需要做的是，将chr1-1（ref）上面的坐标转到querry对应的坐标，类似于liftover的功能。