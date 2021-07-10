# 情景

做生物信息分析，或者其他的一些生物学研究，很多时候我们都需要查找一些基因的信息，而大部分情况都可能会跟人的基因有关，或者最终都要跟人扯上关系，所以如何快速准确的获取基因的信息是一个非常必要掌握的技能。

在众多数据库中，GeneCard是我觉得做得最好也是最全面的数据库，你可以非常方便的输入一个基因名称，然后找到和它所有相关的基因，比如基因与表型间的联系、基因互作蛋白分子、信号通路、临床意义等。

那么问题来了，单个基因很容易通过页面搜索找到，但是有很多基因需要查询的时候该怎么办呢，有没有API接口可以针对数据库查询呢？


# 尝试
## API
找了一圈并没有发现GeneCard的API接口，所以想到一个比较笨的办法就是利用requests 包直接网页端获取页面数据，然后再解析成想要的数据格式。

## 爬虫与反爬虫
其实这个就是一个简单的爬虫程序。
所以基础的代码就很简单：

```python3
def getContentGeneCard(gene):
    api = "https://www.genecards.org/cgi-bin/carddisp.pl?"
    item = "gene={}&keywords={}".format(gene, gene)
    url = api + item

    params = (
        ('gene', gene),
        ('keywords', gene),
    )
    response = requests.get('https://www.genecards.org/cgi-bin/carddisp.pl',  params=params)

    print(url)
    content = response.text
    return content
```

然鹅，事情并没有这么简单，因为像这种网站一般都会有反爬虫的机制。所以一开始你啥也爬不到。如何解决呢？
通过网上查找，发现可以通过获取headers、cookie伪装，就可以很方便的绕过反爬虫了。
具体如何获取heders, 可以参考这篇文章[反 反爬虫](https://blog.csdn.net/qq_43161186/article/details/104382771)

ok,按照这个方法得到了header之后，就可以直接写进程序里面。下面是我的例子，因为每个人用的浏览器和设备不一样，所以需要自行去获得header信息。

![headers](./header.png)

然后稍微修改一下代码，增加headers的参数：

```
response = requests.get('https://www.genecards.org/cgi-bin/carddisp.pl', headers=headers, params=params)
```

然后就可以爬起来了。比如我们爬一下这个基因‘APOB’， 我把它写到一个list里面，按照如下运行程序：

```shell
python3 getInfo_GeneCard.py test.gene.list outdir
```

在outdir 里面有一个命名为APOB的html文件，这个就是我们想要获取的页面信息，但这个html语言的，不是很好看，可以利用buterfly包进行解析

```shell
python3 parse_GeneCard.py test.gene.list outdir
```

得到如下结果：

```text
#### APOB

### Summary
Entrez Gene Summary for APOB Gene
This gene product is the main apolipoprotein of chylomicrons and low density lipoproteins (LDL), and is the ligand for the LDL receptor. It occurs in plasma as two main isoforms, apoB-48 and apoB-100: the former is synthesized exclusively in the gut and the latter in the liver. The intestinal and the hepatic forms of apoB are encoded by a single gene from a single, very long mRNA. The two isoforms share a common N-terminal sequence. The shorter apoB-48 protein is produced after RNA editing of the apoB-100 transcript at residue 2180 (CAA->UAA), resulting in the creation of a stop codon, and early translation termination. Mutations in this gene or its regulatory region cause hypobetalipoproteinemia, normotriglyceridemic hypobetalipoproteinemia, and hypercholesterolemia due to ligand-defective apoB, diseases affecting plasma cholesterol and apoB levels. [provided by RefSeq, Dec 2019]
GeneCards Summary for APOB Gene
APOB (Apolipoprotein B) is a Protein Coding gene.
Diseases associated with APOB include Hypobetalipoproteinemia, Familial, 1 and Hypercholesterolemia, Familial, 2.
Among its related pathways are Signaling by GPCR and Metabolism.
Gene Ontology (GO) annotations related to this gene include binding and heparin binding.
Hypobetalipoproteinemia, Familial, 1
Hypercholesterolemia, Familial, 2
Signaling by GPCR
Metabolism
UniProtKB/Swiss-Prot Summary for APOB Gene
Apolipoprotein B is a major protein constituent of chylomicrons (apo B-48), LDL (apo B-100) and VLDL (apo B-100). Apo B-100 functions as a recognition signal for the cellular binding and internalization of LDL particles by the apoB/E receptor.
APOB_HUMAN,P04114
APOB_HUMAN,P04114
Entrez Gene Summary for APOB Gene
This gene product is the main apolipoprotein of chylomicrons and low density lipoproteins (LDL), and is the ligand for the LDL receptor. It occurs in plasma as two main isoforms, apoB-48 and apoB-100: the former is synthesized exclusively in the gut and the latter in the liver. The intestinal and the hepatic forms of apoB are encoded by a single gene from a single, very long mRNA. The two isoforms share a common N-terminal sequence. The shorter apoB-48 protein is produced after RNA editing of the apoB-100 transcript at residue 2180 (CAA->UAA), resulting in the creation of a stop codon, and early translation termination. Mutations in this gene or its regulatory region cause hypobetalipoproteinemia, normotriglyceridemic hypobetalipoproteinemia, and hypercholesterolemia due to ligand-defective apoB, diseases affecting plasma cholesterol and apoB levels. [provided by RefSeq, Dec 2019]

GeneCards Summary for APOB Gene
APOB (Apolipoprotein B) is a Protein Coding gene.
                                            Diseases associated with APOB include Hypobetalipoproteinemia, Familial, 1 and Hypercholesterolemia, Familial, 2.
                                            Among its related pathways are Signaling by GPCR and Metabolism.
                                            Gene Ontology (GO) annotations related to this gene include binding and heparin binding.
Hypobetalipoproteinemia, Familial, 1
Hypercholesterolemia, Familial, 2
Signaling by GPCR
Metabolism

UniProtKB/Swiss-Prot Summary for APOB Gene
Apolipoprotein B is a major protein constituent of chylomicrons (apo B-48), LDL (apo B-100) and VLDL (apo B-100). Apo B-100 functions as a recognition signal for the cellular binding and internalization of LDL particles by the apoB/E receptor.

                                 APOB_HUMAN,P04114
            ```
            虽然这个格式还有有点丑陋，但不影响我们获得主要信息，如果想要显示的更加整齐美观，可以自行研究下解析的语法和字符串的处理。
            
            

