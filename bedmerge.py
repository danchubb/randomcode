
import pandas as pd
import os
import sys

beds=sys.argv[1]
merged_file=sys.argv[2]

merged_file=merged_file.rstrip(".bed")

#beds=open("/data/scratch/DGE/DUDGE/MOPOPGEN/dchubb/bed_file_merge/bedlist.txt")
#merged_file="/data/scratch/DGE/DUDGE/MOPOPGEN/dchubb/bed_file_merge/merged.bed"

bedtools="/opt/software/applications/bedtools/2.29.2/bin/bedtools"
bedops="/data/scratch/DGE/DUDGE/MOPOPGEN/dchubb/soft/bin/bedops"

#beds="/data/scratch/DGE/DUDGE/MOPOPGEN/dchubb/bed_file_merge/bedlist.txt"

BEDS=open(beds)

bedlists=[]

for i in BEDS.readlines():
    i=i.rstrip()
    b=i.split('/')[-1].rstrip('.bed')
    bedlists.append((b,i))

make=bedops+' --merge'

for i,g in bedlists:
	bedfile=[]
	merge_file=i+".bed"
	merge_sort=bedtools+" sort -i "+g+" | "+bedtools+" merge -i stdin > "+i+"_mergesort.bed"
	os.system(merge_sort)
	print(merge_sort)
	make=make+" "+i+"_mergesort.bed"
make=make+" > "+merged_file+".bed"
print(make)
os.system(make)

intersect="cat "+merge_file 
for i,g in bedlists:
    intersect=intersect+"| "+bedtools+" intersect -a stdin -b "+i+"_mergesort.bed"+" -c "
intersect=intersect+" > "+merged_file+"_support.bed"
os.system(intersect)

allnames=[a[0] for a in bedlists]
table_headers=['chr','start','end']+allnames
combined=pd.read_csv(merged_file+"_support.bed",names=table_headers,sep='\t')
combined['instances']=combined[[a[0] for a in bedlists]].ge(1).agg(sum,axis=1)
combined.to_csv(merged_file+"_final.tsv",sep="\t",index=False)

import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pl
import numpy as np
from sklearn.decomposition import PCA
#x=pd.read_csv("test_final.tsv",sep="\t")
y=combined[allnames].transpose()
pca = PCA(n_components=5)
pca.fit(y.ge(1))
#pca.fit(x[['HT29_Merged_elements', 'CACO2_Merged_elements', 'SW480_Merged_elements', 'HCT116_Merged_elements', 'HCA7_Merged_elements', 'C32_Merged_elements', 'DLD1_Merged_elements', 'CL11_Merged_elements', 'SW948_Merged_elements', 'SW403_Merged_elements', 'instances']].ge(1))
columns = ['pca_%i' % i for i in range(5)]
df_pca = pd.DataFrame(pca.transform(y.ge(1)), columns=columns, index=y.index)
df_pca.reset_index(inplace=True)
df_pca = df_pca.rename(columns = {'index':'samples'})
groups = df.groupby('samples')
fig, ax = pl.subplots()
for name, group in groups: 
    ax.plot(group.pca_0, group.pca_1, marker='o', linestyle='', ms=12, label=name)
ax.legend(bbox_to_anchor=(1,1), loc="upper left")
fig.savefig('pca.pdf',bbox_inches='tight')


#>>> gt=[]
#>>> count=[]
#>>> for i in range(10):
#...  gt.append(i)
#...  count.append((len(combined[combined.instances>i])))
#.
#>>> from matplotlib.ticker import MaxNLocator
#>>> ax.xaxis.set_major_locator(MaxNLocator(10))
#>>> fig=ax.get_figure()
#>>> fig.savefig('running_total.pdf')
#>>> ax=df.plot.line(x='gt',y='count')
#>>> ax.xaxis.set_major_locator(MaxNLocator(10))
#>>> ax.yaxis.set_major_locator(MaxNLocator(10))
#>>> fig=ax.get_figure()
#>>> fig.savefig('running_total.pdf')



#sigs['SampleID']=sigs.Samples.apply(lambda x: x.split('_')[-2]+"_"+x.split('_')[-1])
#heads=list(sigs)

#heads.remove('Samples')
#heads.remove('SampleID')
#sigs_pc=sigs
#sigs_pc[heads]=sigs[heads].div(sigs[heads].sum(axis=1), axis=0)
#classifications=pd.read_csv('/re_gecip/cancer_colorectal/analysisResults/13.driverAnalysis/OncoKB_annotation/CRC_v8_2023_rough_oncotree_classification.tsv',sep="\t")
#sigs_pc=pd.merge(sigs_pc,classifications,left_on='SampleID',right_on='SAMPLE_ID')
#sigs_pc=sigs_pc.rename(columns={"ONCOTREE_CODE": "subtype"})
#del sigs_pc['Samples']
#sigs_pc=sigs_pc[['SampleID','subtype']+heads]
#samples=pd.read_csv("MSS_samples.txt",names=['SampleID'])
#
#sigs_pc=pd.merge(sigs_pc,samples,on='SampleID')
#
#sigs_pc.to_csv("CRC_sample_cov.tsv")
