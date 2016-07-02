# find_position
find the corresponding positions of point (record in BED file, such as DNA methylated sites, SNP or so on) based on the annotation (GTF file). Basically just count the number of reads falling into 5' and 3' UTRs, CDS and introns. Then statistics can be analyzed with R codes. 

### USE:
```perl step1_find_position -r <file.gtf> -f <file.bed> -o <outputName>```

### NOTE: 
1. bed file should hold single site position. 
2. So far this script only tested for gtf file from <b>flybase</b>. If you want to use it for other gtf reference, please check the code and make sure you capture the right name of the transcript id. 
3. outputName is optional. You don't have to assign a name. But better have it or you may get confused in the future :).

When the files were generated from step1, you can use step2_sum.R to figure out whether the distribution and density of those points may be correlated or not with CBI, the codon usage bias index. The code are better run under R studio.   

The CBI reference file attached was made with drosophila r6.07. 
