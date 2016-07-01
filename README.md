# find_position
find the corresponding positions of point (record in BED file, such as DNA methylated sites, SNP or so on) based on the annotation (GTF file). Basically just count the number of reads falling into 5' and 3' UTRs, CDS and introns. 

USE: >perl step1_find_position -r <gtf reference> -f <bed file> -o <optional: the output file name>


So far this script only support gtf file from flybase. If you want to use it for other gtf reference, please check the code and make sure you capture the right name of the transcript id. 

When the files were generated from step1, you can use step2_sum.R to figure out whether the distribution and density of those points may be correlated or not with CBI, the codon usage bias index. The code are better run under R studio.   

The CBI reference file attached was made with drosophila r6.07. 
