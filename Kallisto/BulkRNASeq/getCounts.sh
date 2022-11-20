GTF='gencode.v31.annotation.gtf'
GTF2='gencode.v25.basic.annotation.gtf'

A=$(cat 03_selectSamples_and_featurecounts1.txt)

#featureCounts 
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -t exon -g gene_id -a $GTF -o WolframRerun.txt $(echo $A)
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -t exon -g gene_id -a $GTF2 -o WolframRerun2.txt $(echo $A)
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -f -t transcript -g transcript_id -a $GTF -o WolframRerunTr.txt $(echo $A)
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -f -t transcript -g transcript_id -a $GTF2 -o WolframRerunTr2.txt $(echo $A)
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -t exon -g transcript_id -O -a $GTF -o WolframRerunTr3.txt $(echo $A)
#../../Thorsten//CRUK/subread-2.0.1-Linux-x86_64/bin/featureCounts -p -T 16 -t exon -g transcript_id -O -a $GTF2 -o WolframRerunTr4.txt $(echo $A)

