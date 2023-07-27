#/usr/bin/zsh
#mkdir -p ./Rfam/fasta_files
####下载Rfam数据
#wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/*
#进入http://rfam.xfam.org/search#tabview=tab5，选择microRNA，submit，复制，粘贴保存在Rfam/rfam_miRNAall.txt
cd Rfam
awk 'BEGIN{FS = "\t";OFS="\t"}{print $1}' rfam_miRNAall.txt > rfam_mirid.txt
ls fasta_files > rfam_id.txt
mkdir rfam_nmi
grep -v -f rfam_mirid.txt rfam_id.txt | xargs -I {} cp fasta_files/{} rfam_nmi/
gunzip -r rfam_nmi
cat rfam_nmi/* >rfam_noMIRNA.fa
#bowtie-build rfam_noMIRNA.fa rfam_noMIRNA  创建bowtie索引
