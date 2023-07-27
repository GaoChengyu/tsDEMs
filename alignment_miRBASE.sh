#/usr/bin/zsh
#改变参数 更精细
#先与no miRNA的rfam库比对，将比对不上的部分，认为是没有其他ncRNA的reads

#利用botie将过滤ncRNA的测序数据与重复序列数据库比对，剔除其他重复序列
#再与miRBASE库比对，筛选已知的miRNA 这里用前体miRNA库 #U替换为T再比对
#bowtie-build hairpin.fa hairpin  创建bowtie索引
# $1 测序数据所在目录 $2 索引目录 $3 outputfile dir
index_dir=$(cd $(dirname $2);pwd)
index_base=`basename $2`
cd $1
mkdir $3
for i in `ls *fq*`
do
    nohup bowtie -q -p 20 -n 0 -m 1 --best --strata --un ${3}/${i%%.*}.noMIRBASE.fq -x ${index_dir}/${index_base} $i -S ${3}/${i%%.*}.sam 1>${3}/${i%%.*}.log 2>&1 &
done


