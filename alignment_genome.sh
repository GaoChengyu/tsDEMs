#/usr/bin/zsh
#过滤完的数据与比对不到mirbase的与基因组比对
#bowtie-build genome.fa genome  创建bowtie索引
# $1 测序数据所在目录 $2 索引目录 $3 outputfile dir
index_dir=$(cd $(dirname $2);pwd)
index_base=`basename $2`
cd $1
mkdir $3
for i in `ls *fa*`
do
    nohup bowtie -f -p 20 -n 0 -m 1 --best --strata -x ${index_dir}/${index_base} $i -S ${3}/${i%%.*}.sam 1>${3}/${i%%.*}.log 2>&1 &
done