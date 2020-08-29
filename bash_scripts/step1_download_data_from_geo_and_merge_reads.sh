#step1_download_data_from_geo_and_merge_reads.sh

BASE_DIR=~/PycharmProjects/pipelines_of_DNAm_kinetics
DATA_DIR=$BASE_DIR/data
mkdir -p $DATA_DIR
# !! **** Download THE DATA FROM Google Drive OR GEO ACCORDING TO THE NOTE, MANUALLY UNZIP AND COPY YOUR DATA INTO ~/PycharmProjects/pipelines_of_DNAm_kinetics/data ****!!

cd $DATA_DIR

#Transform and sort data into target format
#0h_rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' All_S1_S6fraction_1hBrdu_0Chase_nascent_merged_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >0h_rep1.bed &

#1h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182523_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep1_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep1.bed &

#1h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182524_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep2_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep2.bed &

#4h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182525_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep1_processed_sort.cov |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >4h_rep1.bed &

#4h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182526_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep2_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >4h_rep2.bed &

#16h rep1
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182527_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep1_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >16h_rep1.bed &

#16h rep2
awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2+1"\t"$3"\t"$5"\t"$6}' GSM2182528_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep2_processed_sort.cov|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >16h_rep2.bed &

python ../python_code/step2_covert_standard_bed_into_reduced_bed.py -input ./GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3.bed -output ./0h_rep2.bed
python ../python_code/step2_covert_standard_bed_into_reduced_bed.py -input ./GSM2642846_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep3.bed -output ./1h_rep3.bed
python ../python_code/step2_covert_standard_bed_into_reduced_bed.py -input ./GSM2642848_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep3.bed -output ./4h_rep3.bed
python ../python_code/step2_covert_standard_bed_into_reduced_bed.py -input ./GSM2642850_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep3.bed -output ./16h_rep3.bed


TMP_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS/TMP
mkdir -p $TMP_DIR

declare -a TIMEPOINTS=('0h' '1h' '4h' '16h') # 
declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

for tp in "${TIMEPOINTS[@]}"; 
do
  echo $tp
  echo $tp*.bed
  f1=$tp"_rep0.bed"
  cat $tp*.bed  | gsort -k 1,1 -k2,2n --parallel=8  -S 50% | bedtools merge -c 4,5 -o sum>$f1 &
done

for tp in "${TIMEPOINTS[@]}"; 
do
  f1=$tp"_rep0.bed"
  echo $f1
  f1_dir=$tp"_rep0"
  mkdir -p $f1_dir
  for chr in "${CHRS[@]}";
  do
      bedextract $chr $f1 >$f1_dir/$chr.bed &
  done 
done

python  ../python_code/step7_construct_AllDat_matrix.py -data_dir . -out_dir ./AllDat

matlab  -nodesktop -nodisplay -nosplash -r "rate_inference_pipeline('../data/AllDat', '../data/Rates', 1, 22)"


mkdir -p unprocessed
for f in *.cov
do
  mv $f unprocessed
done

#0h rep3
awk 'BEGIN {FS="\t"; OFS=","} {print $6"\t"$7+1"\t"$8"\t"$10"\t"$11}' GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3_processed.cov.tile|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >0h_rep2.bed

#1h rep3
awk 'BEGIN {FS="\t"; OFS=","} {print $6"\t"$7+1"\t"$8"\t"$10"\t"$11}' GSM2642846_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep3_processed.cov.tile|gsort -k 1,1 -k2,2n --parallel=8  -S 50% >1h_rep3.bed

for f in *.tile
do
  mv $f unprocessed
done

rate_base_dir=$DATA_DIR/Rates
non_randomT_rates_dir=$rate_base_dir/Rates
randomT_rates_dir=$rate_base_dir/Rates_RandomT
cd $non_randomT_rates_dir
cat *.bed | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > ../non_randomT_rates.bed &

#Partition into different chromosomes

CHR_DIR=$DATA_DIR/CHROMOSOME_SPLITTED
mkdir -p $CHR_DIR

cd $DATA_DIR

for f in *.bed
do
  filename=${f%%.*}
  echo $filename
  FILE_DIR=$CHR_DIR/$filename
  echo $FILE_DIR
  mkdir -p $FILE_DIR
  cp $f $FILE_DIR
  cd $FILE_DIR
  for chr in `bedextract --list-chr $f`; 
  do
      bedextract $chr $f > $chr.bed; 
  done
  rm $f
 cd $DATA_DIR
done



for f in *.sorted.bed
do
  filename=${f%%.*}
  echo $filename
  #gsort -k 1,1 -k2,2n --parallel=8  -S 50%  $f >$filename.sorted.bed
  #rm $f
  mv $filename.sorted.bed $filename.bed 
done
awk 'BEGIN {FS="\t"; OFS=","; } {print $2"\t"$4}' 

awk 'BEGIN {FS="\t"; OFS=","; } {print $1"\t"$2}' 

