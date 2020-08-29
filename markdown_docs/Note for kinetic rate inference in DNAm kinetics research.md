# Note for kinetic rate inference in DNAm kinetics research
## 1. Prepare data
All reads data are collected from Dr. Meissner A’s Lab. Some are available from [GSE82045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82045), see detail below:

### 1.1 Data used
Download reads data from shared [Google Drive](https://drive.google.com/drive/folders/1J41xMGs3H4SmTYseTgb-ILxWRQLdIr61?usp=sharing). Files are listed below: 
1. All_S1_S6fraction_1hBrdu_0Chase_nascent_merged_sort.cov (From Julien, processed Time 0 reads of replicate 1, which is unavailable on GEO, at least  need merge, filter, and process)
2. GSM2182523_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep1_processed_sort.cov (From Julien, Time 1h, replicate 1, also available on GEO, same data with different format, former .cov format, later .bed format)
3. GSM2182524_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep2_processed_sort.cov (Time 1h, rep2, similar with 2)
4. GSM2182525_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep1_processed_sort.cov (Time 4h, rep1, similar with 2)
5. GSM2182526_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep2_processed_sort.cov (Time 4h, rep2, similar with 2)
6. GSM2182527_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep1_processed_sort.cov (Time 16h, rep1, similar with 2)
7. GSM2182528_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep2_processed_sort.cov (Time 16h, rep2, similar with 2)
8. GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3.bed (Time 0h, replicate 3, DOWNLOAD FROM GEO in .bed format)
9. GSM2642846_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep3.bed (Time 1h, replicate 3, from GEO)
10. GSM2642848_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep3.bed (Time 4h, replicate 3, from GEO)
11. GSM2642850_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep3.bed (Time 16h, replicate 3, from GEO)

### 1.2 Data Unused
1. GSM2642847_HUES64_WT_AllCells_1hBrdU_2hChase_nascent_rep3.bed (Time 2h, replicate 3, from GEO data too few not enough data at 2h)
2. GSM2642849_HUES64_WT_AllCells_1hBrdU_8hChase_nascent_rep3.bed (Time 8h, replicate 3, same reason as above)

### 1.3 Processing downloaded data
CODE AVAILABLE at [Google Drive](https://drive.google.com/drive/folders/1ryZGWIVyW0Q8E_WWyfnlR5-bepYlQDJS?usp=sharing) /`step1_download_data_from_geo_and_merge_reads.sh `

1. Create a Project folder,  denoted as BASE_DIR, from now on (in code). 

``` bash
# step1_download_data_from_geo_and_merge_reads.sh

BASE_DIR=~/PycharmProjects/pipelines_of_DNAm_kinetics
DATA_DIR=$BASE_DIR/data
mkdir -p $DATA_DIR
cd $DATA_DIR
```

2. Manually unzip and move all downloaded files into DATA_DIR

``` bash
# !! **** Download THE DATA FROM Google Drive OR GEO ACCORDING TO THE NOTE, MANUALLY UNZIP AND COPY YOUR DATA INTO ~/PycharmProjects/pipelines_of_DNAm_kinetics/data ****!!
```

3. Transform and sort downloaded data into target format  (using bash tool awk)

``` bash
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

```

4. Use `step4_covert_standard_bed_into_reduced_bed.py` to transform standard bed files 

``` bash
python ../python_code/step4_covert_standard_bed_into_reduced_bed.py -input ./GSM2642845_HUES64_WT_AllCells_1hBrdU_0Chase_nascent_rep3.bed -output ./0h_rep2.bed
python ../python_code/step4_covert_standard_bed_into_reduced_bed.py -input ./GSM2642846_HUES64_WT_AllCells_1hBrdU_1hChase_nascent_rep3.bed -output ./1h_rep3.bed
python ../python_code/step4_covert_standard_bed_into_reduced_bed.py -input ./GSM2642848_HUES64_WT_AllCells_1hBrdU_4hChase_nascent_rep3.bed -output ./4h_rep3.bed
python ../python_code/step4_covert_standard_bed_into_reduced_bed.py -input ./GSM2642850_HUES64_WT_AllCells_1hBrdU_16hChase_nascent_rep3.bed -output ./16h_rep3.bed
```


5. Merge the data of replicates  into  one file for each individual time point. e.g. The output for merging `1h_rep1.bed, 1h_rep2.bed,1h_rep3.bed` is `1h_rep0.bed` with the sum reads of the three.

``` bash
declare -a TIMEPOINTS=('0h' '1h' '4h' '16h') # 
declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

for tp in "${TIMEPOINTS[@]}"; 
do
  echo $tp
  echo $tp*.bed
  f1=$tp"_rep0.bed"
  cat $tp*.bed  | gsort -k 1,1 -k2,2n --parallel=8  -S 50% | bedtools merge -c 4,5 -o sum>$f1
done
```

6. Separate each merged bed file by chromosome.
e.g, the output for `0h_rep0.bed` will be: `chr1.bed, chr2.bed chr3.bed … chr22.bed` in the `0h_rep0 ` directory.
``` bash
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
```

7. Construct AllDat Matrix in Malab file format  for each chromosome, e.g `AllDat_chr1.mat`

``` bash
python  ../python_code/step7_construct_AllDat_matrix.py -data_dir . -out_dir ./AllDat
```

8. Rate inference by running `rate_inference_pipeline.m` . This will give the inferred randomT and non-RandomT rates in `../data/Rates/Rates_RandomT` and  `../data/Rates/Rates` directories, and each dir will have `Rates_chr1, Rates_chr2, …, Rates_chr22.mat` files.

``` bash
cd ../matlab_code/
matlab  -nodesktop -nodisplay -nosplash -r "MLEInference('../data/AllDat', '../data/Rates', 1, 22)"
```

9. Convert rates from .mat format into .bed format.

``` bash
cd ../python_code
python step9_convert_rate_from_mat_file_to_bed_file.py -rate_dir ../data/Rates/
```





























