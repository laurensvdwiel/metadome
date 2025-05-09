#!/bin/bash

#SBATCH --job-name=metadomain
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user fferraro@stanford.edu
#SBATCH --nodes=1
#SBATCH --account=smontgom
#SBATCH --partition=batch
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=00:10:00

usage ()
{
  echo 'This scripts builds and annotate protein MetaDomains starting from Gencode, Uniprot, and PFAM input data.'
  echo ''
  echo "Usage : $(basename "$0") -s <pixi.toml> -c <config> -r <clean-up> -t <int> -g <memory>"
  echo '           -s Path to directory containing the manifest, the script, and where the analysis will be run.'
  echo '           -c config.yaml that specifies inputs URLS to download the data from'
  echo '           -r Boolean to specify to remove intermediate files at the end (default=False, alternative True)'
  echo ''  
  exit
}

if [ "$#" -ne 6 ]
then
  usage
fi

while [ "$1" != "" ]; do
case $1 in
        -s )           shift
                       PIXI_MANIFEST_DIR=$1
                       ;;
        -c )           shift
                       RUN_CONFIG=$1
                       ;;
        -r )           shift
                       CLEAN_UP=$1
                       ;;
        * )            usage
                       ;;
    esac
    shift
done

PIXI_MANIFEST=$PIXI_MANIFEST_DIR/pixi.toml


echo -e "Using pixi manifest-path:\n"$PIXI_MANIFEST_DIR
pixi workspace --manifest-path $PIXI_MANIFEST name get
echo ""
echo -e "\nUsing config file paths:\n"$RUN_CONFIG
echo ""
echo "Checking environment..."
echo ""

echo " - Seqkit"
pixi run --manifest-path $PIXI_MANIFEST which seqkit
echo " - MakeBLASTdb"
pixi run --manifest-path $PIXI_MANIFEST which makeblastdb
echo " - BLASTp"
pixi run --manifest-path $PIXI_MANIFEST which blastp
echo " - Python"
pixi run --manifest-path $PIXI_MANIFEST which python

echo "Running pipeline..."
pixi run --manifest-path $PIXI_MANIFEST python $PIXI_MANIFEST_DIR/metadomain_pipeline.py  \
--cores 64  \
--config $RUN_CONFIG \
--working_dir_path $PIXI_MANIFEST_DIR \
--is_for_metadome True 


# Check if DELETE_DATA is set to True and delete data folder if it is
if [ "$CLEAN_UP" = "True" ] || [ "$CLEAN_UP" = "true" ]; then
    echo "Deleting raw and intermediate files..."
    rm -rf data/raw
    rm -rf data/intermediate
    rm -rf data/results/step1_mapping/
    rm -rf data/results/step2_checks/
    rm -rf data/results/step3_mane_annotation/
    rm -rf data/results/step4_pfam/
    rm -rf data/results/step5_single_snv/
    rm -rf data/results/step6_metapositions/
    echo "Removed raw and intermediate files!"
fi