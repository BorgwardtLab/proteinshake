TAG="28DEC2022"
CPUS=500
MEM=$((100*1024/$CPUS)) # 100GB
LOGDIR="./release/$TAG/logs"

mkdir -p $LOGDIR

# labeled datasets
for DATASET in GeneOntologyDataset EnzymeCommissionDataset PfamDataset ProteinProteinInterfaceDataset ProteinLigandInterfaceDataset TMAlignDataset SCOPDataset; do
    sbatch --ntasks 1 --cpus-per-task $CPUS --mem-per-cpu $MEM --time 29-1 -o "$LOGDIR/$DATASET.log" -J $DATASET --wrap "python -m scripts.release --njobs $CPUS --tag $TAG --dataset $DATASET"
done

# RCSB
sbatch --ntasks 1 --cpus-per-task 20 --mem-per-cpu 2500 --time 23:00:00 -o "$LOGDIR/RCSBDataset.log" -J "RCSB" --wrap "python -m scripts.release --njobs 300 --tag $TAG --dataset RCSBDataset"

# AlphaFold
for ORGANISM in arabidopsis_thaliana caenorhabditis_elegans candida_albicans danio_rerio dictyostelium_discoideum drosophila_melanogaster escherichia_coli glycine_max homo_sapiens methanocaldococcus_jannaschii mus_musculus oryza_sativa rattus_norvegicus saccharomyces_cerevisiae schizosaccharomyces_pombe zea_mays swissprot; do
    sbatch --ntasks 1 --cpus-per-task 20 --mem-per-cpu 2500 --time 23:00:00 -o "$LOGDIR/AlphaFoldDataset_$ORGANISM.log" -J $ORGANISM --wrap "python -m scripts.release --njobs 300 --tag $TAG --dataset AlphaFoldDataset --organism $ORGANISM"
done
