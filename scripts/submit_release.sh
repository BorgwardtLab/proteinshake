TAG="27DEC2022"
CPUS=20
MEM=$((50*1024/$CPUS))

# labeled datasets
for DATASET in GeneOntologyDataset EnzymeCommissionDataset PfamDataset ProteinProteinInterfaceDataset ProteinLigandInterfaceDataset TMAlignDataset SCOPDataset; do
    sbatch --ntasks 1 --cpus-per-task $CPUS --mem-per-cpu $MEM --time 10-1 -o "./release/$TAG/logs/$DATASET.log" --wrap "python -m scripts.release --njobs $CPUS --tag $TAG --dataset $DATASET"
done

# RCSB
sbatch --ntasks 1 --cpus-per-task 20 --mem-per-cpu 2500 --time 4:00 -o "./release/$TAG/logs/RCSBDataset.log" --wrap "python -m scripts.release --njobs 300 --tag $TAG --dataset RCSBDataset"

# AlphaFold
for ORGANISM in arabidopsis_thaliana caenorhabditis_elegans candida_albicans danio_rerio dictyostelium_discoideum drosophila_melanogaster escherichia_coli glycine_max homo_sapiens methanocaldococcus_jannaschii mus_musculus oryza_sativa rattus_norvegicus saccharomyces_cerevisiae schizosaccharomyces_pombe zea_mays swissprot; do
    sbatch --ntasks 1 --cpus-per-task 20 --mem-per-cpu 2500 --time 4:00 -o "./release/$TAG/logs/AlphaFoldDataset_$ORGANISM.log" --wrap "python -m scripts.release --njobs 300 --tag $TAG --dataset AlphaFoldDataset --organism $ORGANISM"
done
