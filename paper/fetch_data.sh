

mkdir data
mkdir ref
lftp -e "pget -c -n 10 -O ./ /vol1/run/ERR324/ERR3240157/HG00157.final.cram; quit" ftp://ftp.sra.ebi.ac.uk
lftp -e "pget -c -n 10 -O ./ref /vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa; quit" ftp://ftp.1000genomes.ebi.ac.uk
lftp -e "pget -c -n 10 -O ./ref /vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai; quit" ftp://ftp.1000genomes.ebi.ac.uk

# Convert Cram to Bam
samtools view -@ $(nproc) -b -T ./ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -o ./data/HG00157.final.bam ./data/HG00157.final.cram
samtools index ./data/HG00157.final.bam
