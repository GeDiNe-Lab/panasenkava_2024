# Sample of interests are listed here
folders="CS12_Telencephalon CS12_Diencephalon CS13_prosencephalon CS13_anterior_midbrain CS13_center_midbrain CS13_hindbrain CS14_3_diencephalon CS15_forebrain2 CS15_forebrain"

cd "/home/jules/Documents/phd/Data/literature/Eze_et_al/download/"
for folder in $folders
do
    wget https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/counts/"$folder"/"$folder".mex.tar.gz
    tar -xvzf "$folder".mex.tar.gz
done

# removing archives
rm *.mex.tar.gz