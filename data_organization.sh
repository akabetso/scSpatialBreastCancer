#Directory
data_dir="data"

#untar the files downloaded from Gene Expression Omnnibus
for file in *.gz do; tar -xvzf $file; done

# Now the the files have been untared but are contain in 26 seperate folders
# Two acions are to be taken
# First renamed the files in each folder to ease manipulations in seurat and gzip the files as well 
# So I apply the loop below to find the barcodes key within the file name and rename (move) to barcodes.tsv only. Same thin apply for features and matrix. 
# gzip the files after renaiming them
cd $data_dir
for dir in CID*; do
    cd $dir
    barcodes=$(find -name "*barcodes.tsv*")
    features=$(find -name "*genes.tsv*")
    matrix=$(find -name "*sparse.mtx*")
    mv $barcodes barcodes.tsv
    gzip barcodes.tsv
    mv $features features.tsv
    gzip features.tsv
    mv $matrix matrix.mtx
    gzip matrix.mtx
    ls
    echo "to next file"
    cd ..
done

cd ..