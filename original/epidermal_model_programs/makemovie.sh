#!/bin/bash
outdir=$1
movnam=$2
ni=$3
nf=$4

sed -e "s/INPUT_DIR xxx/INPUT_DIR .\/${outdir}/g" \
    -e "s/OUTPUT xxx/OUTPUT $movnam/g" \
    -e "s/xxx-xxx/$ni-$nf/g" \
    mpeg_para_all > mpeg_para_all_tmp
echo "encode mpg..."
mpeg_encode mpeg_para_all_tmp > mpeg_encode.log
rm -f mpeg_para_all_tmp
