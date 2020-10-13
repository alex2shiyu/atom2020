#!/bin/bash

dst_dir=module
cd ..
if [ -d $dst_dir ];then
    rm -rf $dst_dir
fi
mkdir $dst_dir
cd $dst_dir
for i in `ls ../lib`;
do
    destfile="`basename \"$i\" .f90`"
    aimfile="$destfile"".cpython-36m-x86_64-linux-gnu.so"
    echo "-----------------------------------------------------"
    echo "--   make $destfile   "
    echo "--   aim  $aimfile    "
    echo "-------------------------"
    if [ $destfile != "naivemakefile.sh" -a ! -e $aimfile ]; then
#      f2py --no-lower -m "$destfile"  "$destfile".pyf "$destfile".f90
#      f2py -c "$destfile".pyf "$destfile".f90
        f2py  -m "$destfile" -c '../lib/'"$destfile".f90 
    elif [ $destfile = "naivemakefile.sh" ];then
        echo "omit naivemakefile.sh"
        echo "-----------------------------------------------------"
        echo ""
        echo ""
        echo ""
    else
        echo "$aimfile exists"
        echo "-----------------------------------------------------"
        echo ""
        echo ""
    fi
done

