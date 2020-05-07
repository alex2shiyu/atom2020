#!/bin/bash

for i in `ls ../lib`;
do
    destfile="`basename \"$i\" .f90`"
    echo "make $destfile"
    if [ $destfile != "naivemakefile.sh" ]; then
#      f2py --no-lower -m "$destfile"  "$destfile".pyf "$destfile".f90
#      f2py -c "$destfile".pyf "$destfile".f90
        f2py  -m "$destfile" -c '../lib/'"$destfile".f90 
    else
        echo "omit naivemakefile.sh"
    fi
done

