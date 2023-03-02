
file=$1
exe=$2
rm $exe
pgfortran -Bstatic_pgi -Mcuda -Minfo $file -o $exe
#pgfortran -fast -fastsse -Bstatic_pgi -Mlarge_arrays -mcmodel=medium -tp px -Mcuda=ptxinfo,fastmath,lineinfo -Minfo  -ta=tesla $file -o $exe

#-Bstatic for static libraries
#-Mcuda works when the code has no .cuf or .CUF extension
