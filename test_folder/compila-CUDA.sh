
file=$1
exe=$2
rm $exe
pgfortran -Bstatic_pgi -Mcuda -Minfo -static-nvidia -tp px -Mcuda=ptxinfo,fastmath,lineinfo,cuda11.0 -Minfo  -ta=tesla:cuda11.0,cc6,cc7.5 $file -o $exe
#pgfortran -fast -fastsse -Bstatic_pgi -Mlarge_arrays -mcmodel=medium -tp px -Mcuda=ptxinfo,fastmath,lineinfo -Minfo  -ta=tesla $file -o $exe

#-Bstatic for static libraries
#-Mcuda works when the code has no .cuf or .CUF extension
