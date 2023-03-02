
rm $2
pgfortran -Bstatic -Mcuda $1 -o $2

#-Bstatic for static libraries
#-Mcuda works when the code has no .cuf or .CUF extension
