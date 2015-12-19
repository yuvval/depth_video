To compile mex type:
mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99"  slicmex.c
mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99"  slicomex.c
