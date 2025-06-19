nbits=$1
NUMBETAS=$2

target1=`echo $1 $2 | awk '{printf("create_nbitsLUT%02d_NB%02d",$1,$2)}'`
target2=`echo $1 $2 | awk '{printf("generate_bits_nbitsLUT%02d_NB%02d",$1,$2)}'`

echo $target1
echo $target2

gcc preparation.c -o $target1 -DNUMBITSPREBUSQUEDAS=$nbits -DNUMBETAS=$NUMBETAS -lm -lquadmath -Wall -Wshadow
nvcc random.cu bits-random-ini.cu bits-random-io.cu  bits-random.cu -o $target2 -DNUMBITSPREBUSQUEDAS=$nbits -DNUMBETAS=$NUMBETAS -O3 -lm 
