begin=$(date +%s%N)
./project1D
end=$(date +%s%N)
tottime=$(expr $end - $begin)
echo $(expr $tottime / 1000000) ms
