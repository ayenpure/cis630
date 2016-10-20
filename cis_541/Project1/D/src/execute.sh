begin=$(date +%s)
./project1D
end=$(date +%s)
tottime=$(expr $end - $begin)
echo $tottime
