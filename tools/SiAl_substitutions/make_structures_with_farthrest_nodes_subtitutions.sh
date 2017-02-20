#! /bin/bash
make
if [ -d results ] ; then rm -R results ; fi
if [ -f out.txt ] ; then rm out.txt ; fi
mkdir results
echo "#Config:  Dist_ave   Dist_min_ave  Cost" > out.txt
for i in $(seq 1 50) ; do
    ./farthrest_nodes_subtitutions > out_farthrest_nodes_subtitutions_$i
    mv out.gin out_${i}.gin
    sleep 0.5
    pot=$(tail   -n1 out_farthrest_nodes_subtitutions_${i} | awk '{print $4}')
    dist1=$(tail -n1 out_farthrest_nodes_subtitutions_${i} | awk '{print $2}')
    dist2=$(tail -n1 out_farthrest_nodes_subtitutions_${i} | awk '{print $3}')
    echo $i":" $dist1 $dist2 $pot >> out.txt
    mv out_farthrest_nodes_subtitutions_$i out_${i}.gin results
    echo "config" $i":" ${dist1} ${dist2} ${pot}
done
echo 'salida en out.txt'
rm farthrest_nodes_subtitutions
sort -nk4 out.txt > tmp
mv tmp out.txt
exit 0
