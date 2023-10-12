for k in $(echo 5);do
   rm seg_list.dat
   rm list
   for f in $(ls segment_"$k"_input_af_*sys.pdb);do
      grep ATOM $f |tail -n 1|awk '{print $2}' >>list
   done
   atom=$(sort -n -k 1 list|uniq -c |sort -r -n -k 1 |awk '{print $2}' |head -n 1)
   for i in $(ls segment_"$k"_input_af_*sys.pdb);do
      atoms=$(grep ATOM $i |tail -n 1|awk '{print $2}')
      echo $atoms
      if [ $atoms -eq $atom ]; then
         echo "$i " >> seg_list.dat
      fi
   done

   xtcs=$(cat seg_list.dat |sed 's/_sys.pdb/.rebuilt.xtc/'g|xargs -n 100000000 )
   echo $xtcs
   gmx_mpi trjcat -f $xtcs  -cat -o segment_"$k"_input_af_rebuilt.xtc

done
