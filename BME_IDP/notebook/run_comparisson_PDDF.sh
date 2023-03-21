#!/bin/bash 

for name in $(echo "DSS1"); do
       echo $name
       for j in $(echo 0 1 2 3 4);do

		 python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 71  --SAXS_bw 1
		mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
		mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
        	mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/

        done
done

for name in $(echo "Sic1"); do
        echo $name
        for j in $(echo 0 1 2 3 4);do
        	python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 90  --SAXS_bw 1
	        mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/

        done
done

for name in $(echo "ProTa"); do
        echo $name
        for j in $(echo 0 1 2 3 4);do
        	python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 111  --SAXS_bw 1
                mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	done
done

for name in $(echo "NHE6cmdd"); do
        for j in $(echo 0 1 2 3 4);do
        echo $name
	       python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 116  --SAXS_bw 1
                mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
        done
done


for name in $(echo "ANAC046"); do
        echo $name
        for j in $(echo 0 1 2 3 4);do
	        python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 167  --SAXS_bw 10
	        mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
        done
done

for name in $(echo "GHR-ICD"); do
        echo $name
	 for j in $(echo 0 1 2 3 4);do
        	python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 351  --SAXS_bw 10
	        mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
        done
done

for name in $(echo "aSyn"); do
        echo $name
        for j in $(echo 0 1 2 3 4);do
               python3 PDDF_comparisson.py --master_folder ../../SAXS/PDDF_comparisson/"$name"_distmat --pdb ../../Pesce_dataset_backmaped/"$name"/0_em.pdb --xtc ../../Pesce_dataset_backmaped/"$name"/trjout_bm_em.xtc --sys_n $name --SAXS_PPDF_file SAXS/"$name"_Pr_data_only.out --ev_CB 3 --nresid 140  --SAXS_bw 10
                mkdir ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*txt ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
	        mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*log ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*dat ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
                mv ../../SAXS/PDDF_comparisson/"$name"_distmat/*pdf ../../SAXS/PDDF_comparisson/"$name"_distmat/iter_$j/
        done

done
