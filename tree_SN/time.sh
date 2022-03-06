#!/bin/bash
 
#PBS -q dssc 
#PBS -l nodes=1:ppn=24
#PBS -l walltime=20:00:00 
 
module load openmpi/4.0.3/gnu/9.3.0 

cd $PBS_O_WORKDIR 
avg_time() {
    #
    # usage: avg_time n command ...
    #
    n=$1; shift
    (($# > 0)) || return                   # bail if no command given
    for ((i = 0; i < n; i++)); do
        { time -p "$@" &>/dev/null; } 2>&1 # ignore the output of the command
                                           # but collect time's output in stdout
    done | awk '
        /real/ { real = real + $2; nr++ }
        /user/ { user = user + $2; nu++ }
        /sys/  { sys  = sys  + $2; ns++}
        END    {
                 if (nr>0) printf("real %f\n", real/nr);
                 if (nu>0) printf("user %f\n", user/nu);
                 if (ns>0) printf("sys %f\n",  sys/ns)
               }'
}

export OMP_PLACE=core
export OMP_PROC_BIND=close

############### SN32 10E8 ###############
#echo > TEST_omp_sn_32_10e8out.txt
#for i in {1..24}
#do
#    echo "thread ${i}" >> TEST_omp_sn_32_10e8out.txt
#    avg_time 15 ./tree_sn_32.x 100000000 0 ${i} >> TEST_omp_sn_32_10e8out.txt
#    echo >> TEST_omp_sn_32_10e8out.txt
#done

############### SN16 10E8 ###############
#echo > TEST_omp_sn_16_10e8out.txt
#for i in {1..24}
#do
#    echo "thread ${i}" >> TEST_omp_sn_16_10e8out.txt
#    avg_time 15 ./tree_sn_16.x 100000000 0 ${i} >> TEST_omp_sn_16_10e8out.txt
#    echo >> TEST_omp_sn_16_10e8out.txt
#done

############### SN32 10E7 ###############
#echo > TEST_omp_sn_32_10e7out.txt
#for i in {1..24}
#do
#    echo "thread ${i}" >> TEST_omp_sn_32_10e7out.txt
#    avg_time 20 ./tree_sn_32.x 10000000 0 ${i} >> TEST_omp_sn_32_10e7out.txt
#    echo >> TEST_omp_sn_32_10e7out.txt
#done

############### SN16 10E7 ###############
#echo > omp_sn_16_10e7out.txt
#for i in {1..24}
#do
#    echo "thread ${i}" >> omp_sn_16_10e7out.txt
#    avg_time 20 ./tree_sn_16.x 10000000 0 ${i} >> omp_sn_16_10e7out.txt
#    echo >> omp_sn_16_10e7out.txt
#done

############### SN32 10E9 ###############
#echo > omp_sn_32_10e9out.txt
for i in {15..24}
do
    echo "thread ${i}" >> omp_sn_32_10e9out.txt
    avg_time 6 ./tree_sn_32.x 1000000000 0 ${i} >> omp_sn_32_10e9out.txt
    echo >> omp_sn_32_10e9out.txt
done

############### SN16 10E9 ###############
#echo > omp_sn_16_10e9out.txt
for i in {1..24}
do
    echo "thread ${i}" >> omp_sn_16_10e9out.txt
    avg_time 6 ./tree_sn_16.x 1000000000 0 ${i} >> omp_sn_16_10e9out.txt
    echo >> omp_sn_16_10e9out.txt
done
