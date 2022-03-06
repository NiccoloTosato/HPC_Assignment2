#!/bin/bash
 
#PBS -q dssc 
#PBS -l nodes=1:ppn=24
#PBS -l walltime=10:10:00 
 
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


echo > weak_omp_10e7out.txt
for i in 1 2 4 6 8 10 12 14 16 18 20 22 24
do
    echo "thread ${i}" >> weak_omp_10e7out.txt
    avg_time 15 ./tree_omp.x $((10000000*i)) 0 ${i} >> weak_omp_10e7out.txt
    echo >> weak_omp_10e7out.txt
done

echo > weak_omp_10e8out.txt
for i in 1 2 4 6 8 10 12 14 16 18 20 22 24
do
    echo "thread ${i}" >> weak_omp_10e8out.txt
    avg_time 8 ./tree_omp.x $((100000000*i)) 0 ${i} >> weak_omp_10e8out.txt
    echo >> weak_omp_10e8out.txt
done

#echo > weak_omp_10e9out.txt
#for i in {20..24}
#do
#    echo "thread ${i}" >> weak_omp_10e9out.txt
#    avg_time 5 ./tree_omp.x $((10000000*i))00 0 ${i} >> weak_omp_10e9out.txt
#    echo >> weak_omp_10e9out.txt
#done


