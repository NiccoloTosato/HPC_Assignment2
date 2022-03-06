#!/bin/bash
 
#PBS -q dssc 
#PBS -l nodes=1:ppn=24
#PBS -l walltime=10:00:00 
 
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

###############  2mpi 10E8 ###############
echo > strong_hybrid_2_10e8out.txt
for i in {1..12}
do
    echo "thread ${i}" >> strong_hybrid_2_10e8out.txt
    avg_time 8 mpirun -np 2 --mca btl ^openib  --map-by socket:PE=${i} ./tree_hybrid.x $((100000000)) 0 ${i} >> strong_hybrid_2_10e8out.txt
    echo >> strong_hybrid_2_10e8out.txt
done

###############  4mpi 10E8 ###############
echo > strong_hybrid_4_10e8out.txt
for i in {1..6}
do
    echo "thread ${i}" >> strong_hybrid_4_10e8out.txt
    avg_time 8 mpirun -np 4 --mca btl ^openib --map-by socket:PE=${i} ./tree_hybrid.x $((100000000)) 0 ${i} >> strong_hybrid_4_10e8out.txt
    echo >> strong_hybrid_4_10e8out.txt
done

###############  8mpi 10E8 ###############
echo > strong_hybrid_8_10e8out.txt
for i in {1..3}
do
    echo "thread ${i}" >> strong_hybrid_8_10e8out.txt
    avg_time 8 mpirun -np 8 --mca btl ^openib --map-by socket:PE=${i} ./tree_hybrid.x $((100000000)) 0 ${i} >> strong_hybrid_8_10e8out.txt
    echo >> strong_hybrid_8_10e8out.txt
done



###############  2mpi 10E7 ###############
echo > strong_hybrid_2_10e7out.txt
for i in {1..12}
do
    echo "thread ${i}" >> strong_hybrid_2_10e7out.txt
    avg_time 15 mpirun -np 2 --mca btl ^openib  --map-by socket:PE=${i} ./tree_hybrid.x $((10000000)) 0 ${i} >> strong_hybrid_2_10e7out.txt
    echo >> strong_hybrid_2_10e7out.txt
done

###############  4mpi 10E7 ###############
echo > strong_hybrid_4_10e7out.txt
for i in {1..6}
do
    echo "thread ${i}" >> strong_hybrid_4_10e7out.txt
    avg_time 15 mpirun -np 4 --mca btl ^openib --map-by socket:PE=${i} ./tree_hybrid.x $((10000000)) 0 ${i} >> strong_hybrid_4_10e7out.txt
    echo >> strong_hybrid_4_10e7out.txt
done

###############  8mpi 10E7 ###############
echo > strong_hybrid_8_10e7out.txt
for i in {1..3}
do
    echo "thread ${i}" >> strong_hybrid_8_10e7out.txt
    avg_time 15 mpirun -np 8 --mca btl ^openib --map-by socket:PE=${i} ./tree_hybrid.x $((10000000)) 0 ${i} >> strong_hybrid_8_10e7out.txt
    echo >> strong_hybrid_8_10e7out.txt
done
