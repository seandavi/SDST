#!python3
import subprocess
from multiprocessing import Pool,cpu_count

chroms = ['chr'+str(x+1) for x in list(range(22))]+['chrX','chrY']

def work(cmd):
    print(cmd)
    return subprocess.call(cmd,shell=True)

if __name__=="__main__":
    count = cpu_count()
    pool  = Pool(processes=count)
    cmds = ['module load breakdancer; breakdancer-max -r 4 -h -g /scratch/bdancer.{0}.bed -o {0} -r 4 /scratch/bdancer.cfg > /scratch/bdancer.{0}.txt'.format(chrom) for chrom in chroms]
    cmds += ['module load breakdancer; breakdancer-max -g /scratch/bdancer.trans.bed -h -q 10 -r 4 -t -r 4 /scratch/bdancer.cfg > /scratch/bdancer.trans.txt']
    pool.map(work,cmds)
    subprocess.call('cat /scratch/bdancer.*.txt > /scratch/bdancer.txt',shell=True)
    subprocess.call('cat /scratch/bdancer.*.bed > /scratch/bdancer.bed',shell=True)

















