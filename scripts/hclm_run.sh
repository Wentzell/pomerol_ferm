#!/bin/bash

# job name (default is name of pbs script file)
#PBS -N pomerol_ferm
#
# resource limits: max. wall clock time during which job can be running
#PBS -l walltime=720:00:00
#
#number of nodes and cores per node
#PBS -l nodes=4:ppn=12
#
# path/filename for standard output
#PBS -o $PBS_O_WORKDIR/log/run.log
#
# path/filename for standard error
#PBS -e $PBS_O_WORKDIR/log/run.err
#
# send me mail when job begins/ends/aborts (b/e/a)
#  #PBS -m bea
#
# export all my environment variables to the job
#PBS -V
#
#
# Using PBS - Environment Variables
# When a batch job starts execution, a number of environment variables are
# predefined, which include:
#
#      Variables defined on the execution host.
#      Variables exported from the submission host with
#                -v (selected variables) and -V (all variables).
#      Variables defined by PBS.
#
# The following reflect the environment where the user ran qsub:
# PBS_O_HOST    The host where you ran the qsub command.
# PBS_O_LOGNAME Your user ID where you ran qsub.
# PBS_O_HOME    Your home directory where you ran qsub.
# PBS_O_WORKDIR The working directory where you ran qsub.
#
# These reflect the environment where the job is executing:
# PBS_ENVIRONMENT       Set to PBS_BATCH to indicate the job is a batch job, or
#                 to PBS_INTERACTIVE to indicate the job is a PBS interactive job.
# PBS_O_QUEUE   The original queue you submitted to.
# PBS_QUEUE     The queue the job is executing from.
# PBS_JOBID     The job's PBS identifier.
# PBS_JOBNAME   The job's name.
###

# -----------------------------------------------------

echo "Starting script..."

cd $PBS_O_WORKDIR

NPROCS=`wc -l < $PBS_NODEFILE`

#trap '' USR1 USR2

/opt/openmpi/bin/mpirun -x PATH -x LD_LIBRARY_PATH -machinefile $PBS_NODEFILE -np $NPROCS bin/anderson -U 1.0 -b 20.0 -l -0.7 -t 0.45 -l -0.15 -t 0.34 -l 0.15 -t 0.34 -l 0.7 -t 0.45 --calc2pgf --calcgf --wf 64 --coefftol 1e-8 --reducetol 1e-4 > pom-$PBS_JOBID.log 2> pom-err-$PBS_JOBID.log

echo "Script complete."
