#!/bin/bash
#$ -cwd

rm script_SLIM3_SNP_BP_SLIM3_SIMOVERDOM_chr4_2_trozo.sh.*

if [ $# -ne 7 ]
then
	echo "Usage: $0 <slim3INPUTfile> <N> <M1> <M2> <M3> <REPS> <n>"
	exit 1
fi

#command to run with 1 replicate: qsub script_SLIM3_SNP_BP_SLIM3_SIMOVERDOM.sh INPUT_slim3 20 100 1000 1 001
#command to run with 10 replicates: qsub script_SLIM3_SNP_BP_SLIM3_SIMOVERDOM.sh INPUT_slim3 20 100 1000 10 002

#Set arguments
INPUT=$1
N=$2
M1=$3
M2=$4
M3=$5
REPS=$6
n=$7

#Working directory
WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_SLIM_$INPUT" ]
then
rm -r $WDIR/OUTPUT_SLIM_$INPUT
fi

mkdir -p $WDIR/OUTPUT_SLIM_$INPUT
mkdir -p /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################

cp SNP_BP_SLIM3_2 /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
cp SIMOVERDOM /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
cp $INPUT /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
cp recombination_map_chr4_2trozo_escalado.txt /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
cp shell_F_values_4_2 /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/shell_F_values
cp plink /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
cp gcta-1.94.1 /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID

################################ REPLICATES #############################

for ((r=1; r<=$REPS; r++))
do

###################### run SLIM3 #########################

module load SLiM/3.3.2

START=$(date +%s)
slim $INPUT > slimout
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Slim3 took 		$DIFF seconds" >> timefile

cp slimout slimout$r
cp collect_par_mutations.txt collect_par_mutations$r.txt

#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/slimout$r $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/collect_par_mutations$r.txt $WDIR/OUTPUT_SLIM_$INPUT

###################### SNP_BP_SLIM3_2 #########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)

./SNP_BP_SLIM3_2<<@
-99
$N	N
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_BP took 	$DIFF seconds" >> timefile

cp dataBP.ped dataBP$r.ped
cp dataBP.map dataBP$r.map
cp list_allsnps list_allsnps$r
cp list_qtls list_qtls$r

#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dataBP$r.ped $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dataBP$r.map $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/list_allsnps$r $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/list_qtls$r $WDIR/OUTPUT_SLIM_$INPUT

###################### SIMOVERDOM M1 ########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)
./SIMOVERDOM<<@
0
-99
$N    NINDNP (max 10000)
$M1     NIND (max 10000)
100     classes
1	   Replicates
@

DIFF=$(( $END - $START ))
echo "SIMOVERDOM took 		$DIFF seconds" >> timefile

cat genfile.dat >> genfile.dat.$M1
cp dfilename.dat dfilename$r.$M1.dat
cat data.F >> data.F.$M1
cat len_ROH >> len_ROH.$M1
cp summaryID summaryID$r.$M1
cp summaryM summaryM$r.$M1
cp summaryV summaryV$r.$M1
cp summaryR summaryR$r.$M1
cp summarySqE summarySqE$r.$M1
cp summaryB summaryB$r.$M1

cat summaryID$r.$M1 >> SUMMARYID.$M1
cat summaryM$r.$M1 >> SUMMARYM.$M1
cat summaryV$r.$M1 >> SUMMARYV.$M1
cat summaryR$r.$M1 >> SUMMARYR.$M1
cat summarySqE$r.$M1 >> SUMMARYSqE.$M1
cat summaryB$r.$M1 >> SUMMARYB.$M1

cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/data.F.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/len_ROH.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/genfile.dat.$M1 $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dfilename$r.$M1.dat $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYID.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYM.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYV.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYR.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYSqE.$M1 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYB.$M1 $WDIR/OUTPUT_SLIM_$INPUT


#eliminar los documentos $r una vez se han pasado
rm dfilename$r.$M1.dat
rm summaryID$r.$M1
rm summaryM$r.$M1
rm summaryV$r.$M1
rm summaryR$r.$M1
rm summarySqE$r.$M1
rm summaryB$r.$M1




###################### SIMOVERDOM M2 ########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)
./SIMOVERDOM<<@
0
-99
$N    NINDNP (max 10000)
$M2     NIND (max 10000)
100     classes
1       Replicates
@

DIFF=$(( $END - $START ))
echo "SIMOVERDOM took 		$DIFF seconds" >> timefile

cat genfile.dat >> genfile.dat.$M2
cp dfilename.dat dfilename$r.$M2.dat
cat data.F >> data.F.$M2
cat len_ROH >> len_ROH.$M2
cp summaryID summaryID$r.$M2
cp summaryM summaryM$r.$M2
cp summaryV summaryV$r.$M2
cp summaryR summaryR$r.$M2
cp summarySqE summarySqE$r.$M2
cp summaryB summaryB$r.$M2

cat summaryID$r.$M2 >> SUMMARYID.$M2
cat summaryM$r.$M2 >> SUMMARYM.$M2
cat summaryV$r.$M2 >> SUMMARYV.$M2
cat summaryR$r.$M2 >> SUMMARYR.$M2
cat summarySqE$r.$M2 >> SUMMARYSqE.$M2
cat summaryB$r.$M2 >> SUMMARYB.$M2

cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/data.F.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/len_ROH.$M2 $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/data.map $WDIR/OUTPUT_SLIM_$INPUT/chromosome$r.map
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/data.ped $WDIR/OUTPUT_SLIM_$INPUT/chromosome$r.ped
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/genfile.dat.$M2 $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dfilename$r.$M2.dat $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYID.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYM.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYV.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYR.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYSqE.$M2 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYB.$M2 $WDIR/OUTPUT_SLIM_$INPUT


#eliminar los documentos $r una vez se han pasado
rm dfilename$r.$M2.dat
rm summaryID$r.$M2
rm summaryM$r.$M2
rm summaryV$r.$M2
rm summaryR$r.$M2
rm summarySqE$r.$M2
rm summaryB$r.$M2




###################### SIMOVERDOM M3 ########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)
./SIMOVERDOM<<@
0
-99
$N    NINDNP (max 10000)
$M3     NIND (max 10000)
100     classes
1       Replicates
@

DIFF=$(( $END - $START ))
echo "SIMOVERDOM took 		$DIFF seconds" >> timefile

cat genfile.dat >> genfile.dat.$M3
cp dfilename.dat dfilename$r.$M3.dat
cat data.F >> data.F.$M3
cat len_ROH >> len_ROH.$M3
cp summaryID summaryID$r.$M3
cp summaryM summaryM$r.$M3
cp summaryV summaryV$r.$M3
cp summaryR summaryR$r.$M3
cp summarySqE summarySqE$r.$M3
cp summaryB summaryB$r.$M3

cat summaryID$r.$M3 >> SUMMARYID.$M3
cat summaryM$r.$M3 >> SUMMARYM.$M3
cat summaryV$r.$M3 >> SUMMARYV.$M3
cat summaryR$r.$M3 >> SUMMARYR.$M3
cat summarySqE$r.$M3 >> SUMMARYSqE.$M3
cat summaryB$r.$M3 >> SUMMARYB.$M3

cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/data.F.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/len_ROH.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/genfile.dat.$M3 $WDIR/OUTPUT_SLIM_$INPUT
#cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dfilename$r.$M3.dat $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYID.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYM.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYV.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYR.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYSqE.$M3 $WDIR/OUTPUT_SLIM_$INPUT
cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/SUMMARYB.$M3 $WDIR/OUTPUT_SLIM_$INPUT

##cp -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/* $WDIR/OUTPUT_SLIM_$INPUT

#eliminar los documentos $r una vez se han pasado
rm dfilename$r.$M3.dat
rm summaryID$r.$M3
rm summaryM$r.$M3
rm summaryV$r.$M3
rm summaryR$r.$M3
rm summarySqE$r.$M3
rm summaryB$r.$M3

##############################################

#BORRAR LO DEL SLIM
rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/slimout$r 
rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/collect_par_mutations$r.txt 

rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dataBP$r.ped 
rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/dataBP$r.map 
rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/list_allsnps$r 
rm /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/list_qtls$r 

done

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/pilarSLIM$INPUT/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
