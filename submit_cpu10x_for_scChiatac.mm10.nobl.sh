#!/usr/bin/env bash
set -eo pipefail

# One DNA FASTQ prefix and one RNA FASTQ prefix per row.
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash $0 FASTQ_PREFIX_TABLE [CONFIG_FILE]" >&2
    exit 1
fi
fqlist=$(realpath "$1")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
config_file=${2:-${CHAIR_CONFIG:-}}
if [[ -n "$config_file" ]]; then
    config_file=$(realpath "$config_file")
    source "$config_file"
fi

: "${ARCref:?Set ARCref in the configuration}"
: "${GEXref:?Set GEXref in the configuration}"
: "${ATACref:?Set ATACref in the configuration}"
: "${fasta:?Set fasta in the configuration}"
: "${blacklist:?Set blacklist in the configuration}"
: "${scom:?Set scom to the CPU/UCSC/MACS2 command prefix}"
: "${dedup:?Set dedup to the GATK MarkDuplicates command}"
: "${juicertool:?Set juicertool to the Juicer tools command}"

export genome=${genome:-mm10}
export gmac=${gmac:-mm}
export shiftsize=${shiftsize:--75}
export macsq=${macsq:-0.001}
export peakext=${peakext:-150}
export peak_signal_min=${peak_signal_min:-1.5}
export addLinkerSeq=${addLinkerSeq-'-A ACGTGATATTTCACGACTCT'}
export minMapScore=${minMapScore:-30}
export mapq=${mapq:-20}
export selfbp=${selfbp:-1}
export bulkSelfbp=${bulkSelfbp:-8000}
export minPETsignif=${minPETsignif:-3}
export chemistry=${chemistry:-ARC-v1}

export NTHREAD=${NTHREAD:-48}
export NTHREAD10x=${NTHREAD10x:-$NTHREAD}
export GB=${GB:-200}
export highGB=${highGB:-320}
export cbdictGB=${cbdictGB:-400}
export cbNthread=${cbNthread:-24}
export hrs=${hrs:-21}
export QUEUE=${QUEUE:-cpu}
export smallQUEUE=${smallQUEUE:-$QUEUE}
export highpart=${highpart:-$QUEUE}
export LINKER_GB=${LINKER_GB:-15}
export MAPPING_GB=${MAPPING_GB:-40}
export MERGE_GB=${MERGE_GB:-36}
export TABLE_GB=${TABLE_GB:-8}
export SUMMARY_GB=${SUMMARY_GB:-2}
export DOWNSTREAM_GB=${DOWNSTREAM_GB:-$highGB}

export pairlabel=singlelinker.paired
export singlabel=singlelinker.single
export nonelabel=none
export dumpfolder=${dumpfolder:-dump_intermediate}
export tmppath="$SCRIPT_DIR/codes"
export python=${python:-python3}
export bedcmd=${bedcmd:-bedtools}
export samcmd=${samcmd:-samtools}
export Rcmd=${Rcmd:-Rscript}
export bgzip=${bgzip:-bgzip}
export tabix=${tabix:-tabix}
export pigz=${pigz:-'pigz -f'}
export bamCoverage=${bamCoverage:-bamCoverage}
export arc10x=${arc10x:-'cellranger-arc count'}
export gex10x=${gex10x:-'cellranger count'}
export atac10x=${atac10x:-'cellranger-atac reanalyze'}
export cpuprog=${cpuprog:-"$scom cpu"}

export statpy_noTag="$tmppath/statBarcode_givenCBlist.py"
export addCBpy="$tmppath/addCBZB_givenCBlist.py"
export injectCBreadID="$tmppath/useCBtag_injectReadID.pl"
export pickUpy="$tmppath/pickingUxxU.xyh.py"
export statpy_bamTag="$tmppath/statCB_PET.py"
export buildCBpy="$tmppath/generateCBlist_for_cpuRead.v2.py"
export mergeTable="$tmppath/mergeCBstats.1genome.arg.R"
export bamtofrag="$tmppath/simply_bam2frags_singleEnd.xyh.py"
export sortfrag="$tmppath/sort_frags.R"
export sumsh="$tmppath/extract_summary_cpu10x_1genome.sh"
export filtersmallsize="$tmppath/filterinsertsizeless1000.xyh.py"
export Sbam2fragment="$tmppath/singletonBam2fragment.xyh.py"
export CountBarcodePET="$tmppath/CountBarcodePET.py"

export CHAIR_MIN_RNA=${CHAIR_MIN_RNA:-100}
export CHAIR_MIN_ATAC=${CHAIR_MIN_ATAC:-100}
export CHAIR_MIN_PET=${CHAIR_MIN_PET:-50}
export CHAIR_MAX_MT=${CHAIR_MAX_MT:-20}
export CHAIR_MAX_RIBO=${CHAIR_MAX_RIBO:-50}
export CHAIR_DOUBLET_RATE=${CHAIR_DOUBLET_RATE:-0.054}
export CHAIR_SEED=${CHAIR_SEED:-1234}
export CHAIR_RNA_DIMS=${CHAIR_RNA_DIMS:-30}
export CHAIR_ATAC_DIMS=${CHAIR_ATAC_DIMS:-30}
export CHAIR_RNA_RESOLUTION=${CHAIR_RNA_RESOLUTION:-0.8}
export CHAIR_ATAC_RESOLUTION=${CHAIR_ATAC_RESOLUTION:-0.8}
export CHAIR_WNN_RESOLUTION=${CHAIR_WNN_RESOLUTION:-0.5}

mkdir -p "${OUTPUT_DIR:-$PWD}"
maindir=$(realpath "${OUTPUT_DIR:-$PWD}")
die() { echo "$*" >&2; exit 1; }
function odir {
echo "odir cmd -- $1 "
outputDir=$1
if [ ! -d "$outputDir" ]; then
    echo "Create $outputDir"
    mkdir -p "$outputDir"
fi
}

function submitjob {
    local slurmfile=$1
    shift
    local dependency="" id
    for id in $*; do
        case ":$dependency:" in
            *":$id:"*) ;;
            *) dependency="${dependency:+$dependency:}$id" ;;
        esac
    done
    if [[ -n "$dependency" ]]; then
        sbatch --parsable --dependency="afterok:$dependency" "$slurmfile"
    else
        sbatch --parsable "$slurmfile"
    fi | cut -d';' -f1
}

function submitjob2ok { submitjob "$@"; }
function submitjobOptional { submitjob "$@"; }
function submitjobOptional8 { submitjob "$@"; }

export dedupIDs1='' #for dependent job that need multiple previous jobIDs
export dedupIDs2='' #for dependent job that need multiple previous jobIDs
export jidMerge='' #for dependent job that need multiple previous jobIDs
function storeDedupJOBid1 {
echo "storeDedupJOBid1 cmd -- $1 "
    export dedupIDs1="$dedupIDs1 $1"
}
function storeDedupJOBid2 {
echo "storeDedupJOBid2 cmd -- $1 "
    export dedupIDs2="$dedupIDs2 $1"
}
function storeMergeJOBid {
echo "storeMergeJOBid cmd -- $1 "
    export jidMerge="$jidMerge $1"
}


export storedIDtsv='' #for dependent job that need multiple previous jobIDs
function storeJOBid {
echo "storeJOBid cmd -- $1 "
    export storedIDtsv="$storedIDtsv $1"
}


#--------------------
# START job functions
#------------------------------------------------------------------------------------------------------------

function cpuLinker {
echo "cpuLinker cmd -- $1 "
    fqpref=$1
    export RUN=$(basename $fqpref)
    nf=$(/bin/ls ${fqpref}_*L*_R1*.fastq.gz |wc -l)
    r1=$(/bin/ls ${fqpref}_*L*_R1*.fastq.gz)
    r3=$(/bin/ls ${fqpref}_*L*_R3*.fastq.gz)

export fastq1=${RUN}_R1.fastq.gz
export fastq2=${RUN}_R3.fastq.gz

export cpudir=$maindir/$RUN/bulkChiatac
jobname=${RUN}.linker
pbs=${jobname}.slurm


cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD  # number of threads
#SBATCH --mem ${LINKER_GB}G # memory pool for all cores # XYH some Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
#cd \$SLURM_SUBMIT_DIR

cd "$cpudir"
date


#Taking care of multiple lanes fastq files
if [ "$nf" -gt  1 ]; then
   echo "Merging $nf fastq files to ${fastq1}:"
   /bin/ls ${fqpref}_*L*_R1*.fastq.gz 
   cat ${fqpref}_*L*_R1*.fastq.gz > $fastq1
   echo "Merging $nf fastq files to ${fastq2}:"
   /bin/ls ${fqpref}_*L*_R3*.fastq.gz
   cat ${fqpref}_*L*_R3*.fastq.gz  > $fastq2
else
    echo "Soft-linking $r1 file to $fastq1"
    ln -s -f $r1 $fastq1
    echo "Soft-linking $r3 file to $fastq2"
    ln -s -f $r3 $fastq2
fi

# perform linker detection and generation of different category of fastq files
echo "Linker detection on: ${RUN} " 

$cpuprog stag -W -T 18 -t ${NTHREAD} $addLinkerSeq -O ${RUN} $fastq1 $fastq2 
echo "--- linker detection completed ----" 

#Get the stat
# $cpuprog stat -s -p -T 18 -t ${NTHREAD} ${RUN}.cpu  1>${RUN}.stat # xyh
$cpuprog stat -s -p -T 18 -t 1 ${RUN}.cpu  1>${RUN}.stat
echo "--- statistics done  ----" 

echo "--- pigziping   ----" 
date 
$pigz -p ${NTHREAD} ${RUN}.singlelinker.paired.fastq 
$pigz -p ${NTHREAD} ${RUN}.singlelinker.single.fastq 
$pigz -p ${NTHREAD} ${RUN}.none.fastq 
$pigz -p ${NTHREAD} ${RUN}.conflict.fastq 
$pigz -p ${NTHREAD} ${RUN}.tied.fastq 

mv ${RUN}.tied.* ${RUN}.conflict.* $dumpfolder

echo "$jobname completed:"
date 


ENDHERE

jidLinker=$(submitjob $pbs)

echo "Passing jobID $jidLinker of $jobname"
echo "cpuLinker $jobname $jidLinker $pbs.log" >&2
buildCBlist "$fqpref" "$jidarc:$jidLinker"
cpuMapping $fqpref $pairlabel $jidLinker
cpuMapping $fqpref $singlabel $jidLinker
cpuMapping $fqpref $nonelabel $jidLinker
echo ""
}

function passLinker {
echo "passLinker cmd -- $1 "
#RUN THIS if cpuLinker is completed OK!
    fqpref=$1
export RUN=$(basename $fqpref)
export cpudir=$maindir/$RUN/bulkChiatac
echo "Linker must have been detected under $cpudir"

cpuMapping $fqpref $pairlabel
sleep 7
cpuMapping $fqpref $singlabel
sleep 7
cpuMapping $fqpref $nonelabel
echo ""
}


#------------------------------------------------------------------------------------------------------------
function mapARC10x {
echo "mapARC10x cmd -- $1 -- $2 "
#run 10x pipeline for multiome

export fqpref=$1
export fqRNApref=$2
export RUN=$(basename $fqpref)
export gexRUN=$(basename $fqRNApref)
fqDir=$(dirname $fqpref)
gexDir=$(dirname $fqRNApref)


jobname=${RUN}.arc
pbs=${jobname}.slurm

export wdir=$maindir/$RUN
libcsv=$wdir/$jobname.csv
echo 'fastqs,sample,library_type' > "$libcsv"
echo "${fqDir},${RUN},Chromatin Accessibility" >> "$libcsv"
echo "${gexDir},${gexRUN},Gene Expression" >> "$libcsv"

echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n ${NTHREAD10x}  # number of threads
#SBATCH --mem ${GB}G  # memory pool for all cores # XYH high Memory
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$wdir"
# Software commands are supplied by the configuration.

date

echo "ARC 10x mapping for $RUN" 
echo "Reference genome: $ARCref" 

$arc10x --reference $ARCref --id 10xARC --libraries=$libcsv --localcores=$NTHREAD10x --jobmode=local # xyh
# $arc10x --reference $ARCref --id 10xARC --libraries=$libcsv --localcores=$NTHREAD10x --jobmode=slurm # xyh
# $arc10x --reference $ARCref --id 10xARC --libraries=$libcsv --localcores=$NTHREAD10x --jobmode=slurm --maxjobs 99 # xyh

date 

ENDHERE

jidarc=$(submitjob $pbs)

export jidarc
echo "mapARC10x $jidarc $pbs.log" >&2
# The linker stage schedules the CB dictionary after both inputs are ready.
echo "."
}

#------------------------------------------------------------------------------------------------------------
function buildCBlist {
echo "buildCBlist cmd -- $1 -- $2 "
export fqpref=$1
export RUN=$(basename $fqpref)
local prevJobid=$2

jobname=${RUN}.buildCBpy
pbs=${jobname}.slurm

export wdir=$maindir/$RUN
echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n 1  # number of threads
#SBATCH --mem ${cbdictGB}G  # memory pool for all cores # XYH high Memory
#SBATCH --partition=${highpart}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$wdir"
# Software commands are supplied by the configuration.

date 

echo "Running "$buildCBpy" " 

$python "$buildCBpy" $maindir/$RUN/10xARC/outs/atac_possorted_bam.bam $maindir/$RUN/bulkChiatac/${RUN}_R1.fastq.gz $RUN

echo "CB dictionary built:" 
/bin/ls -ltr ${RUN}*.pkl 
date 

ENDHERE

jidpkl=$(submitjob $pbs $prevJobid)

export jidpkl
echo "buildCBlist $jidpkl $pbs.log" >&2
echo "Recorded jobID $jidpkl of $jobname"
echo ""
}


#------------------------------------------------------------------------------------------------------------
function mapGEX10x {
#run 10x pipeline for scRNA
echo "mapGEX10x cmd -- $1 -- $2 "

export fqpref=$1
export fqRNApref=$2
export gexRUN=$(basename $fqRNApref)
export RUN=$(basename $fqpref)
gexDir=$(dirname $fqRNApref)

jobname=${gexRUN}.gex
pbs=${jobname}.slurm

export wdir=$maindir/$RUN

echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD10x  # number of threads
#SBATCH --mem ${GB}G # memory pool for all cores # XYH high Memory
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$wdir"

date 

echo "GEX 10x mapping for $RUN" 
echo "Reference genome: $GEXref" 
$gex10x --transcriptome $GEXref --id 10xGEX --sample ${gexRUN} --fastqs=$gexDir --localcores=$NTHREAD10x --localmem=$GB --include-introns --chemistry=$chemistry
$bamCoverage -b 10xGEX/outs/possorted_genome_bam.bam -o ${RUN}.GEX.RPKM.bedgraph --binSize 5 --outFileFormat bedgraph --normalizeUsing RPKM -p $NTHREAD10x
$pigz -p $NTHREAD10x ${RUN}.GEX.RPKM.bedgraph

ENDHERE

jidgex=$(submitjob "$pbs")

export jidgex
runid=$jidgex
echo "mapGEX10x $runid $pbs.log" >&2
echo "$pbs job is submitted..."
echo "."
}

#------------------------------------------------------------------------------------------------------------
function cpuMapping {
echo "cpuMapping cmd -- $1 -- $2 -- $3 "
export fqpref=$1
export RUN=$(basename $fqpref)
export LABEL=$2 #either pairlabel, singlabel, or nonelabel
local prevJobid=$3 #sbatch dependency

export wdir=$maindir/$RUN

jobname=${RUN}.map.$LABEL
pbs=${jobname}.slurm

echo "Submitting $pbs under $PWD"
export cpudir=$maindir/$RUN/bulkChiatac

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD  # number of threads
#SBATCH --mem ${MAPPING_GB}G # memory pool for all cores # XYH some Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$cpudir"

export genome=$genome
echo "Reference genome:"
/bin/ls ${fasta}*

echo  START  ${jobname} cpu memaln .. 
date
$cpuprog memaln -T $minMapScore -t ${NTHREAD} $fasta ${RUN}.${LABEL}.fastq.gz 1>${RUN}.${LABEL}.sam
$pigz -p ${NTHREAD} ${RUN}.${LABEL}.sam 

mv ${RUN}.${LABEL}.fastq.gz $dumpfolder

echo "$jobname completed:" 
date 
ENDHERE

jidMap=$(submitjob $pbs $prevJobid)
echo "cpuMapping $jobname $jidMap $pbs.log" >&2
echo "Passing jobID $jidMap of $jobname"
cpuPairing $fqpref $LABEL $jidMap
statCB_RAW $fqpref $LABEL $jidMap
}

function passMapping {
echo "passMapping cmd -- $1 "
#RUN THIS if cpuLinker & cpuMapping are completed OK!
    fqpref=$1
export RUN=$(basename $fqpref)
export cpudir=$maindir/$RUN/bulkChiatac
echo "Mapping must have been finished under $cpudir"

cpuPairing $fqpref $pairlabel
sleep 7 
cpuPairing $fqpref $singlabel
sleep 7 
cpuPairing $fqpref $nonelabel
sleep 7 
statCB_RAW $fqpref $pairlabel
sleep 7 
statCB_RAW $fqpref $singlabel
sleep 7 
statCB_RAW $fqpref $nonelabel
echo "Exit from passMapping"
sleep 7 
}



#------------------------------------------------------------------------------------------------------------
function cpuPairing {
echo "cpuPairing cmd -- $1 -- $2 -- $3 "
export fqpref=$1
export RUN=$(basename $fqpref)
export LABEL=$2 #either pairlabel, singlabel, or nonelabel
local prevJobid=$3
export cpudir=$maindir/$RUN/bulkChiatac

UUbam=${RUN}.${LABEL}.UU.bam
UxxUbam=${RUN}.${LABEL}.UxxU.bam
jobname=${RUN}.pairing.$LABEL
pbs=${jobname}.slurm
echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD  # number of threads
#SBATCH --mem ${MAPPING_GB}G # memory pool for all cores # XYH some Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$cpudir"
echo  STARTED ${jobname}  
date 
# $cpuprog pair -s $selfbp -S -t ${NTHREAD} -q $minMapScore ${RUN}.${LABEL}.sam.gz 1>${RUN}.${LABEL}.stat.xls # xyh
$cpuprog pair -s $selfbp -S -t ${NTHREAD} -q $mapq ${RUN}.${LABEL}.sam.gz 1>${RUN}.${LABEL}.stat.xls 
/bin/ls -ltr *.sam* 
echo "$jobname completed:" 
date 

ENDHERE
jidPairing=$(submitjob $pbs $prevJobid)
echo "cpuPairing $jobname $jidPairing $pbs.log" >&2
echo "Passing jobID $jidPairing of $jobname"
dedupCB  $fqpref $LABEL $jidPairing
}

function passPairing {
echo "passPairing cmd -- $1 "
#if cpuPairing for all LABEL types (pairlabel, singlabel, nonelabel) were done!
export fqpref=$1
export RUN=$(basename $fqpref)
export cpudir=$maindir/$RUN/bulkChiatac

echo "Pairing were completed previously, continue to dedupCB "
dedupCB  $fqpref $pairlabel 
sleep 7 
dedupCB  $fqpref $singlabel 
sleep 7 
dedupCB  $fqpref $nonelabel 
echo "Exit from passPairing job submission function..."
sleep 7 
}

#------------------------------------------------------------------------------------------------------------
function cpuCluster {
echo "cpuCluster cmd -- $1 -- $2 "
export fqpref=$1
export RUN=$(basename $fqpref)
local prevJobid=$2

export cpudir=$maindir/$RUN/bulkChiatac

local cis_prefix="${RUN}.e500.clusters.cis.chiasig"
local juicefile=${RUN}.e500.juice
jobname=${RUN}.cpuCluster.${pairlabel}
pbs=${jobname}.slurm

nrBAM=${RUN}.${pairlabel}.UU.F2048.CB.nr.nsorted.bam
echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD  # number of threads
#SBATCH --mem ${MAPPING_GB}G # memory pool for all cores # XYH some Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$cpudir"

echo  'Computing span of paired tags started...' 
date 

# deduplicated span
echo  Computing span of paired tags UU nr... 
date
$cpuprog span -s $bulkSelfbp -t ${NTHREAD} $nrBAM  1>${RUN}.$pairlabel.UU.nr.span.xls

# cluster tags
echo  STARTED ${RUN} cpu clustering... self ligation removed from juicer...   
$cpuprog cluster -s $bulkSelfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${RUN}.e500 -j -x -v 1 -g ${nrBAM} 
echo  ENDED ${RUN} cpu clustering...

echo  Creating BE2 file and clean with $blacklist 
zcat ${cis_prefix}.gz | awk '{ if ( \$7 > 1 ) print }' | $bedcmd pairtobed -type neither -a stdin -b ${blacklist}  > ${cis_prefix}.BE2

$scom ChiaSig -c $minPETsignif -t $NTHREAD -p -m $bulkSelfbp ${cis_prefix}.BE2

zcat ${juicefile}.gz | sed "s/chr//g" > ${juicefile}
/bin/rm ${juicefile}.gz
$samcmd view -H ${nrBAM}| grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? | sed 's/chr//' > ${RUN}.length.txt
sort -k1 -V ${RUN}.length.txt | grep -v M >  ${RUN}.chrlength.txt
grep M ${RUN}.length.txt >> ${RUN}.chrlength.txt
rm ${RUN}.length.txt
$juicertool pre --threads $NTHREAD -r 2500000,1000000,500000,250000,100000,50000,25000,10000 -k KR,GW_KR ${juicefile} $RUN.hic ${RUN}.chrlength.txt
gzip -f ${juicefile} # xyh
gzip -f ${cis_prefix}.BE2 & # xyh
gzip -f *chiasig # xyh


wait
echo "$jobname completed:" 
date 

mv ${cis_prefix}.BE2.gz ${RUN}.*.BE2.*.chiasig.gz $dumpfolder
ENDHERE

jidChiasig=$(submitjob $pbs $prevJobid)
echo "cpuCluster $jobname $jidChiasig $pbs.log" >&2
echo 'Chiasig will be called.'
echo '..'
#echo "Passing jobID $jidChiasig of $jobname"
#bulk_chiapet_Stat $fqpref $jidChiasig


}

#--------------------------------------------------------------------------------------------------------
function dedupCB {
echo "dedupCB cmd -- $1 -- $2 -- $3 "

export fqpref=$1
export RUN=$(basename $fqpref)
export LABEL=$2 #either pairlabel, singlabel, or nonelabel
local prevJobid=$3
export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac
export pkl=$maindir/$RUN/${RUN}.CB_forCPU.pkl

jobname=${RUN}.dedupCB.$LABEL
pbs=${jobname}.slurm
UUbam=${RUN}.${LABEL}.UU.bam
UxxUbam=${RUN}.${LABEL}.UxxU.bam
picksam=${RUN}.${LABEL}.UxxU.picked.sam.gz

echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $cbNthread  # number of threads
#SBATCH --mem ${GB}G # memory pool for all cores # XYH high Memory
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$cpudir"
# Software commands are supplied by the configuration.

date

echo  STARTED ${jobname}   
# $samcmd view -F 2048 -F 4 -h $UxxUbam |$python "$pickUpy" | $pigz -p $cbNthread  > $picksam # xyh
$samcmd view -F 2048 -F 4 -h $UxxUbam | $python "$pickUpy" ${mapq} | $pigz -p $cbNthread  > $picksam
# Add CB & ZB tags to UU.bam to for generating fragments.tsv
if [ '$LABEL' == '$pairlabel' ]; then 
    echo "Get CB in the $UUbam"  
    date  
    $samcmd  view -F 2048 -@ $cbNthread -h $UUbam | $python "$addCBpy" $pkl > ${RUN}.${pairlabel}.UU.F2048.CB.sam
    echo "CB injected to ${RUN}.${pairlabel}.UU.F2048.CB.sam "
    echo "count nrPET CB in ${RUN}.${pairlabel}.UU.F2048.CB.sam"
    date
    $samcmd view  ${RUN}.${pairlabel}.UU.F2048.CB.sam | $python "$statpy_bamTag"  > ${RUN}.${pairlabel}.UU.F2048.CB.PET.tsv

    echo "Continue to sort & dedup ..."
    date
    grep -v noCB ${RUN}.${pairlabel}.UU.F2048.CB.sam | $samcmd sort -@ $cbNthread -o ${RUN}.${pairlabel}.UU.F2048.CB.csorted.bam
    $dedup -M ${RUN}.${pairlabel}.UU.F2048.CB.nr.MarkDuplicates.txt -O ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam  -I ${RUN}.${pairlabel}.UU.F2048.CB.csorted.bam

    echo "Query sort then count nrPET CB in ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam"
    date
    $samcmd sort -n -@ $cbNthread -o ${RUN}.${pairlabel}.UU.F2048.CB.nr.nsorted.bam ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam  #nsorted for cpu cluster then ChiaSig
    $samcmd view -@ $cbNthread ${RUN}.${pairlabel}.UU.F2048.CB.nr.nsorted.bam | $python "$statpy_bamTag"  > ${RUN}.${pairlabel}.UU.F2048.CB.nrPET.tsv

    echo "Generating interaction ${RUN}.${pairlabel}.UU.F2048.CB.nr.bedpe"
    date
    # $samcmd sort -n -@ $cbNthread ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam | $samcmd view -h | perl "$injectCBreadID" | $bedcmd  bamtobed -bedpe |sort -k7,7 > ${RUN}.${pairlabel}.UU.F2048.CB.nr.bedpe  # xyh
    $samcmd view -h ${RUN}.${pairlabel}.UU.F2048.CB.nr.nsorted.bam | perl "$injectCBreadID" | $bedcmd  bamtobed -bedpe |sort -k7,7 > ${RUN}.${pairlabel}.UU.F2048.CB.nr.bedpe

    echo "Get CB in the $UxxUbam"
    date
    $samcmd view -h $picksam | $python "$addCBpy" $pkl | grep -v noCB  | $samcmd sort -@ $cbNthread -o ${RUN}.${LABEL}.UxxU.F2048.CB.csorted.bam
    $dedup -M ${RUN}.${pairlabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt -O ${RUN}.${pairlabel}.UxxU.F2048.CB.nr.bam  -I ${RUN}.${pairlabel}.UxxU.F2048.CB.csorted.bam

    echo "Get CBstats  $UxxUbam"
    date
    $samcmd view $picksam  | $python "$statpy_noTag" $pkl > ${RUN}.${LABEL}.UxxU.F2048.CB.tsv
    $samcmd view ${RUN}.${LABEL}.UxxU.F2048.CB.nr.bam | $python "$statpy_noTag" $pkl  > ${RUN}.${LABEL}.UxxU.F2048.CB.nr.tsv 

    $pigz -p $cbNthread ${RUN}.${pairlabel}.UU.F2048.CB.sam
    $pigz -p $cbNthread ${RUN}.${pairlabel}.UU.F2048.CB.nr.bedpe
    echo "Interaction CB:"
    ls -ltr *bedpe*
fi

if [ '$LABEL' == '$singlabel' ]; then 
    echo "Get CB in the $UxxUbam"
    date
    $samcmd view -h $picksam | $python "$addCBpy" $pkl | grep -v noCB  | $samcmd sort -@ $cbNthread -o ${RUN}.${LABEL}.UxxU.F2048.CB.csorted.bam
    echo "Dedup starts "
    date
    $dedup -M ${RUN}.${singlabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt -O ${RUN}.${singlabel}.UxxU.F2048.CB.nr.bam  -I ${RUN}.${singlabel}.UxxU.F2048.CB.csorted.bam
    echo "Get CBstats  $UxxUbam"
    date
    $samcmd view $picksam  | $python "$statpy_noTag" $pkl > ${RUN}.${LABEL}.UxxU.F2048.CB.tsv
    $samcmd view ${RUN}.${LABEL}.UxxU.F2048.CB.nr.bam | $python "$statpy_noTag" $pkl  > ${RUN}.${LABEL}.UxxU.F2048.CB.nr.tsv 
fi

if [ '$LABEL' == '$nonelabel' ]; then 
    echo "Get CB in the $UUbam"
    date  
    $samcmd  view -F 2048 -@ $cbNthread -h $UUbam | $python "$addCBpy" $pkl > ${RUN}.${nonelabel}.UU.F2048.CB.sam
    echo "CB injected to ${RUN}.${nonelabel}.UU.F2048.CB.sam "
    echo "count nrPET CB in ${RUN}.${nonelabel}.UU.F2048.CB.sam"
    date
    $samcmd view ${RUN}.${nonelabel}.UU.F2048.CB.sam | $python "$statpy_bamTag" > ${RUN}.${nonelabel}.UU.F2048.CB.PET.tsv

    echo "Continue to sort & dedup ..."
    date
    grep -v noCB ${RUN}.${nonelabel}.UU.F2048.CB.sam | $samcmd sort -@ $cbNthread -o ${RUN}.${nonelabel}.UU.F2048.CB.csorted.bam
    $dedup -M ${RUN}.${nonelabel}.UU.F2048.CB.nr.MarkDuplicates.txt -O ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam  -I ${RUN}.${nonelabel}.UU.F2048.CB.csorted.bam
    $samcmd view -F 2048 -F 4 ${RUN}.${LABEL}.UU.F2048.CB.nr.bam | $python "$statpy_noTag" $pkl  > ${RUN}.${LABEL}.UU.F2048.CB.nr.tsv
    $samcmd view -F 2048 -F 4 ${RUN}.${LABEL}.UU.F2048.CB.nr.bam | $python "$statpy_bamTag" > ${RUN}.${nonelabel}.UU.F2048.CB.nr.PET.tsv
    $samcmd  view -F 2048 -@ $cbNthread  $UUbam | $python "$statpy_noTag" $pkl  > ${RUN}.${LABEL}.UU.F2048.CB.tsv 

    echo "Query sort then count nrPET CB in ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam"
    date
    $samcmd sort -n -@ $cbNthread -o ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.bam ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam
    $python ${filtersmallsize} ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.bam ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted
    $samcmd view -@ $cbNthread ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam | $python "$statpy_bamTag"  > ${RUN}.${nonelabel}.UU.F2048.CB.pet.nr.PET.tsv

    echo "Generating interaction ${RUN}.${nonelabel}.UU.F2048.CB.nr.bedpe"
    date
    $samcmd view -h ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam | perl "$injectCBreadID" | $bedcmd  bamtobed -bedpe | sort -k7,7 > ${RUN}.${nonelabel}.UU.F2048.CB.nr.bedpe

    echo "Get CB in the $UxxUbam"
    date
    $samcmd view -h $picksam | $python "$addCBpy" $pkl | grep -v noCB  | $samcmd sort -@ $cbNthread -o ${RUN}.${LABEL}.UxxU.F2048.CB.csorted.bam
    $dedup -M ${RUN}.${nonelabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt -O ${RUN}.${nonelabel}.UxxU.F2048.CB.nr.bam  -I ${RUN}.${nonelabel}.UxxU.F2048.CB.csorted.bam

    echo "Get CBstats  $UxxUbam"
    date
    $samcmd view $picksam  | $python "$statpy_noTag" $pkl > ${RUN}.${LABEL}.UxxU.F2048.CB.tsv
    $samcmd view ${RUN}.${LABEL}.UxxU.F2048.CB.nr.bam | $python "$statpy_noTag" $pkl  > ${RUN}.${LABEL}.UxxU.F2048.CB.nr.tsv 

    $pigz -p $cbNthread ${RUN}.${nonelabel}.UU.F2048.CB.sam
    $pigz -p $cbNthread ${RUN}.${nonelabel}.UU.F2048.CB.nr.bedpe
    echo "Interaction CB:"
    ls -ltr *bedpe*
fi


echo "$jobname completed:"
date
ENDHERE

jidCBdedup=$(submitjobOptional8 $pbs $prevJobid $jidpkl)
echo "dedupCB $jobname $jidCBdedup $pbs.log" >&2
echo "jobID $jidCBdedup submitted for $jobname"
if [ "$LABEL" == "$pairlabel" ]; then
    echo "Saving jobID:"
    storeJOBid $jidCBdedup
    echo $storedIDtsv
    echo "Calling cpuCluster job and will wait for jobID $jidCBdedup to complete."
    cpuCluster $fqpref $jidCBdedup
fi

storeDedupJOBid1 $jidCBdedup
echo "dedupCB jobids: $dedupIDs1"
echo '..'


}

#------------------------------------------------------------------------------------------------------------
function statCB_RAW {
echo "statCB_RAW cmd -- $1 -- $2 -- $3 "

export fqpref=$1
export RUN=$(basename $fqpref)
export LABEL=$2 #either pairlabel, singlabel, or nonelabel
local prevJobid=$3
export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac
export pkl=$maindir/$RUN/${RUN}.CB_forCPU.pkl

jobname=${RUN}.statCB_RAW.$LABEL
pbs=${jobname}.slurm

echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of threads
#SBATCH --mem ${cbdictGB}G # memory pool for all cores # XYH high Memory
#SBATCH --partition=${highpart}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "$cpudir"
# Software commands are supplied by the configuration.

date
SECONDS=0

echo  STARTED $jobname 
$samcmd view ${RUN}.${LABEL}.sam.gz | $python "$statpy_noTag"  $pkl  > ${RUN}.${LABEL}.statCB.totalPET.tsv 
echo "$jobname completed in \$SECONDS secs."
date


ENDHERE

jidnoTag=$(submitjobOptional8 $pbs $prevJobid $jidpkl)
echo "statCB_RAW $jobname $jidnoTag $pbs.log" >&2
echo "$pbs queued and save jobID: $jidnoTag"
storeJOBid $jidnoTag
echo $storedIDtsv
echo "."
}

#------------------------------------------------------------------------------------------------------------

function combine_tsv {

echo "combine_tsv cmd -- $1 -- $2 -- $3 -- $4 -- $5 -- $6 "

export fqpref=$1
export RUN=$(basename $fqpref)

local j1=$2
local j2=$3
local j3=$4
local j4=$5
local j5=$6
local j6=$7
local j7=$8

export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac

jobname=${RUN}.combine_tsv
pbs=${jobname}.slurm
echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n 1  # number of threads
#SBATCH --mem ${TABLE_GB}G # memory pool for all cores # XYH little Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
# module load R-3.6.1
cd "$cpudir"

date
SECONDS=0

echo  STARTED $jobname 
$Rcmd "$mergeTable" ${RUN} # xyh change Rscript
echo "$jobname completed in \$SECONDS secs."
date

ENDHERE

jid_tsv=$(submitjobOptional8  $pbs $j1 $j2 $j3 $j4 $j5 $j6 $j7)

export jid_tsv
echo "combine_tsv $jobname $jid_tsv $pbs.log" >&2
echo "$pbs queued and wait for $j1 $j2 $j3 $j4 $j5 to finish."
echo '...'

}


#------------------------------------------------------------------------------------------------------------

function mergeFinalBAM {
echo "mergeFinalBAM cmd -- $1 -- $2 -- $3 -- $4 "

export fqpref=$1
export RUN=$(basename $fqpref)
local j1=$2
local j2=$3
local j3=$4

export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac
export jobdir=$maindir/$RUN/slurm
export partTHREAD=$(( (NTHREAD + 1) / 2 ))

jobname=${RUN}.mergeFinalBAM
pbs=${jobname}.slurm
tsvfile=${RUN}.fragments.tsv
finalSuf=$RUN.F2048.CBreadID.nr.final

echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n $NTHREAD  # number of threads
#SBATCH --mem ${MERGE_GB}G # memory pool for all cores # XYH some Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
cd "${cpudir}"

date 
SECONDS=0

echo "Get one tag from ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam" 
date 
# $samcmd view -@ $NTHREAD -hb -f 64 -o ${RUN}.${nonelabel}.UU.F2048.f64.CB.nr.bam ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam # xyh 
# ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam
# ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam
$samcmd view -H ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam | grep '^@SQ' | cut -f 2-3 |grep -v '_chr' | sed s?SN:?? | sed s?LN:??   > ${genome}.chrlength.txt
# $python "$Sbam2fragment" ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam | $samcmd sort -@ $NTHREAD -o ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam.tmp -
$samcmd sort -@ $NTHREAD -o ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam.tmp ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam
$samcmd view -@ $NTHREAD -Sb -f 64 ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam | $samcmd sort -@ $NTHREAD -o ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam.tmp -

echo  STARTED $jobname 
# $samcmd merge -@ $NTHREAD -o ${RUN}.tmp.bam ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam ${RUN}.${pairlabel}.UxxU.F2048.CB.nr.bam ${RUN}.${singlabel}.UxxU.F2048.CB.nr.bam ${RUN}.${nonelabel}.UxxU.F2048.CB.nr.bam ${RUN}.${nonelabel}.UU.F2048.f64.CB.nr.bam # xyh
$samcmd merge -@ $NTHREAD -o ${RUN}.tmp.bam ${RUN}.${pairlabel}.UU.F2048.CB.nr.bam ${RUN}.${pairlabel}.UxxU.F2048.CB.nr.bam ${RUN}.${singlabel}.UxxU.F2048.CB.nr.bam ${RUN}.${nonelabel}.UxxU.F2048.CB.nr.bam ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.singleton.bam.tmp ${RUN}.${nonelabel}.UU.F2048.CB.nr.nsorted.pet.bam.tmp
echo  "Merging nr.bam files  done in \$SECONDS secs."  
date 

$bedcmd intersect -v -nonamecheck -a ${RUN}.tmp.bam -b $blacklist > ${RUN}.nobl.bam

SECONDS=0
echo  "Injecting CB to readID ...  "  
$samcmd view -h -@ ${partTHREAD} ${RUN}.nobl.bam  | perl "$injectCBreadID" | $samcmd view -@ ${partTHREAD} -o ${finalSuf}.bam -
echo  "done in \$SECONDS secs."  

#Generate bigwig for visualization
$bedcmd genomecov -ibam ${finalSuf}.bam -bg > ${finalSuf}.bedgraph
$scom bedSort ${finalSuf}.bedgraph  ${finalSuf}.sorted.bedgraph

# $samcmd view -H ${RUN}.${nonelabel}.UU.F2048.CB.nr.bam | grep '^@SQ' | cut -f 2-3 |grep -v '_chr' | sed s?SN:?? | sed s?LN:??   > ${genome}.chrlength.txt # xyh
$scom bedGraphToBigWig  ${finalSuf}.sorted.bedgraph  ${genome}.chrlength.txt ${finalSuf}.bw
$bamCoverage -b ${finalSuf}.bam -o ${RUN}.ATAC.RPKM.bedgraph --binSize 5 --outFileFormat bedgraph --normalizeUsing RPKM -p $NTHREAD
echo -e  "Done converting to bigwig n" 
date 

#Generating fragment file
SECONDS=0
echo "Starting to get fragments.tsv from ${finalSuf}.bam" 
date  
# perl "$bamtofrag" --read_file ${finalSuf}.bam --output_file ${tsvfile} --samtools_path $toolsbin # xyh
$python "$bamtofrag" ${finalSuf}.bam ${tsvfile}
# perl "$bamtofrag" --read_file ${finalSuf}.bam --output_file ${tsvfile}.tmptmp --samtools_path $toolsbin # xyh
# mv ${tsvfile}.tmptmp.len ${tsvfile}.len # xyh
# awk -v OFS="\t" '{if(\$2<0){\$2=0}print \$0}' ${tsvfile}.tmptmp > ${tsvfile} # xyh
# rm ${tsvfile}.tmptmp # xyh
echo  "done in \$SECONDS secs."  
date

$bgzip -f ${tsvfile}.len # xyh

$Rcmd --vanilla "$sortfrag" ${tsvfile} 

$bgzip -f ${tsvfile} # xyh
$tabix -f -p bed ${tsvfile}.gz # xyh

echo "Call peaks using ${RUN}.forpeak.rmbl.bed"  
$bedcmd bamtobed -i ${finalSuf}.bam | $bedcmd intersect -v -nonamecheck -a stdin -b $blacklist > ${RUN}.forpeak.rmbl.bed
date
$scom macs2 callpeak --keep-dup all --nomodel --shift $shiftsize  --extsize $peakext  -B --SPMR -t ${RUN}.forpeak.rmbl.bed  -f BED -g $gmac -n $RUN  --qvalue $macsq 
echo macs2 finished 
date 
gzip -f ${RUN}_summits.bed # xyh

$scom bedClip ${RUN}_treat_pileup.bdg ${genome}.chrlength.txt ${RUN}_treat_pileup.clip.bdg
# sort -k1,1 -k2,2n ${RUN}_treat_pileup.clip.bdg > ${RUN}_treat_pileup.clip.sorted.bdg # xyh
$scom bedSort ${RUN}_treat_pileup.clip.bdg ${RUN}_treat_pileup.clip.sorted.bdg
$scom bedGraphToBigWig ${RUN}_treat_pileup.clip.sorted.bdg ${genome}.chrlength.txt  ${RUN}.treat_pileup.NDP.bw
echo ENDED ${RUN} coverage density generated 

echo "$jobname is complete." 

date 

$pigz -p ${NTHREAD}  *forpeak.rmbl.bed
$pigz -p ${NTHREAD}  *.bedgraph
$pigz -p ${NTHREAD}  *.bdg

ENDHERE

jm=$(submitjobOptional8  $pbs $j1 $j2 $j3)
echo "combine_tsv $jobname $jm $pbs.log" >&2
export jidmerge=$jm
storeMergeJOBid $jm
echo "$pbs queued and wait for $j1 $j2 $j3 to finish."
echo '...'

}


function summaryJob {
echo "summaryJob cmd -- $1 -- $2 -- $3 -- $4 "

export fqpref=$1
export RUN=$(basename $fqpref)
local j1=$2
local j2=$3
local j3=$4

export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac

jobname=${RUN}.summaryJob
pbs=${jobname}.slurm
echo "Submitting $pbs under $PWD"

cat << ENDHERE > $pbs
#!/usr/bin/env bash

#SBATCH -N 1 # number of nodes
#SBATCH -n 1  # number of threads
#SBATCH --mem ${SUMMARY_GB}G # memory pool for all cores # XYH little Memory
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -eo pipefail
SECONDS=0

echo  STARTED $jobname 
date

#Summary including peaks
cd ${cpudir} 
"$sumsh" $RUN $genome

echo "Organize files under $PWD"
#Clearing files
# /bin/rm ${RUN}.*.tmp.bam ${RUN}*clip.*bdg* ${RUN}.*.bedgraph*  ${RUN}.*.forpeak.* ${RUN}.*picked.sam* ${RUN}.nobl.bam # xyh
####xyh /bin/rm ${RUN}.*.tmp.bam ${RUN}*clip.*bdg* ${RUN}.*.forpeak.* ${RUN}.*picked.sam*
####xyh mv ${RUN}.cpu ${RUN}.*.bam ${RUN}.*.sam.gz ${RUN}.*.CB.*tsv ${RUN}.*.totalPET.tsv ${genome}.*.txt ${RUN}.*.rmbl.bed* ${RUN}*summits.bed* ${RUN}*.bdg* $dumpfolder
####xyh mv $dumpfolder/*final.bam .
cd $dumpfolder
#/bin/rm *.bam


echo FINISHED $jobname 
date

ENDHERE

jid_clean=$(submitjobOptional8  $pbs $j1 $j2 $j3)

export jid_clean
echo "summaryJob $jobname $jid_clean $pbs.log" >&2
echo "$pbs queued and wait for $j1 $j2 $j3 to finish."
echo '...'

}

function cleanfiles {
echo "cleanfiles cmd -- $1 "

#No need to submit slurm job for this one, just run manually at the end after checking everything is completed.
export fqpref=$1
export RUN=$(basename $fqpref)

export wdir=$maindir/$RUN
export cpudir=$maindir/$RUN/bulkChiatac

cd "$cpudir"
echo "Organize files under $PWD"

#Clearing files
mv *.cpu *.bam *.sam.gz ${RUN}.*.CB.*tsv $dumpfolder

/bin/rm ${RUN}.*.tmp.bam ${RUN}*clip.*bdg* ${RUN}.*.bedgraph*  ${RUN}.*.forpeak.* ${RUN}.*picked.sam*
mv *.bam *.CB.sam* *.tsv ${genome}.*.txt  $dumpfolder
mv $dumpfolder/*final.bam .
cd $dumpfolder
# /bin/rm *.bam

echo 'Clearing done...'

}


#------------------------------------------------------------------------------------------------------------




function preparePETStats {
    local RUN=$(basename "$1")
    local cpudir=$maindir/$RUN/bulkChiatac
    local pbs=$maindir/$RUN/slurm/${RUN}.PETstats.slurm
    shift
    cat << ENDHERE > "$pbs"
#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem ${TABLE_GB}G
#SBATCH --partition=${smallQUEUE}
#SBATCH --output="${pbs}.out"
#SBATCH --error="${pbs}.err"
set -eo pipefail
cd "$cpudir"
# Reuse the barcode/name separation and PET counter from the author workflow.
zcat ${RUN}.${pairlabel}.UU.F2048.CB.nr.bedpe.gz ${RUN}.${nonelabel}.UU.F2048.CB.nr.bedpe.gz | \
    perl -lane '@tmp=split(/:/,\$F[6],2); print join("\t",@F[0..5],@tmp,@F[7..9])' | \
    $python "$CountBarcodePET" /dev/stdin ${RUN}.hq.CB.nr.bedpe.stat
ENDHERE
    jidpetstats=$(submitjob "$pbs" "$@")
    export jidpetstats
}

function reanalyzeATAC {
    local RUN=$(basename "$1")
    local wdir=$maindir/$RUN
    local cpudir=$wdir/bulkChiatac
    local pbs=$wdir/slurm/${RUN}.ATACreanalyze.slurm
    shift
    cat << ENDHERE > "$pbs"
#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n ${NTHREAD10x}
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --output="${pbs}.out"
#SBATCH --error="${pbs}.err"
set -eo pipefail
cd "$wdir"
awk '\$7+0 > $peak_signal_min' "$cpudir/${RUN}_peaks.narrowPeak" | grep -v -i chrM | cut -f 1-3 > "$cpudir/${RUN}.peak"
$atac10x --id 10xATAC --peaks "$cpudir/${RUN}.peak" --fragments "$cpudir/${RUN}.fragments.tsv.gz" --reference "$ATACref" --localcores ${NTHREAD10x}
ENDHERE
    jidatac=$(submitjob "$pbs" "$@")
    export jidatac
}

function analyzeCells {
    local RUN=$(basename "$1")
    local wdir=$maindir/$RUN
    local pbs=$wdir/slurm/${RUN}.callcell.slurm
    shift
    cat << ENDHERE > "$pbs"
#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n ${NTHREAD}
#SBATCH --mem ${DOWNSTREAM_GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --output="${pbs}.out"
#SBATCH --error="${pbs}.err"
set -eo pipefail
export CHAIR_MIN_RNA=$CHAIR_MIN_RNA CHAIR_MIN_ATAC=$CHAIR_MIN_ATAC CHAIR_MIN_PET=$CHAIR_MIN_PET
export CHAIR_MAX_MT=$CHAIR_MAX_MT CHAIR_MAX_RIBO=$CHAIR_MAX_RIBO CHAIR_DOUBLET_RATE=$CHAIR_DOUBLET_RATE CHAIR_SEED=$CHAIR_SEED
export CHAIR_RNA_DIMS=$CHAIR_RNA_DIMS CHAIR_ATAC_DIMS=$CHAIR_ATAC_DIMS
export CHAIR_RNA_RESOLUTION=$CHAIR_RNA_RESOLUTION CHAIR_ATAC_RESOLUTION=$CHAIR_ATAC_RESOLUTION CHAIR_WNN_RESOLUTION=$CHAIR_WNN_RESOLUTION
$Rcmd --vanilla "$tmppath/callcell_xingyu.R" "$RUN" "$wdir/10xGEX/outs/raw_feature_bc_matrix.h5" "$wdir/10xATAC/outs/raw_peak_bc_matrix.h5" "$wdir/bulkChiatac/${RUN}.fragments.tsv.gz" "$wdir/bulkChiatac/${RUN}.hq.CB.nr.bedpe.stat" "$wdir/CallCell" "$genome"
ENDHERE
    jidcells=$(submitjob "$pbs" "$@")
    export jidcells
}

#------- loop & submit

#Below shows the example for complete run from the begining:
#Input arg:  fqlist
#Format of fqlist is one library per row consists of  2 columns
#Modify below if needed to run specific tasks
#(Don't worry about slurm job IDs [$jid_tsv, $jidMerge, $dedupIDs], they will be populated automatically)
while read -r dna rna extra; do
   [[ -z "$dna" || "$dna" == \#* ]] && continue
   [[ -n "$rna" && -z "$extra" ]] || die "Expected two FASTQ prefixes per row"
   dna=$(realpath "$dna")
   rna=$(realpath "$rna")
   dedupIDs1=""; dedupIDs2=""; storedIDtsv=""; jidMerge=""
   jidpkl=""; jidgex=""; jidChiasig=""; jidmerge=""
echo "read -- $dna -- $rna "
   RUN=$(basename $dna)   
   export wdir=$maindir/$RUN
   export cpudir=$wdir/bulkChiatac
   export jobdir=$wdir/slurm
   cd "$maindir"
   odir $RUN
   odir "$jobdir"
   odir "$cpudir" && odir "$cpudir/$dumpfolder" 
   cd "$jobdir"
   mapARC10x "$dna" "$rna"
   mapGEX10x "$dna" "$rna"
   cpuLinker "$dna"
   mergeFinalBAM "$dna" $dedupIDs1
   combine_tsv "$dna" $dedupIDs1 $storedIDtsv
   preparePETStats "$dna" $dedupIDs1
   reanalyzeATAC "$dna" "$jidmerge"
   analyzeCells "$dna" "$jidgex" "$jidatac" "$jidpetstats"
   summaryJob "$dna" "$jid_tsv" "$jidmerge" "$jidChiasig"
   echo "All jobs for $dna submitted"
   date
   echo ''
   echo '-------------------------------------------'
   echo ''
done < "$fqlist" 



