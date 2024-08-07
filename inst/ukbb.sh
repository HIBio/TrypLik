#!/bin/bash

# set -e -u -o pipefail

# enable interrupting all the background jobs at once
# trap 'kill $(jobs -p)' INT

# install ripgrep
$(apt-get install ripgrep)
echo "ripgrep version: $(rg --version)"

# install bwa
$(git clone https://github.com/lh3/bwa.git && cd bwa && make && cd ..)
bwadir="$(pwd)/bwa"
echo "BWA installed to $bwadir"
echo "$(bwa/bwa)"

# install TrypLik
$(Rscript --vanilla -e 'install.packages("/mnt/project/TrypLik_0.1.0.tar.gz", repos = NULL)')

# get consensus file and script location
inst="$(Rscript --vanilla -e 'cat(system.file(package = "TrypLik"))')"
echo "consensus file and scripts should be in $inst"

if [[ -f "$consensusdir/consensus.fa.amb" ]]; then
  echo "consensus file already indexed"
else
  echo "Processing consensus..."
  ./bwa/bwa index $inst/consensus.fa
fi

# check samtools
echo "using samtools:"
echo "$(samtools version)"
# echo "with htslib var"
# echo "$(HTS_PATH=/usr/lib/x86_64-linux-gnu/htslib samtools version)"

export HTS_PATH=/usr/lib/x86_64-linux-gnu/htslib

# allow permissions to scripts
$(chmod a+x $inst/*.sh)

# nCores=$(getconf _NPROCESSORS_ONLN)
# Number of LOGICAL CPUs (includes those reported by hyper-threading cores)
# Linux: Simply count the number of (non-comment) output lines from `lscpu -p`,
# which tells us the number of *logical* CPUs.
logicalCpuCount=$([ $(uname) = 'Darwin' ] &&
                       sysctl -n hw.logicalcpu_max ||
                       lscpu -p | egrep -v '^#' | wc -l)

# Number of PHYSICAL CPUs (cores).
# Linux: The 2nd column contains the core ID, with each core ID having 1 or
#        - in the case of hyperthreading - more logical CPUs.
#        Counting the *unique* cores across lines tells us the
#        number of *physical* CPUs (cores).
physicalCpuCount=$([ $(uname) = 'Darwin' ] &&
                       sysctl -n hw.physicalcpu_max ||
                       lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)
echo "Logical cores available: $logicalCpuCount"
echo "Physical cores available: $physicalCpuCount"
nWorkers=${3-4}
nthreads=${4-4}
echo "using $nWorkers workers each with $nthreads threads"

## input dir
indir="$1"
echo "Using $indir as input dir"

## output dir
outdir="${2-.}"
echo "Using $outdir as output dir"

reldir="$(pwd)"
echo "reldir: $reldir"
cd "$indir"

## WHICH DIRS TO PROCESS - START WITH JUST ONE
dirs="${5-10}/"
#dirs=$(ls -d */)

task() { # $1 = idWorker, $2 = asset
  echo "Worker $1: Asset '$2' START!"
  echo "Processing $2"
  # sleep $(( ($RANDOM % 4) + 3 ))
  $inst/preprocess_single.sh $2 $reldir/$outdir/$d $inst $nthreads
  $inst/count_tryptase_fast.sh $reldir/$outdir/${2/.cram/.sam} $reldir/$outdir/$d
  echo "    Worker $1: Asset '$2' DONE!"
}

worker() { # $1 = idWorker
  echo "Worker $1 GO!"
  idAsset=0
  for asset in "${listAssets[@]}"; do # ONLY RUN THE FIRST TEN FILES? USE [@]:0:10
    # split assets among workers (using modulo); each worker will go through
    # the list and select the asset only if it belongs to that worker
    (( idAsset % nWorkers == $1 )) && task $1 "$asset"
    (( idAsset++ ))
  done
  echo "    Worker $1 ALL DONE!"
}

# for d in $(ls -d */); do
for d in $dirs; do
  echo "Working on $d"
  if [[ -d "$reldir/$outdir/$d" ]]; then
    echo "output dir already exists"
  else
    echo "Creating $reldir/$outdir/$d"
    mkdir -p "$reldir/$outdir/$d"
  fi

  # array of assets
  listAssets=( $(ls $d*.cram) )
  for (( idWorker=0; idWorker<nWorkers; idWorker++ )); do
    # start workers in parallel, use 1 process for each
    worker $idWorker &
  done
  wait # until all workers are done

done

echo "Cleaning up bwa"
## clean up temp installs
cd
rm -fr $bwadir
rm -fr $outdir/bwa

echo "*** Aggregating data ***"
cd
cd "$reldir/$outdir/$dirs"
for i in *counts.txt ; do echo $i; cat $i; done > "../results_${dirs/\//}.txt"
# tar -czvf "results_${dirs/\//}.tar.gz" *.txt
rm -f *.counts.txt
cd ..
rm -fr $dirs

echo "*** ALL DONE ***"

exit 0
