#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 4 ]]; then
    echo "Usage: bash prepare_viewer_inputs.sh SAMPLE_DIR VIEWER_DIR GENES_GTF CONFIG" >&2
    exit 1
fi
sample_dir=$(realpath "$1")
mkdir -p "$2"
viewer_dir=$(realpath "$2")
gtf=$(realpath "$3")
config=$(realpath "$4")
source "$config"
code_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
sample=$(basename "$sample_dir")
: "${python:=python3}" "${samcmd:=samtools}" "${bedcmd:=bedtools}"
cells="$sample_dir/CallCell/data/final_barcodes.txt"
rna_bam="$sample_dir/10xGEX/outs/possorted_genome_bam.bam"
atac="$sample_dir/bulkChiatac/$sample.fragments.tsv.gz"
peaks="$sample_dir/bulkChiatac/$sample.peak"
for path in "$cells" "$rna_bam" "$atac" "$peaks" "$gtf"; do
    [[ -s "$path" ]] || { echo "Required input is missing or empty: $path" >&2; exit 1; }
done
read_table() { case "$1" in *.gz) gzip -cd "$1";; *) cat "$1";; esac; }

# RNA interval export reuses the author's alignment-to-interval script.
$samcmd view -h -F 4 "$rna_bam" | \
    awk 'NR==FNR{bc[$1]=1;next} /^@/{print;next} {for(i=12;i<=NF;i++) if($i~/^CB:Z:/ && substr($i,6) in bc){print;break}}' "$cells" - | \
    $samcmd view -b -o "$viewer_dir/RNA.barcode.bam" -
$python "$code_dir/bam2rna_intervals.py" "$viewer_dir/RNA.barcode.bam" "$viewer_dir/RNA.intervals.tsv.gz"
read_table "$atac" | awk 'NR==FNR{bc[$1]=1;next} $4 in bc' "$cells" - | gzip -c > "$viewer_dir/ATAC.intervals.tsv.gz"
cp "$peaks" "$viewer_dir/ATAC.peaks.bed"

# The same barcode/name separation is used by the pipeline PET-statistics stage.
shopt -s nullglob
paired_pet="$sample_dir/bulkChiatac/$sample.singlelinker.paired.UU.F2048.CB.nr.shifted.bedpe.gz"
[[ -s "$paired_pet" ]] || paired_pet="$sample_dir/bulkChiatac/$sample.singlelinker.paired.UU.F2048.CB.nr.bedpe.gz"
pet_files=("$paired_pet" "$sample_dir/bulkChiatac/$sample.none.UU.F2048.CB.nr.bedpe.gz")
for path in "${pet_files[@]}"; do [[ -s "$path" ]] || { echo "Missing PET input: $path" >&2; exit 1; }; done
gzip -cd "${pet_files[@]}" | \
    perl -lane '@tmp=split(/:/,$F[6],2); print join("\t",@F[0..5],@tmp,@F[7..9])' | \
    awk 'NR==FNR{bc[$1]=1;next} $7 in bc' "$cells" - | gzip -c > "$viewer_dir/PET.bedpe.gz"

# Fixed-anchor loop preparation and counting reuse the author's existing scripts.
gzip -cd "$viewer_dir/PET.bedpe.gz" | \
    awk '$1!="chrM" && $4!="chrM" && $1!="chrY" && $4!="chrY"' | gzip -c > "$viewer_dir/PET.loop_input.bedpe.gz"
$python "$code_dir/bedpe2ipet.py" "$viewer_dir/PET.loop_input.bedpe.gz" "$viewer_dir/PET.ipet.gz"
$python "$code_dir/loopCallingGivenAnchor.py" "$viewer_dir/ATAC.peaks.bed" "$viewer_dir/PET.ipet.gz" "$viewer_dir/PET.loops.tsv.gz"

# Gene BED uses the eight-column layout read by ChAIR-Viewer.
read_table "$gtf" | awk '$0!~/^#/ && ($3=="exon" || $3=="UTR")' | \
    perl -lane '$id="#";$id=$1 if(/transcript_id "(.+?)"/);$name="#";$name=$1 if(/gene_name "(.+?)"/);print join("\t",$F[0],$F[3]-1,$F[4],$id,$name,0,$F[6],$F[2])' > "$viewer_dir/genes.bed"

rna_cov="$sample_dir/$sample.GEX.RPKM.bedgraph"
[[ -s "$rna_cov" ]] || rna_cov="$rna_cov.gz"
atac_cov="$sample_dir/bulkChiatac/$sample.ATAC.RPKM.bedgraph.gz"
[[ -s "$atac_cov" ]] || atac_cov="${atac_cov%.gz}"
for path in "$rna_cov" "$atac_cov"; do [[ -s "$path" ]] || { echo "Missing RPKM track: $path" >&2; exit 1; }; done
{
    printf 'track\tpath\n'
    printf 'genes\tgenes.bed\nrna\tRNA.intervals.tsv.gz\natac\tATAC.intervals.tsv.gz\npeaks\tATAC.peaks.bed\npet\tPET.bedpe.gz\nloops\tPET.loops.tsv.gz\n'
    printf 'rna_rpkm\t%s\natac_rpkm\t%s\n' "$rna_cov" "$atac_cov"
} > "$viewer_dir/inputs.tsv"
echo "Viewer manifest: $viewer_dir/inputs.tsv"
