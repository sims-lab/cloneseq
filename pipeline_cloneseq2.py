#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRNASeq2 pipeline to extract barcodes from sequencing data.
"""

import sys
import os
import pathlib
import pandas
import pysam

from ruffus import *

from cgatcore import experiment as E
from cgatcore import iotools as IOTools
from cgatcore import pipeline as P

P.get_parameters('pipeline.yml')

@follows(mkdir("vector.dir"))
@merge(
    input=P.PARAMS["filename_vector_fasta"],
    output="vector.dir/vector.fa"
)
def index_vector_sequence(infile, outfile):
    """Convert vector sequene (viral genome + inserts) to BWA file format."""
    statement = (
        f"cp {infile} {outfile} && "
        f"samtools faidx {outfile} && "
        f"bwa index {outfile} "
        f"> {outfile}.log ")
    return P.run(statement)


@transform(
    input=index_vector_sequence,
    filter=suffix(".fa"),
    output=".bed.gz"
)
def build_motif_bed(infile, outfile):
    """Create a BED file with the positions of the variable bases (i.e. 'N') in the barcode."""
    with pysam.FastxFile(infile) as inf:
        vector = next(inf)

    contig = vector.name
    positions = [x for x, c in enumerate(vector.sequence) if c == "N"]

    with IOTools.open_file(P.snip(outfile, ".gz"), "w") as outf:
        for pos in positions:
            outf.write("{}\t{}\t{}\n".format(contig, pos, pos + 1))
    pysam.tabix_index(P.snip(outfile, ".gz"), preset="bed")


@merge(
    input=index_vector_sequence,
    output="vector.dir/barcode.fa"
)
def build_motif_fasta(infile, outfile):
    """Extract the barcode sequence."""
    with pysam.FastxFile(infile) as inf:
        vector_sequence = next(inf).sequence

    start = vector_sequence.find("N")
    end = vector_sequence.rfind("N") + 1

    with IOTools.open_file(outfile, "w") as outf:
        outf.write(">barcode\n{}\n".format(vector_sequence[start:end]))

@follows(mkdir("fastq.dir"))
@collate(
    input=P.PARAMS["input_fastq_glob"],
    filter=formatter(r"(?P<SAMPLE>[^/]+).fastq.[12].gz"),
    output=["fastq.dir/{SAMPLE[0]}.matched.fastq.1.gz",
            "fastq.dir/{SAMPLE[0]}.matched.fastq.2.gz",
            "fastq.dir/{SAMPLE[0]}.unmatched.fastq.1.gz",
            "fastq.dir/{SAMPLE[0]}.unmatched.fastq.2.gz"],
    add_inputs=(index_vector_sequence,))
def filter_read_data(infiles, outfiles):
    """Filter reads for occurrence of the barcode motif."""

    fastq1, fastq2 = sorted([os.path.abspath(x[0]) for x in infiles])
    vector_fasta = infiles[0][1]
    outprefix = IOTools.snip(outfiles[0], ".matched.fastq.1.gz")

    filter_options = P.get_params().get("kmer_filtering", {})
    filtering_kmer_size=filter_options.get("kmer_size", 10)
    filtering_min_kmer_matches=filter_options.get("min_kmer_matches", 10)

    statement = (
        f"cgat fastqs2fastqs "
        f"--input-filename-fasta={vector_fasta} "
        f"--method=filter-by-sequence "
        f"--force-output "
        f"--output-filename-pattern={outprefix}.%%s "
        f"--filtering-kmer-size={filtering_kmer_size} "
        f"--filtering-min-kmer-matches={filtering_min_kmer_matches} "
        f"{fastq1} "
        f"{fastq2} "
        f"> {outprefix}.log "
    )

    return P.run(statement)


@merge(
    input=filter_read_data,
    output="filtering_summary.tsv"
)
def summarize_filtering(infiles, outfile):
    tables = " ".join([P.snip(x[0], ".matched.fastq.1.gz") + ".log" for x in infiles])
    statement = (
        f"cgat combine-tables "
        f"--log={outfile}.log "
        f"--cat sample "
        f"--regex-filename=\".dir/([^/.]+)\" "
        f"{tables} "
        f"> {outfile}"
    )
    return P.run(statement)


@follows(mkdir("mapped.dir"))
@collate(
    input=filter_read_data,
    filter=formatter(
        r"fastq.dir/(?P<SAMPLE>[^/]+).matched.fastq.[12].gz"),
    output="mapped.dir/{SAMPLE[0]}.bam",
    add_inputs=(index_vector_sequence,),
)
def align_vector_sequences(infiles, outfile):
    """Align reads containing the barcode."""
    fastq_files = sorted([x for x in infiles[0][0] if ".matched" in x])
    vector_fasta = infiles[0][1]
    read_group = P.snip(os.path.basename(fastq_files[0]), ".fastq.1.gz")
    # high clipping penalties, low gap penalties
    # high match score to encourage positions mapping to
    # the non-barcode sequences
    fastq = " ".join(fastq_files)
    statement = (
        f"bwa mem -k14 -W20 -r10 -A1 -B1 -O2 -E1 -L100 -U10 "
        f"-R '@RG\\tID:{read_group}\\tSM:{read_group}' "
        f"{vector_fasta} {fastq} "
        f"2> {outfile}.bwa.err "
        f"| samtools sort - "
        f"2> {outfile}.sort.err "
        f"| samtools view -bS "
        f"2> {outfile}.view.err "
        f"> {outfile} && "
        f"samtools index {outfile}"
    )
    return P.run(statement)


@follows(mkdir("merged_bam.dir"))
@collate(
    input=align_vector_sequences,
    filter=formatter(r"mapped.dir/(?P<SAMPLE>[^/]+)-lane(\S+).bam"),
    output="merged_bam.dir/{SAMPLE[0]}.bam",
)
def merge_lanes(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        "samtools merge -f {outfile} {infiles} "
        "2> {outfile}.merge.err && "
        "samtools index {outfile} ".format(**locals()))
    return P.run(statement)


@follows(mkdir("deduplicated.dir"))
@transform(
    input=[align_vector_sequences, merge_lanes],
    filter=formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
    output="deduplicated.dir/{SAMPLE[0]}.bam",
)
def remove_duplicates(infile, outfile):
    job_memory = "4G"
    picard_opts = '-Xmx{} -XX:+UseConcMarkSweepGC'.format(job_memory)
    statement = (
        f"picard {picard_opts} "
        f"MarkDuplicates "
        f"INPUT={infile} "
        f"OUTPUT={outfile} "
        f"METRICS_FILE={outfile}.metrics "
        f">& {outfile}.log && "
        f"samtools index {outfile}"
    )
    return P.run(statement, job_memory="unlimited")


@follows(mkdir("samtools_bamstats.dir"))
@transform(
    input=remove_duplicates,
    filter=formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
    output="samtools_bamstats.dir/{SAMPLE[0]}.tsv",
)
def compute_samtools_bamstats(infile, outfile):
    statement = f"samtools stats {infile} > {outfile}"
    return P.run(statement)


@follows(mkdir("samtools_depth.dir"))
@transform(
    input=remove_duplicates,
    filter=formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
    output="samtools_depth.dir/{SAMPLE[0]}.tsv",
)
def compute_samtools_depth(infile, outfile):
    statement = (
        f"echo -e \"contig\\tpos\\tdepth\" > {outfile} && "
        f"samtools depth {infile} "
        f">> {outfile} ")
    return P.run(statement)


@merge(
    input=compute_samtools_depth,
    output="samtools_depth.tsv"
)
def summarize_samtools_depth(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        f"cgat combine-tables "
        f"--log={outfile}.log "
        f"--cat sample "
        f"--regex-filename=\".dir/([^/.]+)\" "
        f"{infiles} "
        f"> {outfile} "
    )
    return P.run(statement)


@follows(mkdir("cgat_bamstats.dir"))
@transform(
    input=remove_duplicates,
    filter=formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
    output="cgat_bamstats.dir/{SAMPLE[0]}.tsv",
)
def compute_cgat_bamstats(infile, outfile):
    statement = (
        f"cgat bam2stats {infile} "
        f"--force-output "
        f"--output-filename-pattern={outfile}.%%s "
        f"--log={outfile}.log "
        f"> {outfile} ")
    return P.run(statement)


@merge(
    input=compute_cgat_bamstats,
    output="cgat_bamstats.tsv"
)
def summarize_cgat_bamstats(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        "cgat combine-tables "
        "--log={outfile}.log "
        "--cat sample "
        "--regex-filename=\".dir/([^/.]+)\" "
        "{infiles} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


@follows(mkdir("clone_codes.dir"))
@transform(
    input=remove_duplicates,
    filter=regex(".dir/([^/]+).bam"),
    output=r"clone_codes.dir/\1.tsv",
    add_inputs=(index_vector_sequence, build_motif_bed),
)
def extract_clone_codes_pileup(infiles, outfile):
    """Determine the exact barcode sequence present in each cell."""
    bamfile, fafile, bedfile = infiles

    statement = (
        f"cgat bam-pileup2tsv "
        f"--input-bed={bedfile} "
        f"--reference-fasta={fafile} "
        f"--min-base-quality=0 "
        f"--method=barcode "
        f"--log={outfile}.pileup.log "
        f"{bamfile} "
        f"> {outfile}.pileup "
    )
    P.run(statement, job_memory="8G")

    full_df = pandas.read_csv(outfile + ".pileup", sep="\t")

    with IOTools.open_file(outfile, "w") as outf:
        headers = full_df.consensus_support.describe().index
        full_df["consensus"] = full_df.consensus.fillna("N")
        consensus_sequence = "".join(full_df.consensus)

        eval_df = full_df.copy()
        median_consensus_depth = eval_df.consensus_counts.median()

        # modules to recover partial bar-codes
        # 1.: bases missing from the bar-code if otherwise median depth is high
        # => ignore these bases for depth and other QC metrics
        if median_consensus_depth > 2:
            deleted_barcode_bases = full_df[(full_df.consensus == "N") & (full_df.depth == 0)]
            eval_df = eval_df[~eval_df.pos.isin(deleted_barcode_bases.pos)]
            deleted_barcode_bases = deleted_barcode_bases.index
        else:
            deleted_barcode_bases = ""

        outf.write("\t".join(map(str, ["barcode", "ndeleted_barcode_bases", "deleted_barcode_bases"] +
                                 ["support_{}".format(x) for x in headers] +
                                 ["counts_{}".format(x) for x in headers] +
                                 ["offcounts_{}".format(x) for x in headers])) + "\n")
        outf.write("\t".join(map(str, [
            consensus_sequence,
            len(deleted_barcode_bases),
            ",".join(map(str, deleted_barcode_bases))] +
                                 eval_df.consensus_support.describe().tolist() +
                                 eval_df.consensus_counts.describe().tolist() +
                                 eval_df.offconsensus_counts.describe().tolist())) + "\n")

@merge(
    input=extract_clone_codes_pileup,
    output="clone_codes_pileup.tsv"
)
def summarize_clone_codes_pileup(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        f"cgat combine-tables "
        f"--log={outfile}.log "
        f"--cat sample "
        f"--regex-filename=\".dir/([^/.]+)\" "
        f"{infiles} "
        f"> {outfile}"
    )
    return P.run(statement)

# # NB: the multialignment branch seemed dead in Andreas' version of the pipeline.
# # TODO: summarize_clone_codes_mali does not handle empty TSV files when concatenating tables
# @follows(mkdir("mali_clone_codes.dir"))
# @transform(
#     input=remove_duplicates,
#     filter=regex(".dir/([^/]+).bam"),
#     output=r"mali_clone_codes.dir/\1.tsv",
#     add_inputs=(index_vector_sequence, build_motif_bed, build_motif_fasta),
# )
# def extract_clone_codes_mali(infiles, outfile):
#     """Multi-align regions of the reads containing the barcode (+ an
#     anchor, as there may be deletions or mutations within the barcode
#     sequence).
#     """
#     bamfile, fafile, bedfile, vectorfile = infiles

#     statement = (
#         f"cgat bam2fasta "
#         f"--input-bed-file={bedfile} "
#         f"--reference-fasta={fafile} "
#         f"--merge-intervals "
#         f"--output-filename-pattern={outfile}.%%s "
#         f"--barcode-fasta={vectorfile} "
#         f"--log={outfile}.log "
#         f"{bamfile} "
#         f"2> {outfile}.err "
#         f"> {outfile}"
#     )
#     P.run(statement)


# @merge(
#     input=extract_clone_codes_mali,
#     output="clone_codes_mali.tsv"
# )
# def summarize_clone_codes_mali(infiles, outfile):
#     infiles = " ".join(infiles)
#     statement = (
#         f"cgat combine-tables "
#         f"--log={outfile}.log "
#         f"--cat sample "
#         f"--regex-filename=\".dir/([^/.]+)\" "
#         f"{infiles} "
#         f"> {outfile}"
#     )
#     return P.run(statement)


@merge(
    input=[summarize_clone_codes_pileup,
           summarize_cgat_bamstats,
           summarize_filtering],
    output="summary.tsv"
)
def summarize_all(infiles, outfile):
    """merge all summary tables into a single table."""
    fn_clonecodes, fn_bamstats, fn_filtering = infiles
    df_clonecodes = pandas.read_csv(fn_clonecodes, sep="\t", na_values="na")
    df_bamstats = pandas.read_csv(fn_bamstats, sep="\t")
    df_filterstats = pandas.read_csv(fn_filtering, sep="\t")

    # aggregate the per-lane filter-stats and append to filter stats
    df_filterstats["lane"] = df_filterstats["sample"].str.extract("(lane\d+)", expand=False)
    extra = df_filterstats.loc[~df_filterstats["lane"].isnull(), :].copy()
    extra["sample"] = extra["sample"].str.extract("(.*)-lane\d+", expand=False)
    extra = extra.groupby(["sample"]).sum().reset_index()
    df_filterstats = pandas.concat([df_filterstats, extra], axis=0)
    E.info("adding lane-combined metrics to filter_stats: extra={}, filterstats={}".format(
        len(extra), len(df_filterstats)))

    df_filterstats["percent_matched"] = 100.0 * df_filterstats["matched"] / df_filterstats["input"]

    # normalize bamstats df
    df_bamstats = pandas.pivot_table(df_bamstats, index="sample", columns="category", values="counts").reset_index()

    # sanitize clone codes df
    df_clonecodes["barcode"] = df_clonecodes.barcode.fillna("NNNNNNNNN")
    df_clonecodes["counts_min"] = df_clonecodes.counts_min.fillna(0)
    df_clonecodes["support_min"] = df_clonecodes.support_min.fillna(0)

    # merge
    df_merged = df_clonecodes.merge(df_bamstats, left_on="sample", right_on="sample")\
                             .merge(df_filterstats, left_on="sample", right_on="sample")
    E.info("number of rows bamstats: {}".format(len(df_bamstats)))
    E.info("number of rows clonecodes: {}".format(len(df_clonecodes)))
    E.info("number of rows filterstats: {}".format(len(df_filterstats)))
    E.info("number of rows merged: {}".format(len(df_merged)))
    df_merged.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
