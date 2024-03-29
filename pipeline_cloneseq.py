"""===========================
Pipeline template
===========================

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_cloneseq.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""

import sys
import os
import ruffus
import pandas
import pysam

import cgatcore.experiment as E
from cgatcore import iotools as IOTools
from cgatcore import pipeline as P


def index_vector_sequence(infile, outfile):
    """Convert vector sequene (viral genome + inserts) to BWA file format."""
    statement = (
        "cp {infile} {outfile} && "
        "samtools faidx {outfile} && "
        "bwa index {outfile} "
        "> {outfile}.log ".format(**locals()))
    return P.run(statement)


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


def build_motif_fasta(infile, outfile):
    """Extract the barcode sequence."""
    with pysam.FastxFile(infile) as inf:
        vector_sequence = next(inf).sequence

    start = vector_sequence.find("N")
    end = vector_sequence.rfind("N") + 1

    with IOTools.open_file(outfile, "w") as outf:
        outf.write(">barcode\n{}\n".format(vector_sequence[start:end]))


def filter_read_data(infiles, outfiles):
    """Filter reads for occurrence of the barcode motif."""

    fastq1, fastq2 = sorted([os.path.abspath(x[0]) for x in infiles])
    vector_fasta = infiles[0][1]
    outprefix = IOTools.snip(outfiles[0], ".matched.fastq.1.gz")

    filter_options = P.get_params().get("kmer_filtering", {})

    statement = (
        "cgat fastqs2fastqs "
        "--input-filename-fasta={vector_fasta} "
        "--method=filter-by-sequence "
        "--force-output "
        "--output-filename-pattern={outprefix}.%%s "
        "--filtering-kmer-size={filtering_kmer_size} "
        "--filtering-min-kmer-matches={filtering_min_kmer_matches} "
        "{fastq1} "
        "{fastq2} "
        "> {outprefix}.log ".format(
            filtering_kmer_size=filter_options.get("kmer_size", 10),
            filtering_min_kmer_matches=filter_options.get("min_kmer_matches", 10),
            **locals()))

    return P.run(statement)


def summarize_filtering(infiles, outfile):

    tables = " ".join([P.snip(x[0], ".matched.fastq.1.gz") + ".log" for x in infiles])
    statement = (
        "cgat combine-tables "
        "--log={outfile}.log "
        "--cat sample "
        "--regex-filename=\".dir/([^/.]+)\" "
        "{tables} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


def align_vector_sequences(infiles, outfile):
    """Align reads containing the barcode."""
    fastq_files = sorted([x for x in infiles[0][0] if ".matched" in x])
    vector_fasta = infiles[0][1]
    read_group = P.snip(os.path.basename(fastq_files[0]), ".fastq.1.gz")
    # high clipping penalties, low gap penalties
    # high match score to encourage positions mapping to
    # the non-barcode sequences
    statement = (
        "bwa mem -k14 -W20 -r10 -A1 -B1 -O2 -E1 -L100 -U10 "
        "-R '@RG\\tID:{read_group}\\tSM:{read_group}' "
        "{vector_fasta} {fastq} "
        "2> {outfile}.bwa.err "
        "| samtools sort - "
        "2> {outfile}.sort.err "
        "| samtools view -bS "
        "2> {outfile}.view.err "
        "> {outfile} && "
        "samtools index {outfile}".format(
            fastq=" ".join(fastq_files),
            **locals()))
    return P.run(statement)


def merge_lanes(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        "samtools merge -f {outfile} {infiles} "
        "2> {outfile}.merge.err && "
        "samtools index {outfile} ".format(**locals()))
    return P.run(statement)


def remove_duplicates(infile, outfile):

    job_memory = "4G"
    # picard_opts = '-Xmx{} -XX:+UseParNewGC -XX:+UseConcMarkSweepGC'.format(job_memory)
    picard_opts = '-Xmx{} -XX:+UseConcMarkSweepGC'.format(job_memory)

    statement = (
        "picard {picard_opts} "
        "MarkDuplicates "
        "INPUT={infile} "
        "OUTPUT={outfile} "
        "METRICS_FILE={outfile}.metrics "
        ">& {outfile}.log && "
        "samtools index {outfile}".format(**locals()))
    return P.run(statement, job_memory="unlimited")


def compute_samtools_bamstats(infile, outfile):
    statement = (
        "samtools stats {infile} "
        "> {outfile} ".format(**locals()))
    return P.run(statement)


def compute_samtools_depth(infile, outfile):
    statement = (
        "echo -e \"contig\\tpos\\tdepth\" > {outfile} && "
        "samtools depth {infile} "
        ">> {outfile} ".format(**locals()))
    return P.run(statement)


def summarize_samtools_depth(infiles, outfile):
    infiles = " ".join(infiles)
    statement = (
        "cgat combine-tables "
        "--log={outfile}.log "
        "--cat sample "
        "--regex-filename=\".dir/([^/.]+)\" "
        "{infiles} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


def compute_cgat_bamstats(infile, outfile):
    statement = (
        "cgat bam2stats {infile} "
        "--force-output "
        "--output-filename-pattern={outfile}.%%s "
        "--log={outfile}.log "
        "> {outfile} ".format(**locals()))
    return P.run(statement)


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


def extract_clone_codes_pileup(infiles, outfile):
    """Determine the exact barcode sequence present in each cell."""
    bamfile, fafile, bedfile = infiles

    statement = (
        "cgat bam-pileup2tsv "
        "--input-bed={bedfile} "
        "--reference-fasta={fafile} "
        "--min-base-quality=0 "
        "--method=barcode "
        "--log={outfile}.pileup.log "
        "{bamfile} "
        "> {outfile}.pileup ".format(**locals()))
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


def summarize_clone_codes_pileup(infiles, outfile):

    infiles = " ".join(infiles)
    statement = (
        "cgat combine-tables "
        "--log={outfile}.log "
        "--cat sample "
        "--regex-filename=\".dir/([^/.]+)\" "
        "{infiles} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


def extract_clone_codes_mali(infiles, outfile):
    """Multi-align regions of the reads containing the barcode (+ an
    anchor, as there may be deletions or mutations within the barcode
    sequence).
    """
    bamfile, fafile, bedfile, vectorfile = infiles

    statement = (
        "cgat bam2fasta "
        "--input-bed-file={bedfile} "
        "--reference-fasta={fafile} "
        "--merge-intervals "
        "--output-filename-pattern={outfile}.%%s "
        "--barcode-fasta={vectorfile} "
        "--log={outfile}.log "
        "{bamfile} "
        "2> {outfile}.err "
        "> {outfile}".format(**locals()))
    P.run(statement)


def summarize_clone_codes_mali(infiles, outfile):

    infiles = " ".join(infiles)
    statement = (
        "cgat combine-tables "
        "--log={outfile}.log "
        "--cat sample "
        "--regex-filename=\".dir/([^/.]+)\" "
        "{infiles} "
        "> {outfile}".format(**locals()))
    return P.run(statement)


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


def main(argv=None):
    if argv is None:
        argv = sys.argv

#    def add_input(parser):
#        parser.add_option(
#            "--filename-vector-fasta",
#            dest="filename_vector_fasta",
#            default=None,
#            type="string",
#            help="filename of vector sequence in fasta format "
#            "[default=%default].")
#
#        parser.add_option(
#            "--input-fastq-glob",
#            dest="input_fastq_glob",
#            default=None,
#            type="string",
#            help="glob expression for paired-end read data in fastq format "
#            "[default=%default].")
#
#    options, args = P.parse_commandline(argv,
#                                        config_file="pipeline.yml",
#                                        callback=add_input)
    options, args = P.parse_commandline(argv, config_file="pipeline.yml")

#    if options.filename_vector_fasta is None:
#        raise ValueError("please specify --vector-fasta")
#
#    if options.input_fastq_glob is None:
#        raise ValueError("please specify --input-fastq-glob")

#    if options.config_file:
#        P.get_parameters(options.config_file)
#    else:
#        sys.exit(P.main(options, args))
    P.get_parameters("pipeline.yml")

    pipeline = ruffus.Pipeline("cgatflow-cloneseq")

    task_index_vector_sequence = pipeline.merge(
        task_func=index_vector_sequence,
#        input=options.filename_vector_fasta,
        input=P.get_params().get("filename_vector_fasta"),
        output="vector.dir/vector.fa").mkdir("vector.dir")

    task_build_motif_bed = pipeline.transform(
        task_func=build_motif_bed,
        input=task_index_vector_sequence,
        filter=ruffus.suffix(".fa"),
        output=".bed.gz")

    task_build_motif_fasta = pipeline.merge(
        task_func=build_motif_fasta,
        input=task_index_vector_sequence,
        output="vector.dir/barcode.fa")

    task_filter_read_data = pipeline.collate(
        task_func=filter_read_data,
#        input=options.input_fastq_glob,
        input=P.get_params().get("input_fastq_glob"),
        filter=ruffus.formatter(r"(?P<SAMPLE>[^/]+).fastq.[12].gz"),
        output=["fastq.dir/{SAMPLE[0]}.matched.fastq.1.gz",
                "fastq.dir/{SAMPLE[0]}.matched.fastq.2.gz",
                "fastq.dir/{SAMPLE[0]}.unmatched.fastq.1.gz",
                "fastq.dir/{SAMPLE[0]}.unmatched.fastq.2.gz"],
        add_inputs=(task_index_vector_sequence,),
        ).mkdir("fastq.dir")

    task_summarize_filtering = pipeline.merge(
        task_func=summarize_filtering,
        input=task_filter_read_data,
        output="filtering_summary.tsv")

    task_align_vector_sequences = pipeline.collate(
        task_func=align_vector_sequences,
        input=task_filter_read_data,
        filter=ruffus.formatter(
            r"fastq.dir/(?P<SAMPLE>[^/]+).matched.fastq.[12].gz"),
        output="mapped.dir/{SAMPLE[0]}.bam",
        add_inputs=(task_index_vector_sequence,),
        ).mkdir("mapped.dir")

    task_merge_lanes = pipeline.collate(
        task_func=merge_lanes,
        input=task_align_vector_sequences,
        filter=ruffus.formatter(r"mapped.dir/(?P<SAMPLE>[^/]+)-lane(\S+).bam"),
        output="merged_bam.dir/{SAMPLE[0]}.bam",
        ).mkdir("merged_bam.dir")

    task_remove_duplicates = pipeline.transform(
        task_func=remove_duplicates,
        input=[task_align_vector_sequences, task_merge_lanes],
        filter=ruffus.formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
        output="deduplicated.dir/{SAMPLE[0]}.bam",
        ).mkdir("deduplicated.dir")

    task_compute_cgat_bamstats = pipeline.transform(
        task_func=compute_cgat_bamstats,
        input=task_remove_duplicates,
        filter=ruffus.formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
        output="cgat_bamstats.dir/{SAMPLE[0]}.tsv",
        ).mkdir("cgat_bamstats.dir")

    task_summarize_cgat_bamstats = pipeline.merge(
        task_func=summarize_cgat_bamstats,
        input=task_compute_cgat_bamstats,
        output="cgat_bamstats.tsv")

    task_compute_samtools_bamstats = pipeline.transform(
        task_func=compute_samtools_bamstats,
        input=task_remove_duplicates,
        filter=ruffus.formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
        output="samtools_bamstats.dir/{SAMPLE[0]}.tsv",
        ).mkdir("samtools_bamstats.dir")

    task_compute_samtools_depth = pipeline.transform(
        task_func=compute_samtools_depth,
        input=task_remove_duplicates,
        filter=ruffus.formatter(r".dir/(?P<SAMPLE>[^/]+).bam"),
        output="samtools_depth.dir/{SAMPLE[0]}.tsv",
        ).mkdir("samtools_depth.dir")

    task_summarize_samtools_depth = pipeline.merge(
        task_func=summarize_samtools_depth,
        input=task_compute_samtools_depth,
        output="samtools_depth.tsv")

    task_extract_clone_codes_pileup = pipeline.transform(
        task_func=extract_clone_codes_pileup,
        input=task_remove_duplicates,
        filter=ruffus.regex(".dir/([^/]+).bam"),
        output=r"clone_codes.dir/\1.tsv",
        add_inputs=(task_index_vector_sequence, task_build_motif_bed),
        ).mkdir("clone_codes.dir")

    task_extract_clone_codes_mali = pipeline.transform(
        task_func=extract_clone_codes_mali,
        input=task_remove_duplicates,
        filter=ruffus.regex(".dir/([^/]+).bam"),
        output=r"mali_clone_codes.dir/\1.tsv",
        add_inputs=(task_index_vector_sequence, task_build_motif_bed, task_build_motif_fasta),
        ).mkdir("mali_clone_codes.dir")

    task_summarize_clone_codes_pileup = pipeline.merge(
        task_func=summarize_clone_codes_pileup,
        input=task_extract_clone_codes_pileup,
        output="clone_codes_pileup.tsv")

    task_summarize_clone_codes_mali = pipeline.merge(
        task_func=summarize_clone_codes_mali,
        input=task_extract_clone_codes_mali,
        output="clone_codes_mali.tsv")

    task_summarize_all = pipeline.merge(
        task_func=summarize_all,
        input=[task_summarize_clone_codes_pileup,
               task_summarize_cgat_bamstats,
               task_summarize_filtering],
        output="summary.tsv")

    # primary targets
    pipeline.merge(
        task_func=P.EmptyRunner("all"),
        input=[
            task_merge_lanes,
            task_summarize_filtering,
            task_summarize_cgat_bamstats,
            task_compute_samtools_bamstats,
            task_summarize_samtools_depth,
            task_summarize_clone_codes_pileup,
            task_extract_clone_codes_mali,
#            task_summarize_clone_codes_mali,
            task_summarize_all,
        ],
        output="all")

    return P.run_workflow(options, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
