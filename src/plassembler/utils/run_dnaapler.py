from pathlib import Path
from loguru import logger
from plassembler.utils.external_tools import ExternalTool


def run_dnaapler(threads, logdir, outdir):
    """runs dnaapler bulk
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    canu_plasmid_fasta: Path = Path(outdir) / "plasmids_canu.fasta"
    dnaapler_outdir: Path = Path(outdir) / "dnaapler"


    try:
        dnaapler = ExternalTool(
            tool="dnaapler bulk",
            input="",
            output="",
            params=f" -i {canu_plasmid_fasta} -m plasmid -o {dnaapler_outdir} -t {threads}",
            logdir=logdir,
            outfile="",
        )

        ExternalTool.run_tool(dnaapler, to_stdout=False)
    except Exception:
        logger.error(
            f"Dnaapler failed."
        )


def make_blastdb(canu_output_dir, logdir):
    canu_fasta: Path = Path(canu_output_dir) / "canu.contigs.fasta"
    db: Path = Path(canu_output_dir) / "db"

    makeblastdb = ExternalTool(
        tool="makeblastdb",
        input=f"-in {canu_fasta}",
        output=f"-out {db}",
        params="-dbtype nucl ",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(makeblastdb, to_stdout=False)


"""
some of this adapted from dnaapler 

https://github.com/gbouras13/dnaapler

"""


def run_blast(canu_output_dir, threads, logdir):
    canu_fasta: Path = Path(canu_output_dir) / "canu.contigs.fasta"
    blast_output: Path = Path(canu_output_dir) / "blast_output.txt"
    db: Path = Path(canu_output_dir) / "db"

    blast = ExternalTool(
        tool="blastn",
        input=f"-query {canu_fasta}",
        output=f"-out {blast_output}",
        params=f'-db {db} -evalue  1e-05 -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
        outfile="",
    )
    ExternalTool.run_tool(blast, to_stdout=False)


def parse_blast_output(canu_output_dir, threads, logdir):
    canu_fasta: Path = Path(canu_output_dir) / "canu.contigs.fasta"
    blast_output_file: Path = Path(canu_output_dir) / "blast_output.txt"
    db: Path = Path(canu_output_dir) / "db"

    blast = ExternalTool(
        tool="blastn",
        input=f"-query {input}",
        output=f"-out {blast_output_file}",
        params=f'-db {db} -evalue  1e-05 -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
        outfile="",
    )


def process_blast_output(canu_output_dir, outdir):
    """Processes

    :param input: input file
    :param blast_file: blast output file
    :param out_file: output file
    :return: blast_success: bool - whether a BLAST hit with a valid start codon was identified
    """

    blast_output_file: Path = Path(canu_output_dir) / "blast_output.txt"
    canu_fasta: Path = Path(canu_output_dir) / "canu.contigs.fasta"

    # define colnames
    col_list = [
        "qseqid",
        "qlen",
        "sseqid",
        "slen",
        "length",
        "qstart",
        "qend",
        "sstart",
        "send",
        "pident",
        "nident",
        "gaps",
        "mismatch",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    # read in the dataframe from BLAST
    try:
        blast_df = pd.read_csv(
            blast_output_file, delimiter="\t", index_col=False, names=col_list
        )
    except Exception:
        logger.error("There was an issue with parsing the BLAST output file.")

    # if the BLAST input is empty
    if isinstance(blast_df, pd.DataFrame) and blast_df.empty:
        logger.error("There were 0 BLAST hits. This must be a BLAST error.")

    # keep only the columns where the same contig is blasted against each other

    # Filter rows where qseqid is equal to sseqid - the ones fo BLASTing against themselves
    blast_df = blast_df[blast_df["qseqid"] == blast_df["sseqid"]]

    # read in the canu FASTA and save as a dictionary

    fasta_dict = {}
    i = 1
    for record in SeqIO.parse(canu_fasta, "fasta"):
        fasta_dict[record.id] = {
            "count": i,
            "sequence": str(record.seq),
            "dupe": False,  # stores the dupe status
            "start": 1,  # stores the start
            "end": len(record.seq),  # stores the end
        }
        i += 1

    for contig in fasta_dict.keys():
        tmp_df = blast_df[blast_df["qseqid"] == contig]
        # Sort by 'length' column in descending order
        tmp_df_sorted = tmp_df.sort_values(by="length", ascending=False)
        # get rid of the 100% match row
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["qlen"] != tmp_df_sorted["length"]]
        # more than 99% identical
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["pident"] > 99]
        # starts need to be < 100
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["qstart"] < 100]
        num_rows = tmp_df_sorted.shape[0]
        if num_rows == 0:  # where there is no dupe at all
            fasta_dict[contig]["dupe"] = False
            # exit
        else:
            # Get the first row
            first_row = tmp_df_sorted.iloc[0]
            # ensure the match is good
            if (
                first_row["length"] < 500
            ):  # less than 500bp repeat in the top hit - probably Insertion Seq not a real plasmid dupe - otherwise probably legit
                fasta_dict[contig]["dupe"] = False
            else:
                try:
                    # the repeat will be in the longest hit with qstart < 100  (usually 1 or very close to it)
                    # heuristic i need to check i guess
                    # just take until the next repeat element
                    # this has been filtered for prior
                    best_row = tmp_df_sorted.iloc[0]
                    fasta_dict[contig]["dupe"] = True
                    fasta_dict[contig]["start"] = best_row["qstart"]
                    fasta_dict[contig]["end"] = best_row["sstart"]

                    # if the query end is larger than the sstart - there is an overlap
                    # take 1 as the start and then the sstart as the end
                    # otherwise check for  concatenation (within 1000bp))
                    # otherwise exit just the whole plasmid
                    # if best_row["qend"] > best_row["sstart"]:
                    #     fasta_dict[contig]["dupe"] = True
                    #     fasta_dict[contig]["start"] = best_row["qstart"]
                    #     fasta_dict[contig]["end"] = best_row["sstart"]
                    # #elif (best_row["qend"] + 1000) > best_row[
                    #     #"sstart"
                    # #]:  # the longest match is likely to be a duplication
                    # else:
                    #     fasta_dict[contig]["dupe"] = True
                    #     fasta_dict[contig]["start"] = best_row["qstart"]
                    #     fasta_dict[contig]["end"] = best_row["sstart"]
                    # else:
                    #     fasta_dict[contig]["dupe"] = False
                except Exception:
                    logger.error("Flye not found. Please reinstall Plassembler.")

    # Create a list of SeqRecord objects
    records = []
    for entry_id, entry_data in fasta_dict.items():
        subsequence = entry_data["sequence"][
            entry_data["start"] - 1 : entry_data["end"]
        ]
        count = str(entry_data["count"])
        l = entry_data["end"]

        record = SeqRecord(
            seq=Seq(subsequence), id=count, description=f"{count} len={l}"
        )
        records.append(record)

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "plasmids_canu.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


# then BLAST output
# # then figure out for overlaps
# parse all blast hits as a dictionary
# need more than 1 hit (itself)
# if the blast hit is more than 50% and less than 90 % of contig length, take as duplicate
# elif the blast hit is more than 2000bp (lower could be e.g. IS element) and far away (more than 50% length away)
# then assume partial duplication too
# get all dupe regions this way
