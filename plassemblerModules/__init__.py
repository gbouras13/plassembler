# read and set the version

#from .input_commands import get_input, instantiate_dirs, validate_fastq

from .log import write_to_log
from .bam import sam_to_bam, bam_to_mapped_or_unmapped, extract_long_fastq, extract_short_fastq
from .case_one import case_one
from .case_three import case_three
from .cleanup import remove_intermediate_files, move_and_copy_files, touch_file, touch_output_fail_files
from .concat import concatenate_all_fastqs, concatenate_single
from .deduplicate import deduplicate_fastqs
from .depth import get_depth, concatenate_chrom_plasmids, get_contig_lengths, get_contig_circularity, bwa_map_depth_sort, minimap_depth_sort, get_depths_from_bam, collate_depths, combine_depth_dfs
from .extract import extract_chromosome, extract_chromosome_fasta, extract_plasmid_fastas
from .run_flye import run_flye, contig_count
from .input_commands import get_input, instantiate_dirs, validate_fastq
from .mapping import index_fasta, minimap_long_reads, bwa_map_short_reads
from .qc import nanofilt, trim_short_read
from .main import run
from .run_unicycler import run_unicycler
from .run_mash import get_contig_count





__all__ = ['sam_to_bam', 'bam_to_mapped_or_unmapped', 'extract_long_fastq', 'extract_short_fastq',
            'case_one',
            'case_three',
            'remove_intermediate_files', 'move_and_copy_files', 'touch_file', 'touch_output_fail_files',
            'concatenate_all_fastqs', 'concatenate_single',
            'deduplicate_fastqs',
            'get_depth', 'concatenate_chrom_plasmids', 'get_contig_lengths', 'get_contig_circularity', 'bwa_map_depth_sort', 'minimap_depth_sort', 'get_depths_from_bam', 'collate_depths', 'combine_depth_dfs',
            'extract_chromosome', 'extract_chromosome_fasta', 'extract_plasmid_fastas',
            'run_flye', 'contig_count',
            'get_input', 'instantiate_dirs', 'validate_fastq',
            'write_to_log',
            'index_fasta', 'minimap_long_reads', 'bwa_map_short_reads',
            'nanofilt', 'trim_short_read',
            'run',
            'run_unicycler',
            'get_contig_count'
           ]