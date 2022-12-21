# read and set the version

#from .input_commands import get_input, instantiate_dirs, validate_fastq

from .bam import sam_to_bam, bam_to_mapped_or_unmapped, extract_long_fastq, extract_short_fastq
from .case_one import case_one
from .case_three import case_three
from .cleanup import remove_intermediate_files, move_and_copy_files, touch_file, touch_output_fail_files
from .concat import concatenate_all_fastqs, concatenate_single
from .deduplicate import deduplicate_fastqs
from .depth import get_depth, concatenate_chrom_plasmids, get_contig_lengths, get_contig_circularity, bwa_map_depth_sort, minimap_depth_sort, get_depths_from_bam, collate_depths, combine_depth_dfs
from .extract import extract_chromosome, extract_chromosome_fasta, extract_plasmid_fastas
from .flye import run_flye, contig_count
from .input_commands import get_input, instantiate_dirs, validate_fastq
from .logging import write_to_log
from .mapping import index_fasta, minimap_long_reads, bwa_map_short_reads
from .qc import nanofilt, trim_short_read
from .main import run
from .unicycler import run_unicycler



__version__ = "0.1.1"


# #from .processes import write_to_log, run_flye, contig_count, extract_chromosome, extract_chromosome_fasta, extract_plasmid_fastas, nanofilt, trim_short_read, plasmid_assembly,double_mapping_analysis,index_fasta,minimap_long_reads, bwa_map_short_reads, sam_to_bam, bam_to_mapped_or_unmapped, extract_long_fastq, extract_short_fastq, concatenate_fastqs,concatenate_single,deduplicate_fastqs,unicycler, extract_reads_mapping_to_plasmid_and_chromosome, remove_intermediate_files, move_and_copy_files,touch_file,touch_output_fail_files 

# from .processes import run_minimap2, keep_primary_supplementary_mappings_convert_sam, keep_primary_supplementary_mappings_convert_sam_direct, run_flagstat
# from .post_processing import get_total_read_count, parse_bam, pivot_df
# from .main import run
# from .version import __version__

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
            'run_unicycler'
           ]