"""
VarProfiler: A comprehensive toolkit for analyzing sequence variability, 
satellite DNA, and centromeres in genomic data.
"""

__version__ = "0.1.0"
__author__ = "aglabx"
__email__ = "ad3002@gmail.com"

# Import main utilities for convenience
from .utils.utils import (
    reverse_complement,
    read_fasta,
    write_fasta,
    merge_adjacent_regions,
    get_tqdm
)

from .utils.motif_utils import (
    find_motif_positions,
    count_motif_in_sequence,
    mark_motifs_in_sequence,
    load_motifs_from_file
)

from .utils.centromere_utils import (
    score_centromere_likelihood,
    classify_confidence_level
)

from .utils.gff_utils import (
    read_gff_annotations,
    write_gff_header,
    write_gff_feature
)

__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    
    # Basic utilities
    "reverse_complement",
    "read_fasta",
    "write_fasta",
    "merge_adjacent_regions",
    "get_tqdm",
    
    # Motif utilities
    "find_motif_positions",
    "count_motif_in_sequence",
    "mark_motifs_in_sequence",
    "load_motifs_from_file",
    
    # Centromere utilities
    "score_centromere_likelihood",
    "classify_confidence_level",
    
    # GFF utilities
    "read_gff_annotations",
    "write_gff_header",
    "write_gff_feature",
]