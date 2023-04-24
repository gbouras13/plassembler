import plass_class

class Plass:
    """Plassembler Output Class"""

    def __init__(
        self,
        contig_count: int = 1,
        successful_unicycler_recovery: bool = True,
        
        ) -> None:
        """
        Parameters
        --------
        contig_count: int, required
            the number of contigs assembled by flye assembly
        successful_unicycler_recovery: bool, required
            flag determining whether unicycler finished. If false, suggests no plasmid
        """
        pass