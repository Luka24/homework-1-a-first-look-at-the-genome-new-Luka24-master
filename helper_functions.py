from typing import Tuple, Generator, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[int(loc.start):int(loc.end)]

        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, loc.start.position, loc.end.position))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((region.strand, loc.start.position, loc.end.position))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (loc.start.position, loc.end.position, region.strand)
                )
                print("\t", str(ex))
                
    # Some ORFs in paramecium have lenghts not divisible by 3. Remove these
    orfs = [orf for orf in orfs if (orf[2] - orf[1]) % 3 == 0]

    return orfs

def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regions) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI-provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[int(loc.start):int(loc.end)]
        
        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, int(loc.start), int(loc.end)))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Ensure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert codon not in stop_codons, f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine; add it to the ORFs
            orfs.append((region.strand, int(loc.start), int(loc.end)))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (int(loc.start), int(loc.end), region.strand)
                )
                print("\t", str(ex))

    # Some ORFs may not be divisible by 3; remove these
    orfs = [orf for orf in orfs if (orf[2] - orf[1]) % 3 == 0]

    return orfs

def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)

    """
    in_orf = False
    end = -1
    start = -1
    orfs = []
    for i in range(0, len(sequence) - 2, 3):
        if not in_orf and sequence[i:i + 3] in start_codons:
            start = i
            in_orf = True
        elif in_orf and sequence[i:i + 3] in stop_codons:
            end = i+3
            in_orf = False
            orfs.append((start, end))
    return orfs




def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """

    orfs = []
    for offset in range(3):
        orfs.extend([(1,elem[0]+offset,elem[1]+offset) for elem in find_orfs(sequence[offset:], start_codons, stop_codons)])
    seq_comp = sequence[::-1].replace("A", "N").replace("T", "A").replace("N", "T") \
        .replace("C", "N").replace("G", "C").replace("N", "G")
    seq_len = len(sequence)
    for offset in range(3):
        orfs.extend([(-1,seq_len-elem[1]-offset,seq_len-elem[0]-offset) for elem in find_orfs(seq_comp[offset:], start_codons, stop_codons)])
    
    return orfs





def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """
    translations = {
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "C": ["TGT", "TGC"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "F": ["TTT", "TTC"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
        "H": ["CAT", "CAC"],
        "I": ["ATT", "ATC", "ATA"],
        "K": ["AAA", "AAG"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "M": ["ATG"],
        "N": ["AAT", "AAC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "Q": ["CAA", "CAG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "W": ["TGG"],
        "Y": ["TAC", "TAT"]
    }

    def get_key(element_to_find):
        """Helper function to get the amino acid for a given codon."""
        for key, values in translations.items():
            if element_to_find in values:
                return key
        return ""  # Return empty string if codon is not found

    protein_sequence = ""
    for i in range(0, len(seq) - 2, 3):  # Iterate over the sequence in steps of 3
        codon = seq[i:i+3]
        amino_acid = get_key(codon)
        if amino_acid:  # Only add valid amino acids
            protein_sequence += amino_acid
        else:
            print(f"Warning: Invalid codon {codon} at position {i} ignored")

    return protein_sequence

# Example usage
nucleotide_sequence = "AUGGCCUAA"  # Example nucleotide sequence (start, amino acids, stop)
protein_sequence = translate_to_protein(nucleotide_sequence)
print(protein_sequence)  # Output: "MA"



def find_all_orfs_nested(sequence, start_codons, stop_codons):
    """Bonus problem: Find ALL the possible ORF candidates in the sequence using
    the updated definition of ORFs.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    def find_nested(sequence, start_codons, stop_codons):
        start_found = False
        starts = []
        end = -1
        orfs = []
        for i in range(0, len(sequence) - 2, 3):
            if sequence[i:i + 3] in start_codons:
                starts.append(i)
                start_found = True
            elif start_found and sequence[i:i + 3] in stop_codons:
                end = i+3
                start_found = False
                for s in starts:
                    orfs.append((s, end))
                starts = []
        return orfs

    orfs = []
    for offset in range(3):
        orfs.extend([(1,elem[0]+offset,elem[1]+offset) for elem in find_nested(sequence[offset:], start_codons, stop_codons)])
    seq_comp = sequence[::-1].replace("A", "N").replace("T", "A").replace("N", "T") \
        .replace("C", "N").replace("G", "C").replace("N", "G")
    seq_len = len(sequence)
    for offset in range(3):
        orfs.extend([(-1,seq_len-elem[1]-offset,seq_len-elem[0]-offset) for elem in find_nested(seq_comp[offset:], start_codons, stop_codons)])
    
    return orfs
