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
    orfs = []
    in_orf = False
    start_loc = None
    
    # Iterate through codons in the sequence using the `codons` function
    for i, codon in enumerate(codons(sequence)):
        pos = i * 3  # Calculate the starting position of this codon in the sequence
        
        if codon in start_codons and not in_orf:
            # Start of a new ORF
            in_orf = True
            start_loc = pos

        elif codon in stop_codons and in_orf:
            # End of the ORF, add it to list
            orfs.append((start_loc, pos + 3))
            in_orf = False  # Reset the flag for a new ORF

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

    orfs = []  # Initialize orfs as an empty list to store found ORFs
    
    # Find ORFs in the original sequence (three reading frames)
    for frame in range(3):
        frame_sequence = sequence[frame:]  # Frame-specific sequence
        orfs += [(1, start + frame, stop + frame) for start, stop in find_orfs(frame_sequence, start_codons, stop_codons)]

    # Find ORFs in the reverse complement sequence (three reading frames)
    rev_sequence = sequence.reverse_complement()
    for frame in range(3):
        frame_sequence = rev_sequence[frame:]  # Frame-specific sequence
        found_orfs = find_orfs(frame_sequence, start_codons, stop_codons)
        # Adjust positions for reverse complement
        for start_loc, stop_loc in found_orfs:
            # Reverse complement strand is -1, adjust start and stop locations
            orfs.append((-1, len(sequence) - stop_loc, len(sequence) - start_loc))

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
    codon_table = {
        'AUG': 'M',  # Start codon
        'UAA': '',   # Stop codon
        'UAG': '',   # Stop codon
        'UGA': '',   # Stop codon
        'UUU': 'F',  'UUC': 'F',
        'UUA': 'L',  'UUG': 'L',
        'CUU': 'L',  'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I',  'AUC': 'I', 'AUA': 'I',
        'GUU': 'V',  'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'CCU': 'P',  'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T',  'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A',  'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'CAU': 'H',  'CAC': 'H',
        'CAA': 'Q',  'CAG': 'Q',
        'AAU': 'N',  'AAC': 'N',
        'AAA': 'K',  'AAG': 'K',
        'GAU': 'D',  'GAC': 'D',
        'GAA': 'E',  'GAG': 'E',
        'UAU': 'Y',  'UAC': 'Y',
        'UAA': '',   # Stop codon
        'UAG': '',   # Stop codon
        'CAU': 'H',  'CAC': 'H',
        'AAU': 'N',  'AAC': 'N',
        'GAU': 'D',  'GAC': 'D',
        'UGU': 'C',  'UGC': 'C',
        'UGA': '',   # Stop codon
        'CUG': 'L',
        'UGA': '',   # Stop codon
    }

    protein_seq = []
    
    # Convert the sequence to uppercase and ensure it's a multiple of 3
    seq = seq.upper()
    
    # Iterate through the sequence in steps of 3
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        # Translate the codon into an amino acid
        amino_acid = codon_table.get(codon, '')
        
        # If it's a stop codon, break the translation
        if amino_acid == '':
            break
        
        protein_seq.append(amino_acid)

    return ''.join(protein_seq)

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
    orfs = []

    # For positive strand
    for frame in range(3):
        pos = frame
        while pos < len(sequence) - 3:
            codon = sequence[pos:pos+3]
            if codon in start_codons:
                # Found a start codon, now search for stop codon
                for end in range(pos+3, len(sequence) - 2, 3):
                    stop_codon = sequence[end:end+3]
                    if stop_codon in stop_codons:
                        orfs.append((1, pos, end+3))  # Record ORF and move on
                        break  # Exit inner loop, go to next start codon
            pos += 3
    
    # Repeat for reverse complement if needed
    rev_sequence = sequence.reverse_complement()
    for frame in range(3):
        pos = frame
        while pos < len(rev_sequence) - 3:
            codon = rev_sequence[pos:pos+3]
            if codon in start_codons:
                for end in range(pos+3, len(rev_sequence) - 2, 3):
                    stop_codon = rev_sequence[end:end+3]
                    if stop_codon in stop_codons:
                        orfs.append((-1, len(sequence) - (end+3), len(sequence) - pos))
                        break
            pos += 3

    return orfs
