from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import argparse as ap

def my_parser() -> ap.Namespace:
    '''
    Interacts with the user input.

    Returns
    -------
    ap.Namespace
        Has all the inputs given by the user while using the script.
    '''
    parser = ap.ArgumentParser(
        description='Obtaining protein sequence from a given nucletoide sequence'
    )

    parser.add_argument(
        '-s', '--sequence',
        required=True,
        type=str,
        help='Nucleotide sequence'
    )

    parser.add_argument(
        '-a', '--all_frames',
        action='store_true',
        required=False,
        help='If set, the six reading frames will be tested.'
    )

    return parser.parse_args()

def is_rna(seq: str) -> bool:
    '''
    Checks out if there the sequence is DNA or RNA.

    Parameters
    ----------
    seq : str
        Input nucleotide sequence.

    Returns
    -------
    bool
        True if the sequence is RNA, False if is DNA.

    Raises
    ------
    ValueError
        If there are non-canonical nucleotides, then the sequence is not compatible with the program.
    '''
    if 'T' in seq:
        return False
    if any(x in seq for x in {'Z', 'J'}):
        raise ValueError('[ERROR] Your sequence has non-canonical nucleotides and is not compatible with this program.\n')
    
    return True

def finding_orfs(nt_seq: str, all_frames: bool = False) -> dict[int, list[Seq]]:
    '''
    Converts a DNA sequence to aminoacid chains by using possible ORFs and finishing them at
    both the end of the sequence and if a stop codon is found.

    Parameters
    ----------
    nt_seq : str
        Provided sequence.

    Returns
    -------
    orfs_seq : dict[int, list[Seq]]
        A dictionary with the reading frame as key and a Seq object containing the possible chains.
    '''
    if is_rna(nt_seq := nt_seq.upper()):
        dseq = str(dna_seq := Seq(nt_seq).back_transcribe())
    else:
        dna_seq = Seq(dseq := nt_seq)

    if all_frames:
        dcomp = str(dna_comp := dna_seq.reverse_complement())

    start_codon = 'ATG'; stop_codons = {'TAA', 'TAG', 'TGA'}; orfs_seq = {}; L = len(dseq)

    # Reading frames
    for i in range(3):
        orfs_seq.setdefault(i+1, [])
        valid_start = [pos for pos in nt_search(dseq,start_codon)[1:] if pos%3 == i]
        for start_idx in valid_start:
            cand_stops = (pos for pos in range(start_idx+3, L-2, 3) if dseq[pos:pos+3] in stop_codons)
            stop_idx = next(cand_stops, None)
            orf_end_excl = (stop_idx+3) if stop_idx is not None else (L - ((L - start_idx)%3))
            orfs_seq[i+1].append(dna_seq[start_idx:orf_end_excl].translate(to_stop=True))

        if all_frames:
            orfs_seq.setdefault(-(i+1), [])
            rvalid_start = [pos for pos in nt_search(dcomp,start_codon)[1:] if pos%3 == i]
            for start_idx in rvalid_start:
                cand_stops = (pos for pos in range(start_idx+3, L-2, 3) if dcomp[pos:pos+3] in stop_codons)
                stop_idx = next(cand_stops, None)
                orf_end_excl = (stop_idx+3) if stop_idx is not None else (L - ((L - start_idx)%3))
                orfs_seq[-(i+1)].append(dna_comp[start_idx:orf_end_excl].translate(to_stop=True))

    return orfs_seq

def greatest_chain(orfs: dict[int, list[Seq]]) -> list[str]:
    '''
    Given a dictionary of ORF-aminoacid sequences this function finds the largest one.

    Parameters
    ----------
    orfs : dict[int, list[Seq]]
        Dictionary ORF-aminoacid sequences

    Returns
    -------
    best : list[str]
    '''
    best = []; max_len = 0

    for v in orfs.values():
        for e in v:
            actual_len = len(e:=str(e))
            if actual_len > max_len:
                max_len = actual_len
                best = [e]
            elif actual_len == max_len:
                best.append(e)

    return best

def main():
    args = my_parser()
    orfs_seq = finding_orfs(args.sequence, args.all_frames)
    largest_orf = greatest_chain(orfs_seq)
    print(
        f'\nThe largest(s) ORF(s) found in {args.sequence} '
        'was(were):\n'
    )
    for e in largest_orf:
        print(e)

if __name__ == '__main__':
    main()