import SalazarMendez_Pablo_Ejercicio2 as class1

class tRNA(class1.gene):

    def __init__(
            self, 
            name: str = '', 
            length: int = 0, 
            seq: str = '', 
            gc_con: float = 0.0, 
            lecture_frame: int = 1, 
            anticodon: str = 'CAU',
            loops: bool = True,
            discrim_base: str = 'A',
            acceptor_stem: bool = True
    ):
        super().__init__(name, length, seq, gc_con, lecture_frame)
        self.anticodon = anticodon.upper().replace('T', 'U')
        self.loops = loops
        self.discrim_base = discrim_base
        self.acceptor_stem = acceptor_stem
        self.seq = seq.upper().replace('T', 'U')
    
    def paired_codon(self) -> str:
        comp = {
            'A' : 'U',
            'U' : 'A',
            'C' : 'G',
            'G' : 'C'
        }
        return ''.join(comp.get(self.anticodon[i], self.anticodon[i]) for i in range(len(self.anticodon)-1,-1,-1))
    
    def noncanonical_bases(self) -> bool:
        noncanon = {chr(x) for x in range(66,91)} - {'U','C','G'}
        return True if any(c in noncanon for c in self.anticodon) else False

    def codon_translator(self) -> str:
        if (l := len(self.anticodon)) > 3 or l < 3:
            raise ValueError(f"\n[ERROR] Anticodon must have 3 bases length, it has {l}\n")
        
        if self.noncanonical_bases():
            raise ValueError(f'\n[ERROR] The anticodon sequence {self.anticodon} contains noncanonical bases\n')
        
        comp_codon = self.paired_codon()
        codon_aa = {
            'Met': {'AUG'},
            'Phe': {'UUU', 'UUC'},
            'Leu': {'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'},
            'Ile': {'AUU', 'AUC', 'AUA'},
            'Val': {'GUU', 'GUC', 'GUA', 'GUG'},
            'Ser': {'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'},
            'Pro': {'CCU', 'CCC', 'CCA', 'CCG'},
            'Thr': {'ACU', 'ACC', 'ACA', 'ACG'},
            'Ala': {'GCU', 'GCC', 'GCA', 'GCG'},
            'Tyr': {'UAU', 'UAC'},
            'His': {'CAU', 'CAC'},
            'Gln': {'CAA', 'CAG'},
            'Asn': {'AAU', 'AAC'},
            'Lys': {'AAA', 'AAG'},
            'Asp': {'GAU', 'GAC'},
            'Glu': {'GAA', 'GAG'},
            'Cys': {'UGU', 'UGC'},
            'Trp': {'UGG'},
            'Arg': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},
            'Gly': {'GGU', 'GGC', 'GGA', 'GGG'},
            'STOP': {'UAA', 'UAG', 'UGA'}
        }
        for aa, codons in codon_aa.items():
            if comp_codon in codons:
                return aa
            
        raise ValueError(f"\n[ERROR] No aminoacid could match {self.anticodon} anticodon\n")
    
class ncRNA():
    def __init__(
            self, 
            seq: str = '', 
            gc_con: float = 0, 
            sec_struct: str = 'hairpin', 
            cellular_loc: str = 'cytosol',
            seed_seq: str = None,
            complex: str = '',
            circular: bool = False
    ):
        self.length = len(seq)
        self.gc_con = gc_con
        self.seq = seq.upper().replace('T', 'U')
        self.second_struct = sec_struct.lower()
        self.cellular_loc = cellular_loc.lower()
        self.circular = circular
        self.seed_seq = seed_seq
        self.associated_complex = complex

    def ncRNA_type(self) -> str:
        if self.circular:
            return 'circRNA'
        elif self.length > 200:
            return 'lncRNA'
        elif self.seed_seq:
            return 'miRNA'
        elif 'PIWI' in self.complex:
            return 'piRNA'
        else:
            raise ValueError(f"\n[ERROR] No match for the given ncRNA\n")
        
    def seed_sequence(self, ncrna_type: str = ''):
        if ncrna_type != 'miRNA':
            return None
        return self.seq[1:8]
    
    def gc_percentage(self) -> None:
        self.gc_con += (
            (total_gc := self.seq.count('G') + self.seq.count('C')) / (total_gc + self.seq.count('A') + self.seq.count('T'))
        ) * 100