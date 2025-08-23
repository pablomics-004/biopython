import SalazarMendez_Pablo_Ejercicio2 as class1

class tRNA(class1.gene):

    CODON_TO_AA = {
        'AUG': 'Met',
        'UUU': 'Phe', 'UUC': 'Phe',
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr',
        'CAU': 'His', 'CAC': 'His',
        'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn',
        'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp',
        'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys',
        'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
        'UAA': 'stop', 'UAG': 'stop', 'UGA': 'stop'
    }

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
        super().__init__(name=name, length=length, seq=seq, gc_con=gc_con,
                 organism='', lecture_frame=lecture_frame)

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
        return ''.join(comp.get(b, b) for b in self.anticodon[::-1])
    
    def noncanonical_bases(self) -> bool:
        allowed = {'A','T','C','G'}
        return True if any(c not in allowed for c in self.anticodon) else False

    def codon_translator(self) -> str:
        if (l := len(self.anticodon)) > 3 or l < 3:
            raise ValueError(f"\n[ERROR] Anticodon must have 3 bases length, it has {l}\n")
        
        if self.noncanonical_bases():
            raise ValueError(f'\n[ERROR] The anticodon sequence {self.anticodon} contains noncanonical bases\n')
        
        comp_codon = self.paired_codon()
        aa = self.CODON_TO_AA.get(comp_codon, None)
            
        if not aa:
            raise ValueError(f"\n[ERROR] No aminoacid could match {self.anticodon} anticodon\n")
        return aa
    
class Protein(tRNA):
    def __init__(self, name = '', length = 0, seq = '', gc_con = 0, lecture_frame = 1, anticodon = 'CAU', loops = True, discrim_base = 'A', acceptor_stem = True):
        super().__init__(name, length, seq, gc_con, lecture_frame, anticodon, loops, discrim_base, acceptor_stem)
        self.aminoacids = ''

    @staticmethod
    def has_nocanonical(bases: str) -> bool:
        allowed = {'A','T','C','G'}
        return True if any(b not in allowed for b in bases) else False

    @staticmethod
    def anticodon_method(codon: str) -> str:
        comp = {
            'A' : 'U',
            'U' : 'A',
            'C' : 'G',
            'G' : 'C'
        }
        return ''.join(comp.get(b, b) for b in reversed(codon))
    
    def codon_translator(self, anticodon: str) -> str: # Polymorphism
        if (l := len(anticodon)) > 3 or l < 3:
            raise ValueError(f"\n[ERROR] Anticodon must have 3 bases length, it has {l}\n")
        
        if self.has_nocanonical(anticodon):
            raise ValueError(f'\n[ERROR] The anticodon sequence {anticodon} contains noncanonical bases\n')
        
        aa = self.CODON_TO_AA.get(anticodon, None)
            
        if not aa:
            raise ValueError(f"\n[ERROR] No aminoacid could match {anticodon} anticodon\n")
        return aa
    
    def codons_2_aa(self):
        codons = self.first_frame()[::-1]
        self.aminoacids = ''.join(self.codon_translator(ac) for ac in [self.anticodon_method(c) for c in codons])
    
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
        self.gc_con = (
            (total_gc := self.seq.count('G') + self.seq.count('C')) / (total_gc + self.seq.count('A') + self.seq.count('T'))
        ) * 100