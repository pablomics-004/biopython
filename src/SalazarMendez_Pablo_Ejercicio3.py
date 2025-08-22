import SalazarMendez_Pablo_Ejercicio2 as class1

class tRNA(class1.gene):

    def __init__(
            self, 
            name: str = '', 
            length: int = 0, 
            seq: str = '', 
            gc_con: float = 0.0, 
            organism: str = '', 
            lecture_frame: int = 1, 
            anticodon: str = 'CAU',
            loops: bool = True,
            discrim_base: str = 'A',
            acceptor_stem: bool = True
    ):
        super().__init__(name, length, seq, gc_con, organism, lecture_frame)
        self.anticodon = anticodon.upper().replace('T', 'U')
        self.loops = loops
        self.discrim_base = discrim_base
        self.acceptor_stem = acceptor_stem
        self.seq = seq.upper().replace('T', 'U')
    
    def paired_codon(self):
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
        for aa in codon_aa:
            if comp_codon in codon_aa[aa]:
                return aa
            
        raise ValueError(f"\n[ERROR] No aminoacid could match {self.anticodon} anticodon\n")