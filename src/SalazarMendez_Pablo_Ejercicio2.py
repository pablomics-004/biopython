class gene():
    
    def __init__(
            self, 
            name: str, 
            length: int = 0, 
            seq: str = '', 
            gc_con: float = 0.0, 
            organism: str = '', 
            lecture_frame: int = 1, 
    ):
        self.name = name
        self.length = length
        self.seq = seq.upper()
        self.gc_percentage = gc_con
        self.organism = organism
        self.frame = lecture_frame
        self.nt = ''

    def gc_content(self) -> None:
        self.gc_percentage += (
            (total_gc := self.seq.count('G') + self.seq.count('C')) / (total_gc + self.seq.count('A') + self.seq.count('T'))
        ) * 100
    
    def len_calc(self) -> int:
        self.length += len(self.seq)

    def first_frame(self):
        frame_codons = []
        for i in range(0,len(self.seq),3):
            codon = self.seq[i:i+3]
            if len(codon) == 3:
                frame_codons.append(codon)
        return frame_codons
    
    def second_frame(self):
        frame_codons = []
        for i in range(1,len(self.seq),3):
            codon = self.seq[i:i+3]
            if len(codon) == 3:
                frame_codons.append(codon)
        return frame_codons
    
    def third_frame(self):
        frame_codons = []
        for i in range(2,len(self.seq),3):
            codon = self.seq[i:i+3]
            if len(codon) == 3:
                frame_codons.append(codon)
        return frame_codons
    
    def frame(self):
        match self.frame:
            case 1:
                return self.first_frame()
            case 2:
                return self.second_frame()
            case 3:
                return self.third_frame()
            case _:
                raise ValueError('\n[ERROR] No such frame have been found.\n')
    
    def base_type(self) -> None:
        uracil = self.seq.count('U')
        if uracil > 0:
            self.nt = 'RNA'
        else:
            self.nt = 'DNA'
    
    def possible_stop_cods(self):
        self.base_type()
        match self.nt:
            case 'RNA':
                stop = {'UAG', 'UAA', 'UGA'}
            case 'DNA':
                stop = {'TAG', 'TAA', 'TGA'}
        pass

lexA = gene(
    name='lexA',
    seq='ATGAAAGCGTTAACGGCCAGGCAACAAGAGGTGTTTGATCTCATCCGTGATCACATCAGCCAGACAGGTATGCCGCCGACGCGTGCGGAAATCGCGCAGCGTTTGGGGTTCC',
    organism='Escherichia coli K-12',
    lecture_frame=2
)
lexA.gc_content()
print(f'GC Content: {lexA.gc_percentage}%, {lexA.frame} frame: {lexA.second_frame()}')