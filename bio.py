"""
# Name: Reyhaneh Parnian
# Course Name: CSCI 1103
# Section: Project 4
# Date: 12/11/2024
"""

class DNA:
    """
    This class makes and validates DNAs.
    Attribute(s): sequence (A string of nucleotides (A, C, G, T))
    Methods: validate, __str__
    """

    def __init__(self, sequence:str) -> None:
        """
        parameter(s): sequence
        returns: None
        """
        self.sequence = sequence.upper() #storing the sequence in uppercase.

    def validate(self) -> bool:
        """
        Validate if the sequence contains only valid nucleotides (A, C, G, T).
        Returns True if the sequence is valid, False otherwise.
        """
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        for char in self.sequence:
            if char not in valid_nucleotides:
                return False
        return True

    def __str__(self) -> str:
        """
        Return the string representation of the DNA object.
        Returns The nucleotide sequence in uppercase.
        """
        return self.sequence #same thing, but is now all uppercase
    def convert_to_rna(self):
        """
        Converts current DNA into RNA.
        Returns RNA
        """
        self.sequence = self.sequence.replace("T" , "U")
        return RNA(self.sequence)
    def reverse_complement(self) -> str:
        """
        Creates the reverse complement of a given DNA.
        Returns the reverse version of DNA(str)
        """
        list_sequence = list(self.sequence)
        new_dna = []
        for num in range(0,len(self.sequence)):
            if list_sequence[num] == "A":
                list_sequence[num] = "T"
            elif list_sequence[num] == "T":
                list_sequence[num] = "A"
            elif list_sequence[num] == "C":
                list_sequence[num] = "G"
            else:
                list_sequence[num] = "C"
            new_dna.append(list_sequence[num])
        self.sequence = "".join(new_dna)
        return self.sequence[::-1]
    


CODON_DICT = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
              'UGC': 'C', 'UGU': 'C',
              'GAC': 'D', 'GAU': 'D',
              'GAA': 'E', 'GAG': 'E',
              'UUC': 'F', 'UUU': 'F',
              'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
              'CAC': 'H', 'CAU': 'H',
              'AUA': 'I', 'AUC': 'I', 'AUU': 'I',
              'AAA': 'K', 'AAG': 'K',
              'UUA': 'L', 'UUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
              'AUG': 'M',
              'AAC': 'N', 'AAU': 'N',
              'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
              'CAA': 'Q', 'CAG': 'Q',
              'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
              'AGC': 'S', 'AGU': 'S', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
              'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
              'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
              'UGG': 'W',
              'UAC': 'Y', 'UAU': 'Y'}

STOP_CODONS = {'UAG', 'UAA', 'UGA'}  # Explicit list of stop codons


class RNA:
    """
    This class makes and validates RNAs.
    Attribute(s): sequence (A string of nucleotides (A, C, G, U))
    Methods: validate, __str__
    """

    def __init__(self, sequnece:str):
        self.sequence = sequnece.upper()

    def validate(self) -> bool:
        """
        Validate if the sequence contains only valid nucleotides (A, C, G, U).
        Returns True if the sequence is valid, False otherwise.
        """
        valid_nucleotides = {'A', 'C', 'G', 'U'}
        for char in self.sequence:
            if char not in valid_nucleotides:
                return False
        return True

    def __str__(self) -> str:
        """
        Return the string representation of the RNA object.
        Returns The nucleotide sequence in uppercase.
        """
        return self.sequence #same thing, but is now all uppercase
    def to_protein(self, frame: int)->str:
        """
        Translate the RNA string into a protein sequence.
        Argument(s): frame (The reading frame to start translation (0, 1, or 2))
        Returns the resulting protein sequence.
        """
        protein = []
        start_index = frame
        # Iterate through the RNA sequence in chunks of 3
        while start_index + 3 <= len(self.sequence):
            codon = self.sequence[start_index:start_index + 3]
            if codon in STOP_CODONS:  # Check for stop codon
                break
            if codon in CODON_DICT:  # Translate codon to amino acid
                protein.append(CODON_DICT[codon])
            start_index += 3

        return ''.join(protein)
        

def read_dna_file(filename:str):
    """
    Reads the file, constructs DNA objects, converts to RNA, and outputs proteins.
    Argument(s):filename (the path to the DNA sequence file)
    """
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()

                # Skip blank lines
                if not line:
                    continue

                try:
                    frame, sequence = line.split(maxsplit=1)
                    frame = int(frame)

                    # Create DNA object
                    dna = DNA(sequence)
                    if not dna.validate():
                        print(f"Invalid sequence: {sequence}")
                        continue

                    # Convert DNA to RNA and then to protein
                    rna = dna.convert_to_rna()
                    protein = rna.to_protein(frame)
                    print(protein)
                except ValueError:
                    print(f"Invalid sequence: {line}")
    except FileNotFoundError:
        print("File not found.")



    



def demo() -> None:
    read_dna_file("DNASequences.txt")
if __name__ == "__main__":
    demo()
