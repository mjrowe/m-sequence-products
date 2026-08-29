from main import parse_command_line_input, validate_input
import numpy as np
from scipy.ndimage import correlate1d

class Decimations:
    def __init__(self, M_sequence_length: int):
        self.L = M_sequence_length
        self.n = M_sequence_length.bit_length()

        self.find_non_trivial_decimations()
        self.find_gold_decimations()

    def find_non_trivial_decimations(self):
        """ Only integers coprime to the sequence length, 2^n - 1, will form an interesting decimation.
        These phi(2^n-1) decimations give a distinct m-sequence rather than a circular shift of the original. 
        
        However as there are only phi(2^n-1)/n unique m-sequences at a given order, we further remove this 
        multiplicity (powers of two of a primitive element is also primitive with the same polynomial). """

        candidate_integers = np.arange(1, self.L, dtype=int)

        gcd = np.gcd(candidate_integers, self.L)
        coprime_integers = candidate_integers[gcd == 1]

        primitive_element_cycles = [ [i*(2**k)%self.L for k in range(self.n)] for i in coprime_integers]
        representatives = set([min(cycle) for cycle in primitive_element_cycles])

        self.non_trivial_decimations = sorted(list(representatives))
        self.num_primitive_polynomials = len(self.non_trivial_decimations)

    def find_gold_decimations(self):
        candidate_integers = np.arange(1, self.n, dtype=int)
        gold_exponents = [i for i in candidate_integers if (self.n/np.gcd(self.n,i))%2 == 1]

        self.gold_decimations = [2**k + 1 for k in gold_exponents]

    def summarise(self):
        print("*"*80, "SUMMARY", sep="\n")
        print(f"There are {self.num_primitive_polynomials} unique m-sequences at n = {self.n}, produced by decimating the original at:")
        print("d = ", self.non_trivial_decimations)
        print(f"Of these {self.gold_decimations} are Gold decimations.\n" + "*"*80)


def sanity_check(MSeq: np.ndarray):
    auto_correlation = correlate1d(MSeq, MSeq, mode='wrap')
    unique_levels = np.unique(auto_correlation).size

    return unique_levels == 2

if __name__ == "__main__":

    args = parse_command_line_input()
    m_seq = validate_input(args)

    decimations = Decimations(m_seq.length)

    for d in decimations.non_trivial_decimations:
        print(f"DECIMATION: d = {d}")
        print(f"GOLD TYPE: {d in decimations.gold_decimations}")

        new_m_sequence = m_seq.decimate(d)
        if not sanity_check(new_m_sequence):
            raise Exception("This decimation should produce an m-sequence but doesn't!")

        cross_correlation = correlate1d(m_seq.MSeq, new_m_sequence, mode='wrap')
        unique_levels = np.unique(cross_correlation)
        print(f"{unique_levels.size} LEVEL cross-correlation with d = 1: ", unique_levels, "\n")

    decimations.summarise()


    # TODO Sort out docstrings and typing
    # Gold sequences are not working!
    # Bitbucket change gold codes to gold decimations