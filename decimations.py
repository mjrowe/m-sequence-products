""" Decimating (sub-sampling) an m-sequence returns another m-sequence. Only integers coprime 
to the sequence length, L =  2^n - 1, will form an interesting decimation. These phi(2^n-1) 
decimations give a distinct m-sequence rather than a circular shift of the original. """

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
        """ Finds coprime decimations. As there are only phi(2^n-1)/n unique m-sequences at a 
        given order, we further remove the multiplicity arising from powers of two. That is,
        only one representative from the cycle i*{1,2,4, ..., 2^{n-1}} is required; the square 
        of a primitive element is also primitive with the same polynomial. """

        candidate_integers = np.arange(1, self.L, dtype=int)

        gcd = np.gcd(candidate_integers, self.L)
        coprime_integers = candidate_integers[gcd == 1]

        primitive_element_cycles = [self._cycle(i) for i in coprime_integers]
        representatives = set([min(cycle) for cycle in primitive_element_cycles])

        self.non_trivial_decimations = sorted(list(representatives))
        self.num_primitive_polynomials = len(self.non_trivial_decimations)

    def find_gold_decimations(self):
        """ Particular decimations that give nice three-valued cross-correlations with the 
        original sequence. """

        candidate_integers = np.arange(1, self.n, dtype=int)
        all_gold_decimations = [2**i + 1 for i in candidate_integers if (self.n/np.gcd(self.n,i))%2 == 1]

        gold_representatives = set([min(self._cycle(i)) for i in all_gold_decimations])

        self.gold_decimations = sorted(list(gold_representatives))

    def _cycle(self, i: int):
        return [i*(2**k)%self.L for k in range(self.n)]

    def summarise(self):
        print("*"*80, "SUMMARY", sep="\n")
        print(f"There are {self.num_primitive_polynomials} unique m-sequences at n = {self.n}, ",
              "produced by decimating the original at:", sep="")
        print("d = ", self.non_trivial_decimations)
        print(f"Of these {self.gold_decimations} are Gold decimations.\n" + "*"*80)


def sanity_check(MSeq: np.ndarray):
    """ M-sequences have two-level auto-correlations. """

    auto_correlation = correlate1d(MSeq, MSeq, mode='wrap')
    unique_levels = np.unique(auto_correlation).size

    return unique_levels == 2


if __name__ == "__main__":

    args = parse_command_line_input()
    m_seq = validate_input(args)

    decimations = Decimations(m_seq.length)

    for d in decimations.non_trivial_decimations:
        print(f"DECIMATION: d = {d}", 
              f"GOLD TYPE: {d in decimations.gold_decimations}", sep="\n")

        new_m_sequence = m_seq.decimate(d)
        if not sanity_check(new_m_sequence):
            raise Exception("This decimation should produce an m-sequence but doesn't!")

        cross_correlation = correlate1d(m_seq.MSeq, new_m_sequence, mode='wrap')
        unique_levels = np.unique(cross_correlation)
        print(f"{unique_levels.size} LEVEL cross-correlation with d = 1: ", unique_levels, "\n")

    decimations.summarise()