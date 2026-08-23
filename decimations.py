from m_sequences import MSequenceFromFile, MSequenceFromScipy
import numpy as np
from scipy.ndimage import correlate1d

#TODO

n = 7

m_seq = MSequenceFromScipy(n, None, None)
print(m_seq.MSeq)
m_seq.find_Zech_logs()
print(m_seq.Zech_logs)

############################

great_common_divisors = np.gcd(range(1,m_seq.length), m_seq.length)

""" Only coprime integers will form an m-sequence when decimated, phi(2^n-1) of these. """
integers_coprime_to_length = np.array(range(1,m_seq.length))[great_common_divisors == 1]

""" However there are only phi(2^n-1)/n unique m-sequences, powers of 2 correspond to the same sequence. """
primitive_elements_cycles = [[i*(2**k)%m_seq.length for k in range(n)] for i in integers_coprime_to_length]
unique_decimations = set([min(cycle) for cycle in primitive_elements_cycles])

print(unique_decimations)

for d in unique_decimations:
    new_m_sequence = m_seq.decimate(d)

    auto_correlation = correlate1d(new_m_sequence, new_m_sequence, mode='wrap')
    cross_correlation = correlate1d(m_seq.MSeq, new_m_sequence, mode='wrap')

    print(f"Decimation d = {d}")
    print("Autocorrelation: ", np.unique(auto_correlation))
    print("Cross-correlation with original: ", np.unique(cross_correlation))
    print(f"{np.unique(cross_correlation).size} LEVEL\n")