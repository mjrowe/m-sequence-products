import numpy as np
from scipy.ndimage import correlate1d
from scipy.signal import max_len_seq

class MSequence:
    def __init__(self, MSeq: np.ndarray):
        self.MSeq = MSeq
        self.length = len(MSeq)
        self.base_peak = self._find_peak(correlate1d(MSeq, MSeq, mode='wrap'))

    def find_shift(self, shifted_array: np.ndarray):
        """ Finds the offset between a shifted sequence and the base sequence for the class. """
        derived_impulse = correlate1d(self.MSeq, shifted_array, mode='wrap')
        derived_peak = self._find_peak(derived_impulse)

        return (self.base_peak - derived_peak)%self.length

    def find_Zech_logs(self):
        self.Zech_logs = np.array([self.find_shift(self.MSeq*self._shift(self.MSeq, i)) 
                                        for i in range(1, self.length)])

    @staticmethod
    def _shift(arr: np.ndarray, shift: int):
        return np.roll(arr, shift)

    @staticmethod
    def _find_peak(impulse: np.ndarray):
        return np.argmax(abs(impulse))


class MSequenceFromFile(MSequence):
    def __init__(self, path: str):
        super().__init__(self.load(path))

    def load(self, path: str):
        """ Ensure M-sequence is {-1,1} not {0,1}. Note for nbits >= 32 there will
        be an overflow error due to the int32. """
        file_data = np.loadtxt(path, dtype=np.int32)

        return 2*file_data - 1 if 0 in np.unique(file_data) else file_data


class MSequenceFromScipy(MSequence):
    def __init__(self, order: int, state: str, taps: str):
        super().__init__(self.load(order, state, taps))

    def load(self, order: int, state: str, taps: str):
        """ Gets scipy to generate an M-sequence with the provided arguments. """
        state = list(map(int, state.split(','))) if state is not None else None
        taps = list(map(int, taps.split(','))) if taps is not None else None

        m = 2*max_len_seq(order, state=state, taps=taps)[0]-1

        return m.astype(np.int32, copy=False)

# Is typing correct