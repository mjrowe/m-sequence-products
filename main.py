import argparse
from m_sequences import MSequenceFromFile, MSequenceFromScipy


def main():
    """ Entry point for the program """
    args = parse_command_line_input()

    if args.file_path is not None:
        m_seq = MSequenceFromFile(args.file_path)
    elif args.nbits is None:
        raise Exception('Either a file or an integer must be provided as arguments.')
    else:
        m_seq = MSequenceFromScipy(args.nbits, args.state, args.taps)

    m_seq.find_product_offsets()
    print(m_seq.product_offsets)

def parse_command_line_input():
    """ Parse the input. Either a text file containing the M-sequence must be provided,
     or inputs to interface with scipy.signal.max_len_seq(). """

    parser = argparse.ArgumentParser(argument_default=None)
    parser.add_argument('file_path', nargs='?', type=str, 
                        help='Path of text file containing the M-sequence (optional)')
    parser.add_argument('-n', '--nbits', type=int, 
                        help='The order of the M-sequence (required if no path provided, overridden otherwise)')
    parser.add_argument('-s', '--state', type=str,
                        help='Initial state of the shift register')
    parser.add_argument('-t', '--taps', type=str,
                        help='Primitive polynomial coefficients (feedback taps)')

    return parser.parse_args()


if __name__ == "__main__":
    main()

    # Can the n argument be positional as well?