# m-sequence-products

The purpose of this repository is threefold:
1. The Zech logarithms of $\text{GF}(2^n)$ are calculated via products of maximum length (m-) sequences.
2. Gold sequences are created from m-sequence decimations.
3. We give expository notes on how m-sequences arise from finite field (Galois) theory.

## Background

The document *'M-Sequences from Galois Fields'* covers (3) and introduces just enough field theory to understand the trace definition of an order $n$ m-sequence
```math
m[i] = (-1)^{\text{Tr}(\alpha^i)}
```
in terms of some primitive element $\alpha \in \text{GF}(2^n)$. It covers the alternative definitions, polynomials over $`\{0,1\}`$ and linear feedback shift registers, which is often not made clear in existing literature.

Taking the product of an m-sequence with a cyclically shifted version of itself gives the same m-sequence with a curious offset.
```math
m[i] m[i+\tau] = m[i+Z_\alpha(\tau)]
```
These offsets are the Zech logarithms of $\text{GF}(2^n)$ corresponding to the primitive element $\alpha$. This is relevant in the modelling of non-linear systems as a Volterra series
```math
y[n] = \sum_i h_1[i] x[n-i] + \sum_{ij} h_2[i,j] x[n-i] x[n-j] + \cdots
```
which naturally contains sucessive products of the input $x[n]$. If the 'kernels' $h_i = 0$ for $i > 1$ (i.e. linearity) then $h_1$ is the system impulse response. When using m-sequences to extract $h_1$ through cross-correlation, $h_1 = y \star x$ (with $x=m$), any non-linear behaviour will cause unwanted spikes to appear at positions given by the Zech logarithms [1].

## Code
The program `main.py` calculates products of a user supplied m-sequence with all possible shifted versions of itself (i.e. all $\tau$). As the circular autocorrelation of an m-sequence is an impulse (see the notes for a beautiful proof!) one can find the resulting offset $Z_\alpha(\tau)$ by the delay of this peak. The program thus calculates the Zech logarithms of $\text{GF}(2^n)$ for the primitive element that corresponds to the m-sequence. This is not a particularly efficient method for doing this calculation but rather a pedagogical approach.

A second program `decimations.py` generates different m-sequences formed from decimating (sub-sampling) the original sequence. For special (co-prime) decimations taking products with the original sequence generates what is known as Gold sequences. These are the program's output.

## Usage
Given a text file of an m-sequence, one integer per line, the Zech logarithms are produced via
```{sh}
$ python3 main.py mseq.txt
```
If the sequence is written as $`\{0,1\}`$ it will automatically be converted to $`\{-1,1\}`$ (note the order). Alternatively one can use `scipy.signal.max_len_seq()` to generate the sequence from its polynomial taps (see the `scipy` documentation). Use the `-n` (order) `-t` (taps) and `-s` (seed) flags, for example
```{sh}
$ python3 main.py -n 3 -t 1 -s 1,1,1
```
gives the logarithms [5,3,2,6,1,4] of $\text{GF}(8)$ corresponding to the primitive polynomial $X^3 + X^2 + 1$. Note for binary fields the Zech logarithms satisfy $Z_\alpha(\tau) = k \Leftrightarrow Z_\alpha(k) = \tau$.

Note if a text file is provided any flags will be ignored.

## Requirements
Any recent version of `numpy` and `scipy` should work, I used versions 1.26.4 and 1.11.2 respectively.

## References

[1] Matthew Wright, *Comments on 'Aspects of MLS Measuring Systems'*, Journal of the Audio Engineering Society, 43, 1 (1995).

No artificial intelligence tools were used in the creation of this repository, all mistakes are entirely mine!
