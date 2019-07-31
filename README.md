Extremely simple single polynomial quadratic sieve in Python
===

I set out to write an extremely simple, barebones implementation of the: (https://en.wikipedia.org/wiki/Quadratic_sieve)[Quadratic sieve]

Mainly for didactic purposes. I wanted to understand some asymptomatically 
good integer factorization algorithms and the quadratic sieve is far simpler 
than the faster alternative, the general number field sieve (GNFS).

This version of the quadratic sieve is single polynomial and uses Gaussian
Elimination, meaning that it's runtime is far from optimal. Still, it should
be better than trial divison for sufficiently larger integers.

The file `QS.py` is the implementation of the quadratic sieve, you can run it
simply by invoking:

  `python QS.py <number>`

Theoretically, it can handle numbers of arbitrary size, but the largest I tried
with this implementation was about 10^30, which is small in comparision to
what most implementations support and this took a good 30 min to factor.
