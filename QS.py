#!/usr/bin/env python

import sys
import random


'''
This first section of code is an implementation of Dixon's factorization
method using Gaussian Elimination, this is not optimal, but is simpler than
using Block Lanczos by a lot.
'''

# https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
def gauss(M):
  marks = [False] * len(M)
  for j in xrange(len(M[0])):
    for i in xrange(len(M)):
      if M[i][j] == 1:
        marks[i] = True
        for k in xrange(j):
          if M[i][k] == 1:
            for row in xrange(len(M)):
              M[row][k] = (M[row][k] + M[row][j]) % 2
        for k in xrange(j+1, len(M[0])):
          if M[i][k] == 1:
            for row in xrange(len(M)):
              M[row][k] = (M[row][k] + M[row][j]) % 2
        break
  return (marks, M)

# A dependent column is any column with value "1" that exists
# within an unmarked row
def get_dep_cols(row):
  ret = []
  for i in xrange(len(row)):
    if row[i] == 1:
      ret.append(i)
  return ret

# Sum a set of rows in GF(2) 
def row_add(new_row, current):
  ret = current
  for i in xrange(len(M[new_row])):
    ret[i] ^= M[new_row][i]
  return ret

def is_dependent(cols, row):
  for i in cols:
    if row[i] == 1:
      return True
  return False

# Finds a set of linear dependencies which include the row specified
def find_linear_deps(row):
  ret = []
  dep_cols = get_dep_cols(M[row])
  current_rows = [row]
  current_sum = M[row][:]
  for i in xrange(len(M)):
    if i == row:
      continue
    if is_dependent(dep_cols, M[i]):
      current_rows.append(i)
      current_sum = row_add(i, current_sum)
      if sum(current_sum) == 0:
        ret.append(current_rows[:])
  return ret

def testdep(dep):
  x = y = 1
  for row in dep:
    x *= smooth_vals[row][0]
    y *= smooth_vals[row][1]
  return xgcd(x - isqrt(y), N)[0]

'''
This next section of code implements the single polynomial quadratic sieve
for finding smooth numbers about some factor base.
'''

# Euler's totient function given a prime p
def phi(p):
  return p-1

# The Legendre symbol (a/p) determines quadratic residues mod p
def legendre(a, p):
  if a % p == 0:
    return 0
  return pow(a, (p - 1) // 2, p)

# Miller-Rabin probabilistic prime test
def miller(n, trials=5):
  s = 0
  d = n - 1
  if n == 2:
    return True
  if n % 2 == 0:
    return False
  while True:
    quotient, remainder = divmod(d, 2)
    if remainder == 1:
      break
    s += 1
    d = quotient
  for _ in xrange(trials):
    a = random.randrange(2, n)
    composite = True
    if pow(a, d, n) == 1:
      continue
    for i in xrange(s):
      if pow(a, 2**i * d, n) == n-1:
        composite = False 
        break
    if composite:
      return False
  return True

# Newton's method for finding sqrt(N)
def isqrt(n):
  x = n
  y = (x + 1) // 2
  while y < x:
    x = y
    y = (x + n // x) // 2
  return x

# Extended euclidian algorithm
def xgcd(a,b):
  prevx, x = 1, 0
  prevy, y = 0, 1
  while b:
    q, r = divmod(a,b)
    x, prevx = prevx - q*x, x  
    y, prevy = prevy - q*y, y
    a, b = b, r
  return a, prevx, prevy

# https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
def tonelli(n, p):
  q = p - 1
  s = 0
  z = 2
  i = 1
  while q % 2 == 0:
    q //= 2
    s += 1
  if s == 1:
    return pow(n, (p + 1) // 4, p)
  while z < p and p - 1 != legendre(z, p):
    z += 1
  c = pow(z, q, p)
  r = pow(n, (q + 1) // 2, p)
  t = pow(n, q, p)
  m = s
  t2 = 0
  while (t - 1) % p != 0:
    t2 = (t * t) % p
    i = 1
    while i < m:
      if (t2 - 1) % p == 0:
        break
      t2 = (t2 * t2) % p
      i += 1
    b = pow(c, 1 << (m - i - 1), p)
    r = (r * b) % p
    c = (b * b) % p
    t = (t * c) % p
    m = i
  return r

# Create a factor base up to some bound B in which all primes have
# legendre symbol 1, ie. n is a quadratic residues mod p
def create_base(n, B):
  base = []
  i = 2
  while len(base) < B:
    if legendre(n, i) == 1:
      base.append(i)
    i += 1
    while not miller(i):
      i += 1
  return base

# Sieving polynomial function of the form:
#   (Ax + B)^2 - N
def poly(x, a, b, n):
  return ((a * x + b) ** 2) - n

# We are solving for X in this equation specifically: 
#   (Ax + B)^2 - N = 0 mod p.
#
# Note that if a != 1, then gcd(a,p) for all p must be 1 since we will have
# to find the multiplicative modular inverse of a mod p.
def solve(a, b, n):
  start_vals = []
  for p in base:
    ainv = 1
    if a != 1:
      g, ainv, _ = xgcd(a, p)
      assert g == 1
    r1 = tonelli(n, p)
    r2 = (-1 * r1) % p
    start1 = (ainv * (r1 - b)) % p
    start2 = (ainv * (r2 - b)) % p
    start_vals.append([start1, start2])
  return start_vals

'''
This final section of code contains some utility functions and the program
entry point along with some initialization and setup
'''

# Use trial division to produce an exponent vector for n with respect
# to the factor base. Note exponents are in GF(2)
def trial(n, base):
  ret = [0] * len(base)
  if n > 0:
    for i in xrange(len(base)):
      while n % base[i] == 0:
        n //= base[i]
        ret[i] = (ret[i] + 1) % 2
  return ret


N = int(sys.argv[1])
a = 1
b = isqrt(N) + 1
bound = 50
base = create_base(N, bound)
needed = phi(base[-1]) + 1

sieve_start = 0
sieve_stop = 0
sieve_interval = 100000

M = []
smooth_vals = []
start_vals = solve(a, b, N)
seen = set([])

print '=='
print 'Searching for %d integers which are smooth to the factor base' % needed
print '==\n'

# Sieve until we have enough smooth numbers
while len(smooth_vals) < needed:
  sieve_start = sieve_stop
  sieve_stop += sieve_interval
  interval = [poly(x, a, b, N) for x in xrange(sieve_start, sieve_stop)]

  for p in range(len(base)):
    t = start_vals[p][0]

    while start_vals[p][0] < sieve_start + sieve_interval:
      while interval[start_vals[p][0] - sieve_start] % base[p] == 0:
        interval[start_vals[p][0] - sieve_start] /= base[p]
      start_vals[p][0] += base[p]

    # Sieve the values using both start positions if they aren't equal
    if start_vals[p][1] != t:
      while start_vals[p][1] < sieve_start + sieve_interval:
        while interval[start_vals[p][1] - sieve_start] % base[p] == 0:
          interval[start_vals[p][1] - sieve_start] /= base[p]
        start_vals[p][1] += base[p]

  # To form the congruence of squares, we will need the full radical value 
  # which is a*x + b. Note we also need unique relations.
  for i in xrange(sieve_interval):
    if interval[i] == 1:
      x = sieve_start + i
      y = poly(x, a, b, N)
      exp = trial(y, base)
      if not tuple(exp) in seen:
        smooth_vals.append(((a * x) + b, y))
        M.append(exp)
        seen |= set([tuple(exp)])

  print 'After sieving the next interval we have %d' % (len(smooth_vals))

print '\n=='
print 'Running Gaussian Elimination'
print '==\n'

marks, M = gauss(M)
print 'Done!'

print '\n=='
print 'Testing Linear Dependencies'
print '==\n'

for i in range(len(marks)):
  if not marks[i]:
    deps = find_linear_deps(i)
    for dep in deps:
      print dep
      gcd = testdep(dep)
      if gcd != 1 and gcd != N:
        print '\nNon-trivial factor: %d' % gcd
        sys.exit(0)
