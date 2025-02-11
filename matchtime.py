#!/usr/bin/env python3
'''
Prints the time table it would take for the number of trials t at rate r
to get at least 1 match with n characters of base b and a probability p with the formula

    t = log(1-p)/log(1-1/b^n)

See https://github.com/cathugger/mkp224o/issues/27#issuecomment-568291087
'''
import math
import datetime

R = 18_000_000
P = [50, 95, 99]
N = [4, 5, 6, 7, 8]
B = 64

def t(p, n):
    return math.log(1 - p/100) / math.log(1 - 1/B**n)

print('Time to get n-symbol match at rate {R} trials per second with probability of:')
print('n', end='')
for p in P:
    print(f' {p:19d}%', end='')
print()

for n in N:
    print(n, end='')
    for p in P:
        v = t(p, n) / R
        s = str(datetime.timedelta(seconds=round(v)))
        print(f' {s:>20}', end='')
    print()
