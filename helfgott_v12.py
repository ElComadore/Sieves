from bitarray import bitarray
from bitarray.util import zeros
from math import sqrt, floor


# ______Evaluation functions______#
def getNums(S, n=0):
    from numpy import array
    return array(S.search(bitarray('1'))) + n


def countNums(S):
    return S.count()


def equals(A, B):
    from bitarray.util import count_xor
    return not count_xor(A, B)


def compare(A, B, n, printNums=False):
    from bitarray.util import count_xor
    diff_count = count_xor(A, B)
    if diff_count == 0:
        print("The lists are equal")
        return
    print("The lists differ on ", diff_count, " points")
    diff = A ^ B
    A_not_B = A & diff
    B_not_A = B & diff
    if A_not_B.any():
        print("A has ", A_not_B.count(), " points which are not in B")
        if printNums:
            print(getNums(A_not_B, n))
            print('')
    if B_not_A.any():
        print("B has ", B_not_A.count(), " points which are not in A")
        if printNums:
            print(getNums(B_not_A, n))
            print('')
    return


# ______Sieve functions______#
def RemoveNonTwins(S):
    for i in S.itersearch(bitarray('100')):
        S[i:i + 3] = False
    return


def DiophAppr(alpha, Q):
    b = floor(alpha)
    p, q = b, 1
    p_m, q_m = 1, 0
    s = 1
    while q <= Q:
        if alpha == b:
            return p, -s * q_m, q
        alpha = 1 / (alpha - b)
        b = floor(alpha)
        p_p, q_p = b * p + p_m, b * q + q_m
        p_m, q_m = p, q
        p, q = p_p, q_p
        s = -s
    return p_m, s * q, q_m


def SimpleSiev(N):
    P = zeros(N + 1)  # Create bitarray of zeros representing numbers 0 to N
    P[2] = True  # Set 2 to True
    P[3::2] = True  # Set odd numbers >1 to True

    for num in range(3, int(sqrt(N)) + 1, 2):
        if P[num]:
            P[num * num::2 * num] = False
    return P


def SimpleSegSiev(n, Delta, M):
    S = zeros(Delta + 1)
    S.setall(True)
    if n <= 2:
        S[0:2 - n] = False  # Set 0 1 to false if present
    Primes = SimpleSiev(M)  # Get list of primes up to M
    for p in Primes.itersearch(bitarray('1')):
        S[max(-(n // -p), 2) * p - n::p] = False  # Set multiples of p to false
    return S


def SubSegSiev(n, Delta, M):
    S = zeros(Delta + 1)
    S.setall(True)
    S[0:2 - n] = False  # Set numbers 0 and 1 to false if present

    D_s = int(sqrt(M))  # Size of each segment (-1)
    for M_s in range(1, M + 1, D_s + 1):  # Iterate through segments
        Primes = SimpleSegSiev(M_s, D_s, int(sqrt(M_s + D_s)))  # Get primes in segment
        # Primes = SimpleSegSiev(M_s, min(M-M_s, D_s), int(sqrt((M_s+D_s))  # Get primes in segment
        for p in Primes.itersearch(bitarray('1')):
            p = p + M_s
            S[max(-(n // -p), 2) * p - n::p] = False
    return S


def SegSieveList(n, Delta):
    M = int(sqrt(n + Delta))
    S = zeros(Delta + 1)
    S.setall(True)
    if n <= 2:
        S[0:2 - n] = False
    D_s = int(sqrt(M))
    L = readInPrimes("primes1BILL", M+1)
    for M_s in range(1, M + 1, D_s + 1):
        for p in L[M_s:M_s + D_s].itersearch(bitarray('1')):
            p = p + M_s
            S[max(-(n // -p), 2) * p - n::p] = False
    return S


def NewSegSiev(n, Delta, K):
    if n ** (1 / 3) > Delta or Delta > n:
        raise Exception("Delta should be between n^(1/3) and n.")
    if K < 2.5:
        raise Exception("K must be at least 2.5")

    upper = n + Delta
    lower = n - Delta

    M = floor(K * Delta) + 1
    S = SubSegSiev(lower, 2 * Delta, M - 1)
    bound1 = sqrt(upper)
    f = sqrt(Delta / (4 * n))

    while M <= bound1:
        R = floor(M * f)
        m0 = M + R
        alpha0 = n / m0 % 1
        alpha1 = -n / (m0 * m0) % 1
        a, ainv, q = DiophAppr(alpha1, 2 * R)
        c = floor(alpha0 * q + 0.5)
        k = floor(5 * Delta * q / (4 * M))  # int-type needed for range. floor(x/y) is faster than int(x//y)
        bound2 = M + 2 * R + 1
        for j in range(-k - 1, k + 2):
            r0 = -ainv * (c + j) % q
            for m in range(m0 + r0 - ((m0 + r0 - M) // q) * q, bound2, q):  # iterating through elements in the intersection
                n_s = (upper // m) * m
                if n_s >= lower and n_s > m:  # The first criterion is less common than the second
                    S[n_s - lower] = False
        M = bound2
    return S


def intervalSieve(numInts, length):
    name = "intPrim"
    t = True
    for i in range(1, numInts):
        S = NewSegSiev(i*length, length, 2.5)
        if t:
            createPrimeList(S, name)
            t = not t
        else:
            S.remove(S[0])
            appendPrimeList(S, name)
        del S


def plotPrimes(S):
    import matplotlib.pyplot as plt
    import numpy as np
    if len(S) > 10000:
        intLength = len(S)//10000
    else:
        intLength = 1
    pi = []
    points = []
    pi.append(S[0:10].count())
    points.append(10)
    for i in range(10 + intLength, intLength*10000, intLength):
        pi.append(S[i:i+intLength].count() + pi[-1])
        points.append(i)
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    ax.plot(points, pi, label = "sieve")  # Plot some data on the axes.
    x = np.linspace(10, intLength*10000, 10000)
    ax.plot(x, np.divide(x, np.log(x)), label = "x/log(x)")  # Plot some data on the axes.
    ax.set_xlabel('x')  # Add an x-label to the axes.
    ax.set_ylabel('f(x)')  # Add a y-label to the axes.
    ax.set_title("Prime Counter vs. PNT")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    plt.show()
    print(pi[-1])
    wait = input("Press Enter to continue.")


def createPrimeList(primeArray, name):
    bits_per_line = 100
    with open(name, 'w') as f:
        f.write(primeArray[0:bits_per_line].to01())
        for i in range(bits_per_line, len(primeArray), bits_per_line):
            f.write('\n' + primeArray[i:i + bits_per_line].to01())


def readInPrimes(primeList = False, chars = False):
    if not primeList:
        primeList = "primes_in_bits"
    f = open(primeList, "rt")
    if not chars:
        return bitarray(f.read())
    return bitarray(f.read(chars))


def appendPrimeList(primeArray, name):
    bits_per_line = 100
    with open(name, 'a') as f:
        f.write(primeArray[0:bits_per_line].to01())
        for i in range(bits_per_line, len(primeArray), bits_per_line):
            f.write('\n' + primeArray[i:i + bits_per_line].to01())


# ______IO functions for running file as script______
def createconf():  # Creates config file
    import configparser
    config = configparser.ConfigParser(allow_no_value=True)
    config['Values'] = {'; Values for the sieve': None,
                        'n': 1000000,
                        'Delta': 10000,
                        'K': 2.5}
    config['Filenames'] = {'; Names of files to write results to': None,
                           'bitarray': 'primes_in_bits',
                           'statistics': 'stats'}
    with open(_config_file, 'w') as f:
        config.write(f)
    return


def conf():  # Handles config file
    import os.path
    import configparser
    if not os.path.isfile(_config_file):  # Checks if config file exists
        createconf()
        print('Created config.ini\n' +
              'Please edit it and then run this script again')
        quit()

    # Read config file
    config = configparser.ConfigParser()
    try:
        config.read(_config_file)
        global _bit_file, _stats_file, _args
        _bit_file = config['Filenames']['bitarray']
        _stats_file = config['Filenames']['statistics']
        _args['sieve'] = NewSegSiev.__name__
        _args['n'] = int(config['Values']['n'])
        _args['Delta'] = int(config['Values']['Delta'])
        _args['K'] = float(config['Values']['K'])
    except:
        print('config.ini not set correctly!\n' +
              'You can delete it to get a new default version')
    return


def info():  # Creates string with info about function call
    info = ''
    for key, value in _args.items():
        info += str(key) + ' = ' + str(value) + '\n'
    return info


def writeout(pr, S):  # Writes stats to one file and bitarray to another
    import sys
    import pstats
    bits_per_line = 100  # This can be changed to any positive integer
    num_line = 'Found ' + str(countNums(S)) + ' primes!\n'
    cprof = _stats_file + '.cprof'
    pr.print_stats()
    pr.dump_stats(cprof)  # Writes stats but not human readable 
    with open(_stats_file, "w") as f:  # Converts stats to human readable and writes to stats file
        ps = pstats.Stats(cprof, stream=f)
        ps.sort_stats('time')
        ps.print_stats()
        f.write('Found ' + str(countNums(S)) + ' primes!\n\n' + info())
    with open(_bit_file, "w") as f:  # Writes bitarray to bitarray file
        f.write(S[0:bits_per_line].to01())
        for i in range(bits_per_line, len(S), bits_per_line):
            f.write('\n' + S[i:i + bits_per_line].to01())
    return


def main():  # Main function that runs when file is run as script
    import cProfile
    conf()  # Handle config
    print(info())  # Print info about the sieving to be done
    pr = cProfile.Profile()
    pr.enable()
    S = bitarray()
    plotPrimes(readInPrimes("primes1BILL"))
    # intervalSieve(11, 50000000)
    # S = readInPrimes("intPrim")  # Sieve
    # S = SegSieveList(0, 2*_args['n'])
    pr.disable()
    print(len(S))
    print('Found ' + str(countNums(S)) + ' primes!\ncProfile:')  # Print n.o. primes found

    # print(getNums(S, _args['n'] - _args['Delta']))
    writeout(pr, S)  # Write to files

    return


if __name__ == "__main__":  # If file is run as script
    # Globals
    _config_file = 'config.ini'
    _args = {}
    _bit_file = ''
    _stats_file = ''

    main()  # Run
