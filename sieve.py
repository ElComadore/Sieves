import numpy as np


def simpleSieve(N):
    P = []
    for i in range(0, N + 1):
        if i % 2 == 1:
            P.append(1)
        elif i == 2:
            P.append(1)
        else:
            P.append(0)
    P[1] = 0

    m = 3
    n = m * m

    while n <= N:
        if P[m] == 1:
            while n <= N:
                P[n] = 0
                n = n + 2 * m
        m = m + 2
        n = m * m
    return P


def simpleSegSieve(N, delta, M):
    S = []
    for i in range(0, delta + 1):
        S.append(1)
    for i in range(0, 2 - N):
        S[i] = 0

    P = simpleSieve(M)

    for m in range(1, M + 1):
        if P[m] == 1:
            c = m * np.ceil(N / m)
            NPrime = int(np.maximum(c, 2 * m))
            while NPrime <= N + delta:
                S[NPrime - N] = 0
                NPrime = NPrime + m
    return S


def subSegSieve(N, delta, M):
    S = []
    for i in range(0, delta + 1):
        S.append(1)
    for i in range(0, 2 - N):
        S[i] = 0

    deltaPrime = int(np.floor(np.sqrt(M)))
    MPrime = 1

    while MPrime <= M:
        P = simpleSegSieve(MPrime, deltaPrime, int(np.floor(np.sqrt(MPrime + deltaPrime))))
        for p in range(MPrime, int(np.minimum(M, MPrime + deltaPrime))):
            if P[p - MPrime] == 1:
                NPrime = int(np.maximum(p * np.ceil(N / p), 2 * p))
                while NPrime <= N + delta:
                    S[NPrime - N] = 0
                    NPrime = NPrime + p
        MPrime = MPrime + deltaPrime + 1
    return S


def segSieve(n, delta):
    return subSegSieve(n, delta, int(np.floor(np.sqrt(n + delta))))


n1 = 100
d1 = 1500
m1 = 15

asdf = segSieve(n1, d1)

primes = []


for j in range(0, len(asdf)):
    if asdf[j] == 1:
        primes.append(j + n1)
print(primes)