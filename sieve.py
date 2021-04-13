import numpy as np


def simpleSieve(N):  # Original sieve of Eratosthenes
    P = []
    for i in range(0, N + 1):
        if i % 2 == 1:  # Simple little even 'wheel'
            P.append(1)
        elif i == 2:
            P.append(1)  # Making sure two is prime
        else:
            P.append(0)  # All evens are zero
    P[1] = 0

    m = 3  # Start with our first prime
    n = m * m  # looking starting ar its square

    while n <= N:
        if P[m] == 1:
            while n <= N:
                P[n] = 0  # m prime therefore n not
                n = n + 2 * m  # look at the next odd multiple of m; n = m(m + 2)
        m = m + 2  # Next odd number
        n = m * m
    return P


def simpleSegSieve(N, delta, M):  # sieve over a little interval segment
    S = []
    for i in range(0, delta + 1):  # Setting all the number to be prime candidates
        S.append(1)
    for i in range(0, 2 - N):  # In case we have 0 or 1, set them to 0
        S[i] = 0

    P = simpleSieve(M)  # Use simple Era to get our sieving primes

    for m in range(1, M + 1):
        if P[m] == 1:
            c = m * np.ceil(N / m)  # a possible starting point for nPrime
            NPrime = int(np.maximum(c, 2 * m))
            while NPrime <= N + delta:  # looking only in the interval
                S[NPrime - N] = 0  # as nPrime is always some multiple of m prime, get rid of those positions
                NPrime = NPrime + m  # next multiple of m
    return S


def subSegSieve(N, delta, M):  # segmented sieve which calls only the primes it needs each time; "subroutine"
    S = []
    for i in range(0, delta + 1):  # Prime candidates
        S.append(1)
    for i in range(0, 2 - N):  # if 0 or 1 get rid
        S[i] = 0

    deltaPrime = int(np.floor(np.sqrt(M)))  # an optimised interval length to find primes in; m' + delta'
    MPrime = 1  # where does this prime interval start

    while MPrime <= M:
        P = simpleSegSieve(MPrime, deltaPrime, int(np.floor(np.sqrt(MPrime + deltaPrime))))  # getting them primes
        for p in range(MPrime, int(np.minimum(M, MPrime + deltaPrime))):
            if P[p - MPrime] == 1:
                NPrime = int(np.maximum(p * np.ceil(N / p), 2 * p))  # setting n' as multiple of p
                while NPrime <= N + delta:  # composites in the desired interval
                    S[NPrime - N] = 0  # get rid of them
                    NPrime = NPrime + p  # next prime multiple
        MPrime = MPrime + deltaPrime + 1
    return S


def segSieve(n, delta):  # segmented optimal(?) sieve of Era
    return subSegSieve(n, delta, int(np.floor(np.sqrt(n + delta))))  # optimal? maybe


def diophApprox(alpha, Q):  # don't fucking ask me
    b = int(np.floor(alpha))
    p = b
    q = 1
    pMinus = 1
    qMinus = 0
    s = 1

    while q <= Q:
        if alpha == b:
            return [p, -s * qMinus, q]
        alpha = 1 / (alpha - b)
        b = int(np.floor(alpha))
        pPlus = b * p + pMinus
        qPlus = b * q + qMinus

        pMinus = p
        qMinus = q

        p = pPlus
        q = qPlus

        s = -s
    return [int(pMinus), int(s * q), int(qMinus)]


def newSegSieve(n, delta, K):  # segmented sieve with diophantine approx
    SPrime = subSegSieve(n - delta, 2 * delta, K * delta)  # get an approximation of all the primes; some false
    S = []
    n0 = n - delta
    for j in range(-delta, delta + 1):  # format is rather weird, is just to reflect that the interval is n +- delta
        S.append(SPrime[j + delta])  # mirrors our earlier, perhaps faulty calc
    M = int(np.floor(K * delta)) + 1

    while M <= np.sqrt(n + delta):  # gets a bit fucked, only passes if n large and delta smallish; I think
        R = int(np.floor(M * np.sqrt(delta / (4 * n))))
        m0 = M + R

        a1 = n / (m0 * m0) % 1
        a0 = n / m0 % 1
        nu = 5 * delta / (4 * M)

        aaq = diophApprox(a1, 2 * R)
        c = int(np.floor(a0 * aaq[2] + 0.5))
        k = int(np.floor(nu * aaq[2]))

        for j in range(-k - 1, k + 2):
            r0 = -aaq[1] * (c + j) % aaq[2]

            kMin = int(np.ceil((M - m0 - r0) / aaq[2]))  # begin work for dealing with the intersection
            kMax = int(np.floor((M + 2 * R - m0 - r0) / aaq[2]))

            intersection = []

            for i in range(kMin, kMax + 1):  # construct valid elements of the intersection
                intersection.append(m0 + r0 + aaq[2] * i)
            for m in intersection:
                nPrime = int(np.floor((n + delta)) / m) * m  # similar things again with setting n' multiple of m
                if n - delta <= nPrime:
                    if nPrime <= n + delta:
                        if nPrime > m:
                            S[nPrime - n0] = 0  # if in interval and larger than m; kill it!
        M = M + 2 * R + 1  # start point of new M interval; have to make it larger than M + 2R
    return S


def subSegSieveFac(n, delta, M):  # factorise numbers in interval using subroutine of segment sieve
    F = []
    Pi = []
    for i in range(0, delta + 1):
        F.append([])  # making a list of empty lists
        Pi.append(1)  # everything has 1 as a factor
    deltaPrime = int(np.floor(np.sqrt(M)))  # sub interval length for finding primes
    MPrime = 1  # where that interval starts

    while MPrime <= M:  # only using primes less than M
        P = simpleSegSieve(MPrime, deltaPrime, int(np.floor(np.sqrt(MPrime + deltaPrime))))  # get useful primes
        for p in range(MPrime, MPrime + deltaPrime):
            if P[p - MPrime] == 1:
                k = 1  # which power of the prime?
                d = p  # start with the first power!
                while d <= n + delta:
                    nPrime = d * int(np.ceil(n / d))  # same same bb
                    while nPrime <= n + delta:
                        if k == 1:  # when the first power divides
                            F[nPrime - n].append([p, 1])  # add it to the list; [prime, power]
                        else:
                            for i in range(0, len(F[nPrime - n])):  # finds the same prime in the list with small power
                                if F[nPrime - n][i][0] == p:
                                    F[nPrime - n][i][1] = k  # if same update power

                        Pi[nPrime - n] = p * Pi[nPrime - n]  # add that prime factor to the divisor list
                        nPrime = nPrime + d  # next multiple of prime power
                    k = k + 1  # next prime power
                    d = p * d  # p^k
        MPrime = MPrime + deltaPrime  # next sub interval to take primes from
    return [F, Pi]


def segSieveFac(n, delta):  # optimal interval length once again I suppose
    [F, Pi] = subSegSieveFac(n, delta, int(np.floor(np.sqrt(n + delta))))

    for nPrime in range(n, n + delta + 1):  # now properly collecting the primes and there power in F
        if Pi[nPrime - n] != nPrime:
            p0 = int(nPrime / Pi[nPrime - n])
            F[nPrime - n].append([p0, 1])
    return F


def newSegSieveFac(n, delta, K):  # uses the diophantine thing to generate primes and factor instead
    [FPrime, PiPrime] = subSegSieveFac(n - delta, 2 * delta, K * delta)
    F = []
    Pi = []
    n0 = n - delta

    for j in range(-delta, delta + 1):
        F.append(FPrime[j + delta])
        Pi.append(PiPrime[j + delta])

    M = int(np.floor(K * delta)) + 1

    while M <= np.sqrt(n + delta):
        R = int(np.floor(M * np.sqrt(delta / (4 * n))))
        m0 = M + R

        a1 = n / (m0 * m0) % 1
        a0 = n / m0 % 1
        nu = 5 * delta / (4 * M)

        aaq = diophApprox(a1, 2 * R)
        c = int(np.floor(a0 * aaq[2] + 0.5))
        k = int(np.floor(nu * aaq[2]))

        for j in range(-k - 1, k + 2):
            r0 = -aaq[1] * (c + j) % aaq[2]

            kMin = int(np.ceil((M - m0 - r0) / aaq[2]))
            kMax = int(np.floor((M + 2 * R - m0 - r0) / aaq[2]))

            intersection = []

            for i in range(kMin, kMax + 1):
                intersection.append(m0 + r0 + aaq[2] * i)
            for m in intersection:
                nPrime = int(np.floor((n + delta)) / m) * m
                if n - delta <= nPrime:
                    if nPrime <= n + delta:  # very very similar down to here
                        if int(nPrime / Pi[nPrime - n0]) % m == 0:  # now finding more prime divisors (I think)
                            if nPrime % (m * m) == 0:
                                Pi[nPrime - n0] = m * m * Pi[nPrime - n0]
                                F[nPrime - n0].append([m, 2])
                            else:
                                Pi[nPrime - n0] = m * Pi[nPrime - n0]
                                F[nPrime - n0].append([m, 1])
        M = M + 2 * R + 1  # new interval once again

    for nPrime in range(n0, n0 + 2 * delta + 1):  # same meme again
        if Pi[nPrime - n0] != nPrime:
            p0 = int(nPrime / Pi[nPrime - n0])
            F[nPrime - n0].append([p0, 1])
    return F


# TODO Implement maybe one of those summation estimates? Who knows:)

def mobFun(facList):
    for factors in facList:
        if factors[1] > 1:
            return 0
    if facList[0][0] == 1:
        return 1
    return int(np.power(-1, len(facList) % 2))

def mobSumCheck(factorised, bound=0):
    m = 1
    ressdf = 1
    max = 1
    for facList in factorised:
        m = m + mobFun(facList)
        if m > max:
            max = m
        if (m > bound):
            print('Holy fuck')
    print(max)
    return m


n1 = 1000000
d1 = n1 - 2
m1 = 13
K1 = 2.5

# asdf = newSegSieve(n1, d1, K1)
# asdfasdf = segSieveFac(n1, d1)
asdfasdfasdf = newSegSieveFac(n1, d1, K1)

m = mobSumCheck(asdfasdfasdf, np.sqrt(n1))

primes = []

print(m)
