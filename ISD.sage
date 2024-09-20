"""
List of ISD algorithms.

Includes:
 - Prange
 - Lee-Brickell
 - Stern-Dumer over Fq
 - Stern-Dumer over Fq with improvements from [P09]
 - MMT over Fq


The level of aggressivity in the parameters selection is controlled by
the variables conservative and csplit.
"""

conservative=False # If true, then ignore some polynomial factors to be as conservative as possible.
csplit=False # Applies a penalty for being csplit

def Tgauss(n, k):
    """
    Lower bound on the cost Gaussian elimination.
    Favors the attacker.
    """
    return n*(n-k)

def c_split_loss(n, t, c):
    """
    Penalty for having the error regularly split into c blocks
    """
    if csplit:
        return float(log(binomial(n/c, t/c)^c/binomial(n, t), 2))
    else:
        return 0

def Prange(t, k, n, full=False, **kwargs):
    """
    Bit complexity of classic Prange algorithm.

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
    """
    linalg = log(Tgauss(n,k), 2)
    res = float(log(binomial(n, t) / binomial(n - k, t), 2) + linalg)
    if full:
        return res, {}
    else:
        return res


def Lee_brickell(t, k, n, q, p, verbose=False, **kwargs):
    """
    Bit complexity of classic Lee&Brickell algorithm.

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
     - q: size of the field
     - p: Number of errors in the information set
    """
    if verbose:
        print(f"n={n}, k={k}, t={t}, p={p}")

    P_succ = (binomial(k, p)*binomial(n-k, t-p))/binomial(n, t)
    L = binomial(k, p)*(q-1)^p
    T_gauss = Tgauss(n, k)

    C_iter = T_gauss + L
    res = C_iter / P_succ

    return float(log(res, 2))

def Lee_brickell_opt(t, k, n, q, p_min=0, full=False, verbose=False, **kwargs):
    """
    Look for the optimal p in Lee&Brickell algorithm.
    In practice, p=2 is optimal.

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
     - q: size of the field
     - full: Return tuple (BitComplexity, p_optimal)
     - verbose: More verbose
    """
    T_min = Lee_brickell(t, k, n, q, p_min)
    for p in range(p_min+1, t):
        T_cur = Lee_brickell(t, k, n, q, p)
        if T_cur < T_min:
            T_min=T_cur
            p_min = p
    if verbose:
        print(f"Lee-Brickell for t={t}, n={n}, k={k} is optimal for p={p_min} --> {T_min}")
    if full:
        return T_min, {'p_min': p_min}
    else:
        return T_min

def Stern(t, k, n, q, p, ell, verbose=False):
    """
    Bit complexity of q-ary Stern Algorithm.

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
     - q: size of the field
     - p: Number of errors in the information set
     - ell:
    """
    if t == 0:
        return 0

    T_gauss = Tgauss(n, k+ell) # n*(n-k-ell)
    L_1 = binomial(floor((k+ell)/2), floor(p/2)) * (q-1)^(p/2)
    L = (L_1^2)/q^ell


    if conservative:
        C_iter = L*(n-k-ell)*(k+ell) + max(L_1, L)+L_1
    else:
        C_iter = L*(n-k-ell)*(k+ell) + 2*ell*(k+ell)*L_1 + L*ell

    P_succ = (binomial(n-k-ell, t-p)*binomial(floor((k+ell)/2), floor(p/2))^2)/binomial(n, t)
    if P_succ == 0:
        return oo

    res = ((T_gauss + C_iter)*log(q, 2)) / P_succ
    return float(log(res, 2))

def Stern_opt(t, k, n, q=4, p_min=0, ell_min=0, full=False, verbose=False, **kwargs):
    """
    Find the optimal parameters and bit complexity for Stern Algorithm.

    Parameters::
      - t: number of errors
      - k: dimension of the code
      - n: length of the code
      - q: size of the field
      - full: Return tuple (BitComplexity, p_optimal, ell_optimal)
      - verbose: More verbose
    """
    p_min, ell_min = 0, 0
    T_min = Stern(t, k, n, q, p_min, ell_min, verbose=verbose)
    if t < 30:
        p_max, ell_max = 8, 30
    elif t < 80:
        p_max, ell_max = 15, 50
    else:
        p_max, ell_max = 15, 100
    for p in range(p_min+1, min(t+1, p_max)):
        for ell in range(ell_min+1, min(n-k+1, ell_max)): # l can in theory be as big as n-k + p - t but this is good only for large t (say linear in n)
            T_cur = Stern(t, k, n, q, p, ell, verbose=verbose)
            if T_cur < T_min:
                T_min = T_cur
                p_min = p
                ell_min = ell
    Stern(t, k, n, q, p_min, ell_min, verbose=verbose)
    if verbose:
        print(f"Complexity Stern for decoding {t} errors, in an [{n}, {k}]_{q}-code is optimal for p={p_min} l={ell_min} ---> {T_min}")
    if full:
        return T_min, {'p_min': p_min, 'ell_min': ell_min}
    else:
        return T_min


def Optimized_Stern(t, k, n, q, p, ell, verbose=False):
    """
    Bit complexity of q-ary Stern Algorithm with [P09] optimizations

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
     - q: size of the field
     - p: Number of errors in the information set
     - ell:
    """
    if t == 0:
        return 0
    T_gauss = Tgauss(n, k+ell) # n*(n-k-ell)
    L_1 = binomial(floor(k/2), p) * (q-1)^p
    L_2 = binomial(k-floor(k/2), p)*(q-1)^p
    T_lists = (k/2 - p + 1 + L_1 + L_2) * ell
    N_cols = (L_1 * L_2) / (q^ell)
    T_checks = ((q/(q - 1)) * (t - 2*p + 1) * 2*p * (1 + (q-2)/(q-1)) * N_cols)

    Psucc = (L_1 * L_2 * binomial(n-k-ell, t-2*p)/(binomial(n, t)*(q-1)^(2*p)))
    if Psucc < 1e-53:
        return oo

    C_iter = (T_gauss + T_lists + T_checks) * log(q, 2)

    res = C_iter / Psucc
    return float(log(res, 2))

def Optimized_Stern_opt(t, k, n, q=4, p_min=0, ell_min=0, full=False, verbose=False, **kwargs):
    """
    Find the optimal parameters and bit complexity for Stern Algorithm.

    Parameters::
      - t: number of errors
      - k: dimension of the code
      - n: length of the code
      - q: size of the field
      - full: Return tuple (BitComplexity, p_optimal, ell_optimal)
      - verbose: More verbose
    """
    if t < 30:
        p_max, ell_max = 15, 30
    elif t < 80:
        p_max, ell_max = 15, 50
    else:
        p_max, ell_max = 15, 100
    p_min, ell_min = 0, 0
    T_min = Optimized_Stern(t, k, n, q, p_min, ell_min, verbose=verbose)
    for p in range(p_min+1, min(floor(t / 2), floor(k / 2), p_max)):
        for ell in range(ell_min+1, min(n - k, n - k + 2 * p - t, ell_max)): # l can in theory be as big as n-k + 2p - t but this is good only for large t (say linear in n)
            T_cur = Optimized_Stern(t, k, n, q, p, ell, verbose=verbose)
            if T_cur < T_min:
                T_min = T_cur
                p_min = p
                ell_min = ell
    Optimized_Stern(t, k, n, q, p_min, ell_min, verbose=verbose)
    if verbose:
        print(f"Complexity Optimized_Stern for decoding {t} errors, in an [{n}, {k}]_{q}-code is optimal for p={p_min} l={ell_min} ---> {T_min}")
    if full:
        return T_min, {'p_min': p_min, 'ell_min': ell_min}
    else:
        return T_min

def MMT(t, k, n, q, p, ell, verbose=False, **kwargs):
    """
    Bit complexity of q-ary MMT algorithm

    Parameters:
     - t: Number of errors
     - k: Dimension of the code
     - n: Length of the code
     - q: size of the field
     - p: Number of errors in the information set
     - ell:
    """
    if t == 0:
        return 0
    T_gauss = Tgauss(n, k+ell)
    P_succ1 = (binomial(k+ell, p)*binomial(n-k-ell, t-p))/binomial(n, t)
    P_succ2 = binomial(floor((k+ell)/2), floor(p/4))^4/binomial(k+ell, floor(p/2))^2

    P_succ = P_succ1*P_succ2
    if P_succ==0:
        return oo

    r = log(binomial(p, floor(p/2)), q)

    L00 = binomial(floor((k+ell)/2), floor(p/4))*(q-1)^(floor(p/4))
    L02 = (binomial(floor((k+ell)/2), floor(p/4))^2*(q-1)^(floor(p/2)))/q^r
    L = (binomial(floor((k+ell)/2), floor(p/4))^4*(q-1)^p)/(q^(ell+r))

    if conservative:
        C_iter = T_gauss + L*(n-k-ell)*(k+ell) + max(L00, L02, L) + L00
    else:
        C_iter = T_gauss + L*(n-k-ell)*(k+ell) + max(L00, L02, L)*ell*3 + 4*(k+ell)*ell*L00

    res = (C_iter*log(q, 2))/P_succ

    return float(log(res, 2))


def MMT_opt(t, k, n, q=4, p_min=0, ell_min=1, full=False, verbose=False, **kwargs):
    """
    Find the optimal parameters and bit complexity for q-ary MMT Algorithm.

    Parameters::
      - t: number of errors
      - k: dimension of the code
      - n: length of the code
      - q: size of the field
      - full: Return tuple (BitComplexity, p_optimal, ell_optimal)
      - verbose: More verbose
    """
    ell_min=1
    p_min=1
    if ell_min==0:
        ell_min=1
    T_min = MMT(t, k, n, q, p_min, ell_min, verbose=verbose)
    iter=0
    if t < 30:
        p_max, ell_max = 8, 30
    elif t < 160:
        p_max, ell_max = 15, 50
    else:
        p_max, ell_max = 15, 100
    for p in range(1, min(t, p_max)+1):
        for ell in range(1, min(n-k+p-t+1, ell_max)): # l can in theory be as big as n-k + 2p - t but this is good only for large t (say linear in n)
            T_cur = MMT(t, k, n, q, p, ell, verbose=verbose)
            if T_cur < T_min:
                T_min = T_cur
                p_min = p
                ell_min = ell
            else:
                iter+=1

    if verbose:
        print(f"Complexity MMT for decoding {t} errors, in an [{n}, {k}]_{q}-code is optimal for p={p_min} l={ell_min} ---> {T_min}")
    if full:
        return T_min, {'p_min': p_min, 'ell_min': ell_min}
    else:
        return T_min

def BJMM(t, k, n, q, p, ell, eps, verbose=False, **kwargs):
    T_gauss = Tgauss(n, k+ell) # n*(n-k-ell)

    p1 = p/2+eps

    P_succ1 = (binomial(k+ell, p)*binomial(n-k-ell, t-p))/binomial(n, t)
    P_succ2 = binomial(floor((k+ell)/2), floor(p1/2))^4/binomial(k+ell, floor(p1))^2
    P_succ = P_succ1*P_succ2
    if P_succ == 0:
        return oo


    R=0
    for eps1 in range(eps+1):
        R+=binomial(p, floor(p/2)-eps1)*binomial(floor(p/2)+eps1, 2*eps1)*binomial(k+ell-p, eps-eps1)*(q-1)^(2*eps1+(eps-eps1))
    r = log(R, q)

    L00 = binomial(floor((k+ell)/2), floor(p1/2))*(q-1)^(p1/2)
    L02 = (binomial(floor((k+ell)/2), floor(p1/2))^2*(q-1)^(p1))/q^r
    L = (binomial(floor((k+ell)/2), floor(p1/2))^4*(q-1)^(2*p1))/(q^(ell+r))


    if conservative:
        C_iter = T_gauss + L*(n-k-ell)*(k+ell) + max(L00, L02, L) + L00
    else:
        C_iter = T_gauss + L*(n-k-ell)*(k+ell) + max(L00, L02, L)*ell*3 + L00*ell*(k+ell)*4

    res = (C_iter*log(q, 2))/P_succ

    return float(log(res, 2))


def BJMM_opt(t, k, n, q=4, p_min=0, ell_min=0, eps_min=0, full=False, verbose=False, **kwargs):
    """
    Find the optimal parameters and bit complexity for q-ary BJMM Algorithm.
    Never seems to improve on q-ary MMT.

    Parameters::
      - t: number of errors
      - k: dimension of the code
      - n: length of the code
      - q: size of the field
      - full: Return tuple (BitComplexity, p_optimal, ell_optimal)
      - verbose: More verbose
      - eps_min: in practice always 0.
    """
    T_min = BJMM(t, k, n, q, p_min, ell_min, eps_min, verbose=verbose)
    for p in range(p_min+1, min(t+1, 10)):
        for ell in range(ell_min+1, min(n - k, 30)): # l can in theory be as big as n-k + 2p - t but this is good only for large t (say linear in n)
            for eps in range(floor(p/2)+1):
                T_cur = BJMM(t, k, n, q, p, ell, eps, verbose=verbose)
                if T_cur < T_min:
                    T_min = T_cur
                    p_min = p
                    ell_min = ell
                    eps_min = eps
    if verbose:
        print(f"Complexity BJMM for decoding {t} errors, in an [{n}, {k}]_{q}-code is optimal for p={p_min} l={ell_min} eps={eps_min}---> {T_min}")
    if full:
        return T_min, {'p_min': p_min, 'ell_min': ell_min}
    else:
        return T_min


ISDs = {"Prange": Prange,
        "Lee-Brickell": Lee_brickell_opt,
        "Stern": Stern_opt,
        "Optimized_Stern": Optimized_Stern_opt,
        "MMT": MMT_opt,
        # "BJMM": BJMM_opt,
        }
