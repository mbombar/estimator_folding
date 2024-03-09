from scipy.optimize import newton

load("ISD.sage")

A.<X> = QQ[]


def h(q, x):
    """q-ary entropy function."""
    return -x * log(x / (q - 1), q) - (1 - x) * log(1 - x, q)


def GV(R, q=4):
    """
    Compute the GV bound for a given rate: h_q^-1(1-R).

    Args:
     - R: code rate
     - q: size of the field (default is F4)
    """
    return newton(lambda x: h(q, x) - (1 - R), 0.1)

def doom_loss(s, q=4):
    """
    Speed up decoding by factor sqrt(|G|).
    """
    return float(-(s/2)*log(q-1, 2))

def prange_c_split_doom(c, t, s, q=4):
    """
    Compute a lower bound the cost of Prange algorithm in the c_split regime.
    Using theorem 3.6 from projective space paper
    Parameters:
     - c: compression factor
     - t: weight per block
     - s: log of the size of the group, i.e number of variables, |G|=(q-1)^s
    """
    N=(q-1)^s
    w=c*t
    n=c*N
    k=(c-1)*N
    return float(Prange(w, k, n) + doom_loss(s) + c_split_loss(n, w, c))


def estimator_init(q=4, c=4, s=16, target_security=128):
    """Chooses the minimal weight which ensures a good enough security
    level before folding.

    Parameters::

    - q: Size of the field.
    - c: Compression factor.
    - s: Logsize of the group. Equivalently, number of variables.
    - target_security: Minimal number of bits of security.
    """
    t=1
    while prange_c_split_doom(c, t, s, q)<target_security:
        t*=2
    l=t//2
    r=t
    while (l<=r):
        t=(l+r)//2
        if prange_c_split_doom(c, t, s, q) > target_security:
            if prange_c_split_doom(c, t-1, s) <= target_security:
                return t
            r=t-1
        else:
            l=t+1
    t=l
    return t

def min_complexity_fold(t, q=4, c=4):
    """
    Auxiliary function to recover the complexity of the best ISD
    for decoding at distance t close to GV for rate 1-1/c.

    Parametters::

    - t: Number of errors
    - q: Size of the field
    - c: Compression factor
    """
    gv=GV(1-1/c,q)
    s_fold=ceil(log(t/gv, q-1))
    N = (q-1)^s_fold
    Cmin = oo
    for ISD in ISDs.values():
        T = ISD(c*t, (c-1)*N, c*N, q=q)+doom_loss(s_fold, q)
        if Cmin > T:
            Cmin=T
    return float(Cmin)

def refine_t(q=4, c=4, s=16, t=None, target_security=128):
    """
    Takes ISDs in the folded code into account to give a first estimation of t.
    """
    if t is None:
        t=estimator_init(q, c, s, target_security)
    while min_complexity_fold(t, q, c)<target_security:
        t+=1
    return t

def f_0(l, t, q=4):
    """
    Computes the weight enumerator of the [l, l-1]-parity code, up
    to monomial X^t.

    Saves time and memory.

    f_0 = (1/q)*(((1+(q-1)*X)^l)+(q-1)*(1-X)^l)


    Parameters::

     - l: Length of the parity code.
     - t: Bound on the degree.
     - q: Size of the base field.
    """
    P = A(0)
    for k in range(t+1):
        P+=binomial(l, k)*((q-1)^k+(-1)^k*(q-1))*X^k
    return (1/q)*P


def f_1(l, t, q=4):
    """Computes the weight enumerator of the complement of the
    [l,l-1]-parity code, up to monomial X^t.

    Saves time and memory.

    f_1 = ((q-1)/q)*((1+(q-1)*X)^l-(1-X)^l)

    Parameters::

     - l: Length of the parity code.
     - t: Bound on the degree.
     - q: Size of the base field.

    """
    P = A(0)
    for k in range(t+1):
        P+=binomial(l, k)*((q-1)^k-(-1)^k)*X^k
    return ((q-1)/q)*P


def distribution_weight_fold(s, s_H, t, q=4, verbose=False):
    """Compute the probability distribution of the weight of the
    folded error.

    Outputs::

     - List [l0, ..., l_t] such that l_u is the probability that the
       folded error has weight u.


    Parameters::

    - s: Logsize of the initial group G: we have |G|=(q-1)^s.
    - s_H: Logsize of the subgroup H with respect to which the folding is defined.
    - t: Weight of the error, before folding.
    - q: Size of the field

    """
    L=[0]*(t+1)

    fold_s = s-s_H

    n=(q-1)^s
    np = (q-1)^fold_s # size of the folded code
    l = (q-1)^(s_H)  # size of the subgroup

    f0 = f_0(l, t, q)
    f1 = f_1(l, t, q)

    for u in range(t+1):
        if verbose:
            print(f"\tComputing the probability that the folded error has weight {u}")
        P = (f1^u)*(f0^(np-u))
        try:
            Atu = P[t]
        except:
            return P
        L[u] = float((binomial(np, u)*Atu)/(binomial(n, t)*(q-1)^t))
    return L


def distribution_weight_folded_full(s, s_H, t, c, q=4, verbose=False):
    """Compute the probability distribution of the weight of the full
    error.

    This is the c-fold convolution of the probability distribution of
    the folded.

    Outputs::

    - A polynomial Sum p_k X^k such that p_k is the probability that
      the folded error has weight k.


    Parameters::

    - s: Logsize of the initial group G: we have |G|=(q-1)^s.
    - s_H: Logsize of the subgroup H with respect to which the folding is defined.
    - t: Weight of the error, before folding.
    - q: Size of the field.
    """
    L = distribution_weight_fold(s, s_H, t, q, verbose)
    PL = A(L)
    return PL^c




def ISD_vec(c, t, fold_s, q=4, verbose=verbose):
    """
    Prepare the complexity of all ISD for all possible weights.

    Parameters::

    - c: Compression factor.
    - t: Number of errors per block before folding.
    - fold_s: log_size of the folded group
    - q: Size of the field
    """
    n=c*(q-1)^fold_s
    k=(c-1)*(q-1)^fold_s
    d = {}
    for algo in ISDs:
        if verbose:
            print(f"Computing cost of {algo} up to {c*t} errors in an [{n}, {k}]_{q} code...")
        ISD = ISDs.get(algo)
        L=[]
        params = [{"p_min": 0, "ell_min": 0}]
        for w in range(c*t+1):
            res, *params = ISD(w, k, n, q=q, verbose=verbose, full=True, **params[0])
            L.append(res+doom_loss(fold_s)+c_split_loss(n, c*t, c))
        d[algo] = L
        if verbose:
            print(f"...Done")
    return d



def expect_cost(c, t, s, q=4, ISD_l=None, verbose=False, offset=0):
    """
    Computes the expected cost of the folding attack.

    Parameters:

    - c: Compression factor
    - t: Number of errors per block
    - s: Logsize of the group
    - q: Size of the field
    - ISD_l: Complexity of decoding at all c*t possible distance. Saves time if given (default None).
    - verbose: More verbose
    - offset: Sometimes it is a little bit better to fold less. This is the offset from the closest target folding from GV.
    """

    R = 1-1/c
    gv = GV(R, q)

    fold_s = ceil(log((t)/gv, q-1))+offset
    s_H = s-fold_s
    np=(q-1)^fold_s
    l=(q-1)^(s_H)

    Cost_folding = c*(q-1)^s

    if ISD_l is None:
        ISD_l = ISD_vec(c, t, fold_s, q, verbose=verbose)

    ISD_min = list(map(min, zip(*list(ISD_l.values())))) # Best Decoding complexity for each weight

    P = distribution_weight_folded_full(s, s_H, t, c, q, verbose)

    Tmin = sum([C*p for (C, p) in zip(ISD_min, P)])

    C_dict = {}
    for algo in ISD_l:
        C_dict[algo] = sum([C*p for (C, p) in zip(ISD_l[algo], P)])
    Tmin2 = min(C_dict.values())

    ratioComplexityProba = [(2^a+Cost_folding)/b if b!=0 else oo for (a, b) in zip(ISD_min, P)]
    ratio_min = min(ratioComplexityProba)
    w_min = ratioComplexityProba.index(ratio_min)
    Tmin3 = float(log(ratio_min, 2))

    if verbose:
        for algo in C_dict:
            print(f"{algo} --> {C_dict[algo]}")
        print(f"Folded code has length {c*(q-1)^fold_s}")
        print(f"Average running time of best ISD at each weight is {Tmin}")
        print(f"Weight which minimizes ISD/Proba is {w_min} which "\
              f"happens with proba {float(P[w_min])}.\n\t--> Repeat "\
              f"{ceil(1/P[w_min])} (~=2^{ceil(log(1/P[w_min], 2))}) times "\
              "ISD with aborts (each of them costing "\
              f"2^{float(log(2^ISD_min[w_min]+Cost_folding, 2))}) ---> {Tmin3}")
        print(f"There are {gaussian_binomial(s, s_H, q-1)}~=2^{float(log(gaussian_binomial(s, s_H, q-1),2))}  possible subgroups of size {q-1}^{s_H}.")
    return min(Tmin, Tmin2, Tmin3)


def prepare_ISD_list(c, t, s, q, ISD_l=None, verbose=False):
    """
    Compute the cost of all enabled ISD for all possible weights on the folded code.

    Outputs::

    A dictionary Algo --> List of the complexity of decoding with Algo for all weights from 0 to c*t

    Parameters::

    - c: Compression factor
    - t: Number of errors per block
    - s: Logsize of the group
    - q: Size of the field
    - ISD_l: If already precomputed, simply update the list of complexities
    - verbose: More verbose
    """
    R=1-1/c
    gv=GV(R, q)
    fold_s = ceil(log(t/gv, q-1)) # folding as close to GV bound as possible
    n=c*(q-1)^fold_s
    k=(c-1)*(q-1)^fold_s

    if ISD_l is None:
        ISD_l=ISD_vec(c, t, fold_s, q, verbose=verbose)

    for ISD in ISDs:
        algo = ISDs.get(ISD)
        params = [{"p_min": 0, "ell_min": 0}]
        ISD_complexity = ISD_l.get(ISD, [])
        for w in range(len(ISD_complexity), c*t+1):
            res, *params = algo(w, k, n, q=q, verbose=verbose, full=True, **params[0])
            ISD_complexity.append(res+doom_loss(fold_s)+c_split_loss(n, c*t, c))
        ISD_l[ISD] = ISD_complexity
    return ISD_l
    

def find_t(c, s, t=None, security_parameter=128, q=4, verbose=False):
    """This is the core of this script. It finds the required number
    of errors per block to reach the targeted security_level.

    Parameters::

    - c: Compression factor
    - s: Logsize of the group
    - t: A lower bound on the result. Saves time if provided and close to the result.
    - security_parameter: Targeted number of bits of security
    - q: Size of the field
    - verbose: More verbose
    """
    if t is None:
        t = refine_t(q, c, s)
    if verbose:
        print(f"Starting from t={t}.")
    ISD_l = prepare_ISD_list(c, t, s, q, ISD_l=None, verbose=verbose)
    cost=expect_cost(c, t, s, q, ISD_l, verbose) # cost of folding
    while cost < security_parameter:
        if verbose:
            print(f"\nt={t} --> {cost} bits. Too small. Testing t={t+1}.\n\n")
        t+=1
        ISD_l = prepare_ISD_list(c, t, s, q, ISD_l, verbose=verbose)
        cost=expect_cost(c, t, s, q, ISD_l, verbose=verbose)
    print(f"For s={s}, c={c}, we need t={t} for a security estimated to {cost} bits")
    return t
