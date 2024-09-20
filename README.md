# Estimator complexity for the $\mathbb{F}_{4}$ OLEAGE PCG

Sage script to select parameters for the $\mathbb{F}_4$ OLEAGE Pseudorandom Correlation Generator (PCG).
See the [paper](https://eprint.iacr.org/2024/429.pdf) for details.

We also provide a prototype [implementation](https://github.com/sachaservan/FOLEAGE-PCG) of the PCG in C.

## How to run the estimator

 1. Open a sage shell
 2. Load the module:
 
 ```sage
 sage: load("folding.sage")
 ```
 It will automatically load the `ISD.sage` module.
 3. To estimate the cost of a specific instance, run `expect_cost`. Example:
 ```sage
 sage: expect_cost(c=5, t=14, s=16, q=4, verbose=True)
 
 [...]
 
Computing the probability that the folded error has weight 14
Prange --> 181.08545120625692
Lee-Brickell --> 169.5462959218013
Stern --> 175.78460971203532
Optimized_Stern --> 164.2579898159922
MMT --> 158.96462039020085
Folded code has length 3645
Average running time of best ISD at each weight is 158.96462039020085
Weight which minimizes ISD/Proba is 62 which happens with proba 8.722688081689128e-05.
	--> Repeat 11465 (~=2^14) times ISD with aborts (each of them costing 2^142.6586575704365) ---> 156.14352524421696
There are 75628919722004322604209288760~=2^95.93292466855976  possible subgroups of size 3^10.
156.14352524421696
```
4. To get a suitable number of errors to achieve a specific security level, run `find_t`. Example:

``` sage
sage: t=find_t(c=5, s=15, q=4, security_parameter=128, verbose=True)
Starting from t=11.

[...]

t=11 --> 119.29335826621663 bits. Too small. Testing t=12.

Prange --> 158.1763289977271
Lee-Brickell --> 146.6532818400539
Stern --> 153.47915023359465
Optimized_Stern --> 144.02676671966543
MMT --> 138.09212474075946
Folded code has length 1215
Average running time of best ISD at each weight is 138.09212474075946
Weight which minimizes ISD/Proba is 44 which happens with proba 2.4516886148369733e-08.
	--> Repeat 40788215 (~=2^26) times ISD with aborts (each of them costing 2^103.55170809187008) ---> 128.83335709476145
There are 1279025522911365763892449~=2^80.08131933073459  possible subgroups of size 3^10.
For s=15, c=5, we need t=12 for a security estimated to 128.83335709476145 bits
```


## Levels of aggressivity in the parameters selection

Historically, the analysis of attacks against the decoding problem at the GV bound (essentially ISD algorithms), have been focused on optimising the exponent of the complexity. In particular, polynomial factors have often been ignored when designing parameters in code-based cryptography. On the other hand, we consider this to be overly conservative, and we count them by default. If you want to ignore the polynomial factors, you can set the global variable `conservative=True` in `ISD.sage`. Moreover, in our instances of the decoding problem, the errors is actually split into blocks having the same number of erroneous coordinates. No algorithm have been known to take this into account in the range of our code-rates. However, if you want to be even more conservative, you can set the variable `csplit=True` in `ISD.sage` to give a provable lower-bound on putative algorithms taking this structure into account.



## How to add a new ISD algorithm

It suffices to implement the complexity in the `ISD.sage` module,
following the same template as the other algorithms already there.
Once it is implemented, you can add it to the dictionnary `ISDs` at
the end. This is the one that will be exposed in `Folding.sage`. There
is no need to change the main file.
