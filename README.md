# Estimator complexity for F4-OLE

## How to run the estimator

 1. Open a sage shell
 2. Load the module:
 
 ```sage
 sage: load("Folding.sage")
 ```
 It will automatically load the `ISD.sage` module.
 3. To estimate the cost of a specific instance, run `expect_cost`. Example:
 ```sage
 sage: expect_cost(c=5, t=14, s=16, q=4, verbose=True)
 
 [...]
 
Prange --> 169.34343621773516
Lee-Brickell --> 157.80428093327953
Stern --> 149.60824546567432
Optimized_Stern --> 152.51597482747042
MMT --> 140.96761009080828
Folded code has length 3645
Average running time of best ISD at each weight is 140.96761009080828
Weight which minimizes ISD/Proba is 62 which happens with proba 8.722688081689128e-05.
	--> Repeat 11465 (~=2^14) times ISD with aborts (each of them costing 2^124.60204541920282) ---> 138.08691309298325
There are 75628919722004322604209288760~=2^95.93292466855976  possible subgroups of size 3^10.
138.08691309298325
 ```
4. To get a suitable number of errors to achieve a specific security level, run `find_t`. Example:

``` sage
sage: t=find_t(c=5, s=15, q=4, security_parameter=128, verbose=True)
Starting from t=12.

[...]

t=12 --> 111.85718141458838 bits. Too small. Testing t=13.

[...]

t=13 --> 121.41917938528088 bits. Too small. Testing t=14.

[...]

Prange --> 169.34350329465923
Lee-Brickell --> 157.8043467206174
Stern --> 149.60831033537673
Optimized_Stern --> 152.51604028505164
MMT --> 140.96767327320416
Folded code has length 3645
Average running time of best ISD at each weight is 140.96767327320416
Weight which minimizes ISD/Proba is 62 which happens with proba 8.721446046640164e-05.
	--> Repeat 11466 (~=2^14) times ISD with aborts (each of them costing 2^124.60204541920282) ---> 138.08711853485795
There are 103741619611085612124067759~=2^86.42312526537337  possible subgroups of size 3^9.
For s=15, c=5, we need t=14 for a security estimated to 138.08711853485795 bits
```


## How to add a new ISD algorithm

It suffices to implement the complexity in the `ISD.sage` module,
following the same template as the other algorithms already there.
Once it is implemented, you can add it to the dictionnary `ISDs` at
the end. This is the one that will be exposed in `Folding.sage`. There
is no need to change the main file.
