# **FFT**
Implementation of discrete fourier transorm, forward and backward methods using complex numbers.\
complex numbers class was implemented inside.\
written in cpp with minimal modification required to use the methods in c.


# **DEPENDENCIES**
NONE.


# **INSTALLATION**
1- clone the repository.\
2- include the header in the project.



# **USAGE**
1- create an instance of the `DFT` class.\
2- call the method `SetUpTwiddlesLookupTable` to generate the twiddles lookup table, to eleminate the most of multiplication operations in the runtime.\
3- call the function: `Forward` to transfrom from time domain to frequency domain. \
4- call the function: `Backward` to transfrom from frequency domain to time domain. \
5- call the function: `GetPowerSpectral` to get the signal spectral.
