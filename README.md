This repository contains several plotting and other convenience scripts for post-processing MCMC runs from all samplers (MH, CH, NS, PC) currently used in [MontePython](https://github.com/brinckmann/montepython_public).

IMPORTANT: You must have translated your MCMC run into a default MontePython chain via

`python /path/to/montepython_public/montepython/MontePython.py info /path/to/your/MontePython/chain/{PC, NS, CH}`

The main script is `make_all.py` and it will post-process any translated Monte Python chain fully
automatically and will also plot 1D parameter histograms, 2D contours and 2D parameter
triangle plots.

Call it like this:

`python make_all.py /path/to/your/MontePython/chain/ model_name={'arbitrary string'} type_sampler={'MH', 'MN', 'NS', 'PC', 'CH'} chain_is={'2c', '2cosmo', '1c', '1cosmo'}`

Various other plotting-related options can be set directly in the script. All other scripts work in stand-alone mode and you can find their call signature in the header of each script!
