#hexpSim

Example of a Python script to produce simulations of HEX-P observations, for
both the Low Energy and High Energy Telescopes (LET and HET). It uses the
reflection model RELXILL, but it can be easily adapted to other scenarios.

It requires:

- XSPEC (specifically, pyxspec must be working)
- RELXILL: http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/index.html
- HEX-P responses for LET and HET: https://hexp.org/simulation-files/

It produces simulations for a range of black hole spin using relativistic
reflection models. Prints the list of recovered spins and errors.
**WARNING**: This code has not been fully tested. It's meant to provide a
starting point, so use it at your own risk!

Usage: hexpSim.py [options]

Run a series of HEX-P simulations for BH spin measurements with RELXILL


Options:
  -h, --help            show this help message and exit
  -v, --version         show version number
  -e EXPOSURE, --exposure=EXPOSURE
                        source exposure in ks (default: 10)
  -s SCALE, --scale=SCALE
                        scale factor for background exposure w.r.t. source
                        (default: 100)
  -f FLUX, --flux=FLUX  source flux (2-10 keV) in mCrab (default: 10)
  -i ITER, --ier=ITER   Number of simulations (default: 10)

