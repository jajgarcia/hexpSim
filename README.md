# hexpSim

Example of a Python script to produce simulations of Athena+HEX-P observations, using the reflection model RELXILL. 

Requires:

- XSPEC (specifically, pyxspec must be workin)
- RELXILL: http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/index.html
- Athena responses: http://x-ifu-resources.irap.omp.eu/PUBLIC/RESPONSES/CC_CONFIGURATION/
- HEX-P responses: https://hexp.org/simulation-files/

It produces simulations for a range of BH spin using relativistic reflection models. Prints the list of recovered spins and errors.
**WARNING**: This code has not been fully tested. It's meant to provide an starting point, so use it at your own risk!
