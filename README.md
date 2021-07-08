[![hpc ci](https://gitlab.com/cscs-ci/arbor-sim/arbor/badges/master/pipeline.svg)](https://gitlab.com/cscs-ci/arbor-sim/arbor/-/commits/master)
[![basic ci](https://github.com/arbor-sim/arbor/actions/workflows/basic.yml/badge.svg)](https://github.com/arbor-sim/arbor/actions/workflows/basic.yml)
[![gitpod](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/arbor-sim/arbor)
[![docs](https://readthedocs.org/projects/arbor/badge/?version=latest)](https://docs.arbor-sim.org/en/latest/)

# Arbor Library

[Arbor](https://arbor-sim.org) is a library for implementing performance portable network simulations of multi-compartment neuron models.

An installation guide and library documentation are available online at [Read the Docs](http://arbor.readthedocs.io/en/latest/).

[Submit a ticket or start a discussion](https://github.com/arbor-sim/arbor) if you have any questions or want help.

### Citing Arbor

Please refer to [our documentation](https://docs.arbor-sim.org/en/latest/index.html#citing-arbor).

#### Compile bestand Lennart
set -ex
(cd arbor; python3 setup.py install --user --neuroml --gpu=cuda)
python3 /home/lennart@dc.cube.local/max/arbor/example/noise_correlation.py

### compile bestand max
cmake .. -DARB_USE_BUNDLED_LIBS=on -DARB_WITH_PYTHON=on -DPYTHON_EXECUTABLE=$(which python3) -DARB_WITH_MPI=off
cp -r /home/max/Documents/ARBORMASTER/arbor/build/python/arbor /home/max/.local/lib/python3.8/site-packages/.
python ../example/noise_correlation.py