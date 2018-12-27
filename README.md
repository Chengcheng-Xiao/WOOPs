
```
 _    _  _____  ___________
| |  | ||  _  ||  _  | ___ \
| |  | || | | || | | | |_/ /__
| |/\| || | | || | | |  __/ __|
\  /\  /\ \_/ /\ \_/ / |  \__ \
 \/  \/  \___/  \___/\_|  |___/
```
# Wannier Orbital Overlap Population tools (WOOPs)

[THIS PACKAGE IS STILL UNDER DEVELOPMENT]

A post-processing tool written in python to get Wannier Orbital Overlap Population (WOOP), Wannier Orbital Position Population (WOPP)* from [Wannier90](https://github.com/wannier-developers/wannier90) package.


Before getting into things, you might want to check out these two papers:
1.  [Phys. Rev. B 56 12847 (1997)](http://dx.doi.org/10.1039/b918890h)
2.  [Phys. Rev. B 91, 195120 (2015)](http://dx.doi.org/10.1103/PhysRevB.91.195120)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

For the script to work, you need to have an valid installation of `python` (2.7.x or 3.x both work).
Also, `numpy` and `h5py` package are needed, you can install them by pip:
```
pip install h5py numpy
```
or by conda
```
conda install h5py numpy
```
if you use a supercomputer and don't have enough privilege:

1. install anaconda by downloading from [here](https://www.anaconda.com/download/) and upload it to your directory.
2. using queue system to install anaconda by `chmod 755 anaconda*.sh && ./anaconda*.sh`
3. install `numpy` and `h5py` by download them from [here](https://anaconda.org/anaconda/numpy) and [here](https://anaconda.org/anaconda/h5py), upload them as well.
4. manually install package by `conda install numpy*.tar.bz2` and `conda install *h5py*.tar.bz2`

### Installing

Python script does not need manual compilation, so the installation is very easy.

1. Download the script:
```
wget https://github.com/Chengcheng-Xiao/WOOPs/blob/master/WOOPs.py
```

2. give correct permission:
```
chmod 755 WOOPs.py
```

3. You can also put it into your system's PATH (assume you are using `BASH`):
```
echo export PATH=$( pwd )/WOOPs.py:'$PATH' >> .bashrc
```

## Useage
File need for WOOPs:

1. `wannier90_u_AO.mat` (Unitary matrix of atomic orbitals from valence band and conduction interpolation)
2. `wannier90_u_MO.mat` (Unitary matrix of molecular orbitals from only valence band interpolation)
3. `wannier90_r.dat` (Use master branch for compatibility)
4. `input.woops`

Detailed preparation for Wannier90 generated files can be found in Wanneir90's user guid or in the `example` folder.

The format of input.woops follows:
```
readmode   = XX               #[must be True or False]
cal        = XX               #[Things you want to do]
num_mo     = XX               #[number of MO]
num_ao     = XX               #[number of AO]
num_kpts   = XX               #[number of kpts]
cell_param = XX XX XX         #[cell_param]
cell_dim   = XX               #[cell_dim: 0D, 1D, 2D or 3D] NOTE: currently 1D->x, 2D->xy
cprec      = XX             #controls the PRINTING precision of C_matrix, lower means accurate.   default=1e-4
bprec      = XX             #controls the PRINTING precision of WOOP_matrix, lower means accurate default=1e-4
```
*IMPORTANT*: Every tag must be written in lowercase, full length without abbreviation.

`cal` include:

| NAME                   | REQURE                                     |
|:----------------------:|:------------------------------------------:|
| `get_orbital`          | [needs `wannier90*` file]                  |
| `check_completeness`   | [needs `get_AO` and `get_MO`]              |
| `get_c_mat`            | [needs `get_AO` and `get_MO`]              |
| `WOOP` [default]       | [needs `get_C_mat`]                        |
| `get_charge`           | [needs `WOOP`]                             |
| `WOPP` [experimental]  | [needs `WOOP` and `get_r_mat`]             |

When the input file is properly prepared:
```
./WOOPs.py
```
or if you put it in your path:
```
WOOPs.py
```
This may take a few minutes, consider submitting it to the queue system.


These commands will produce a copy of data stored in hdf5 format and a text file with corresponding name for the calculation.

## Running the tests

Go check the description in `example` folder.

## Notes
*WOPP is still under development, the problem now I'm facing is that the phase factor in Wannier90's r matrix are not "correct" (they are arbitrary). As a result, the ferroelectric polarization decomposition cannot be correctly calculated.

## TODO_list
1. Complete a simplified description of WOOP and its capability.

2. Fix phase factor error in WOPP interface.

3. Increase efficiency.

4. Add progress bar.

## How to cite

For the method please cite the following paper in any publications arising from the use of this code:

  J. Bhattacharjee and U. V. Waghmare,
  *Wannier orbital overlap population (WOOP), Wannier orbital position population (WOPP) and the origin of anomalous dynamical charges*,[Phys. Chem. Chem. Phys., 2010, 12, 1564–1570](http://dx.doi.org/10.1039/b918890h)

  K. Terakura, S. Ishibashi,
  *Mechanism of covalency-induced electric polarization within the framework of maximally localized Wannier orbitals*,[Phys. Rev. B 91, 195120 (2015)](http://dx.doi.org/10.1103/PhysRevB.91.195120)


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/Chengcheng-Xiao/WOOPs/tags).

## Authors

* **Chengcheng XIAO** - *Initial work* - [E-mail](iconxicon@me.com)
## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/Chengcheng-Xiao/WOOPs/blob/master/LICENSE.md) file for details
