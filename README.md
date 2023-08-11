## ThreeFieldMagnetoElasticity
The three-field finite element procedure for hard-magnetic soft materials

### Purpose 
This repository provides the source code and the input files for the numerical examples used in the paper draft titled “A three-field based finite element analysis for a class of magnetoelastic materials”. The repository contains the following content:
1. the source code of the three-field finite element procedure for magneticelastic material modeling.
2. the input files for the three numerical examples used in the aforementioned manuscript.

### How to compile
The three-field finite element procedure is implemented in [deal.II](https://www.dealii.org/) (version 9.4.0), which is an open source finite element library. In order to use the code (**main.cc**) provided here, deal.II should be configured with MPI and at least with the interfaces to BLAS, LAPACK, Threading Building Blocks (TBB), and UMFPACK. For optional interfaces to other software packages, see https://www.dealii.org/developer/readme.html.

Once the deal.II library is compiled, for instance, to "~/dealii-9.4.0/bin/", follow the steps listed below:
1. cd SourceCode
2. cmake -DDEAL_II_DIR=~/dealii-9.4.0/bin/  .
3. make debug or make release
4. make

### How to run
1. copy the input files contained in one of the example folders into the folder SourceCode
2. make run

### Reference
This work is based on the constitutive formulation proposed in the paper:
Zhao, R., Kim, Y., Chester, S.A., Sharma, P., Zhao, X., 2019. Mechanics of hard-magnetic soft materials. Journal of the Mechanics and Physics of Solids, 124, 244–263. doi:https://doi.org/10.1016/j.jmps.2018.10.008.

The three-field finite element implementation is developed by using Step-44 of the deal.II library as the basis:
Pelteret, Jean-Paul, & McBride, Andrew. (2012). The deal.II tutorial step-44: Three-field formulation for non-linear solid mechanics. Zenodo. https://doi.org/10.5281/zenodo.439772.
