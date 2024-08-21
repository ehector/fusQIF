The files in this folder can be used to reproduce the simulation results in the paper. It is highly recommended that these be run on a high-performance computing cluster. 

All simulations and analyses were run in R version 4.0.2. Please ensure this version of R is being used when attempting to reproduce results; notably, RcppArmadillo's default version of pinv was changed in version 10.3 and results in different behaviour. Before running any of the files in this folder, please install the fusQIF R package included in the supplement. 

The files are:

- normal_SettingI.R : this file generates one dataset from Setting I of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- normal_SettingII.R : this file generates one dataset from Setting II of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- normal_SettingIII.R : this file generates one dataset from Setting III of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- poisson_SettingI.R : this file generates one dataset from Setting I of the Poisson regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- poisson_SettingII.R : this file generates one dataset from Setting II of the Poisson regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- poisson_SettingIII.R : this file generates one dataset from Setting III of the Poisson regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- normal_SettingI.R : this file generates one dataset from Setting I of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- normal_SettingII.R : this file generates one dataset from Setting II of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.

- normal_SettingIII.R : this file generates one dataset from Setting III of the linear regression simulations and fits the proposed estimator to this dataset. To generate full simulation results, this script should be run 500 times with 500 different seeds.
