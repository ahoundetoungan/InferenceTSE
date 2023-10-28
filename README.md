## InferenceTSE
**Replication codes for *Inference for Two-Stage Extremum Estimators* By Aristide Houndetoungan and Abdoul Maoude**

### Simulations
- `tsesimu.cppfunctions.cpp` contains `C++` functions that are used in the R codes of the simulations.
- `tsesimu.rfunctions.R` contains additional `R` functions that are used in the simulations.
- `simu.dgpA.R` replicates the simulation results for DGP A.
- `simu.dgpB.R` replicates the simulation results for DGP B.
- `simu.dgpC.R` replicates the simulation results for DGP C.

### Application: peer effects on adolescent smoking habits
- `smoking-functions.cpp` contains `C++` functions that are used in the R codes of the application.
- `0_Inschool.do` extracts the part of the data set to be used from the Add Health data set.
- `1_smoking.data.R` prepares the data set to be used.
- `2_smoking.model.R` replicates the peer effect model estimation.