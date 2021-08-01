# Source code for *Dual Dynamic Programming for the Mean Standard Deviation Canadian Traveller Problem*
## Description
The source code includes the following files:
* **main.m**, the entry to the program.
* **simulator_for_realistic_networks.m**, a simulator for simulating observation and road information based on real datas.
* **DDL_DL.m**, function of the method we proposed.
* **Pi_tau.m**, function of the benchmark algorithm $\pi(\tau)$. 
* **AO_star.m**, function of the benchmark algorithm $AO*$. 
* **RAO_star.m**, function of the benchmark algorithm $RAO*$. 
* **Kahni_RSP.m**, function of the benchmark algorithm $DP$. 
* **Zyf_RSP.m**, function of the benchmark algorithm $SG$.
  
## Dataset
We construct the network data into `.xlsx`, with the mean and variance of travel cost responding to each edge.

## Notes
- The default representation of networks is the adjacency matrix.
- The default sampling model is "Gaussian".
- Support custom policy function.
- Contact details: shiruiaa@foxmail.com
