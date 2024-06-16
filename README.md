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
- We construct the network data into `.xlsx`, with the mean and variance of travel cost responding to each edge.
- Please unzip the datasets in the folder **experiment** before you attempt to run the code

## Dependency
- **Graph**:
1. Open the MATLAB 'Add-On Explorer'. You can do this by:  
   Selecting "Home" -> "Get Add-Ons" from the main menu bar **or** typing `addons.open()` in the command window
2. In the 'Add-On Explorer', search for "graph" and find the "MATLAB Graph Package".
3. Click the "Install" button to start downloading and installing the graph package.

## Notes
- The default representation of networks is the adjacency matrix.
- The default sampling model is "Gaussian".
- Support custom policy function.
- Contact details: shiruiaa@foxmail.com
