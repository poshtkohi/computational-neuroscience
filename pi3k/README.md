# Mathematical Modelling of PI3K/Akt Pathway in Microglia

This is a home page regarding the MATLAB source codes for a novel biophysical model for cytosolic Ca<sup>2+</sup> activation of the PI3K/Akt pathway in microglia where Ca<sup>2+</sup> influx is mediated by both P2Y purinergic receptors (P2YR) and P2X purinergic receptors (P2XR).

This respiratory has two models as follows.

## 1. P2Y<sub>12</sub> Model
It develops a new IP<sub>3</sub>R (IP<sub>3</sub> receptor) model along with a complete biophysical model to capture human microglial Ca<sup>2+</sup> data from activation of P2Y<sub>12</sub> receptor to Ca<sup>2+</sup> release from intracellular stores consisting of SERCA and Ca<sup>2+</sup> leak channel.

This model has several files, among the most important of which are explained below:

1. The experimental data used can be found under the directory _data_.
2. The reaction network for P2Y and IP3R models resides in the file _reaction_network_hp2y12.m_.
3. SECRA, IP3R and leak currents are implemented as three separate functions in the files _secra_current.m_, _ip3r_current.m_ and _leak_current.m_.
4. The fitting process which includes a mathematical optimisation function using an evolutionary strategy for numerical simulation of the model differential equations can be found in several files including _cmaes.m_, _fhngen_hp2y12_fit.m_, _loss_function.m_ and _main_hp2y12_fit.m_.
5. The model predictions are implemented in the file _hp2y12_cai_prediction.m_.

## 2. PI<sub>3</sub>K/Akt Model
It builds a new PI<sub>3</sub>K/AKT model that uses a time dependant [Ca<sup>2+</sup>]<sub>i</sub> profile to activate the PI<sub>3</sub>K/pAkt via CaMKII. This model supports that calcium influx via P2YR and P2XR can justify the experimentally observed twin peaks in the pAkt data.

