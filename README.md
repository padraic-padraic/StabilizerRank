# Stabilizer Rank Project

This is the main git repository for my MRes project on the Bravyi-Gosset Algrotihm for simulating quantum circuits, undertaken to fulfill the requirements of my MRes in Delivering Qunantum Technologies at UCL. 

The project is split into several chunks: 

1) **iPython Notebooks**
    These contain prototyping code for generating n-Qubit stabilizer groups, and a small test for bounding the size of approximate decompositions of n-fold magic states in the comp
utational basis. 

2) **'StabiliserRank'**
    This folder contains C-code to implement the pseudo-code descriped by Bravyi & Gosset in the appendices of their [paper](http://arxiv.org/abs/1601.07601). When extended with code to generate the correct sequence of Pauli measurements to implement a given circuit, this will implement the full StabiliserRank method
    
3) **Robustness**
    This folder contains python code implementing the `Robustness' measure of Mark Howard and Earl Campbell (Work not yet published), to bound the decompositions of n-qubit states in to sums of Stabiliser states

4) **Report**
    This folder contains Tex for generating the MRes report, and can be built using `latexmk --shell-escape -pdf report.tex`. 

