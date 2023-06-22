# eco-analytic

This code is used to generate results from the paper "Evaluating The Impact Of Species Specialisation On Ecological Network Robustness Using Analytic Methods" by C Jones, D Zurell and K Wiesner.

The module "eco_module" includes all functions necessary for replicating our results. The file "ecological_robustness_example" implements various functions in order to predict species survival on ecological networks for various extinction scenarios, and plots results similar to figures in our paper. The folder "large_poll_nets" contains several data files for real world plant pollinator interaction networks.

numpy, scipy, networkx and matplotlib packages are required. This code was originally written and run in Python 3.9.
