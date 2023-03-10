# CondensatesDynamicsML
VAMPNet application in deep learning analysis of binding peptides to condensates.

Protein condensates are actively studied in physics (active of active emulsions, polymer physics & soft matter), material science (controlling chemical reactions, desigining new materials from the bottom up) and biology (role in many diseases, possible roles in the origin of life as reaction crucibles). Therefore, analysing condensates dynamics and automating states transitions detection is very useful. I use [VAMPnet](https://github.com/markovmodel/deeptime/tree/master/vampnet) open source Python package to detect binding events in condensate-peptide problem. The dataset was simulated with HOOMD.


![Snapshot](https://github.com/egorys2000/CondensatesDynamicsML/blob/main/img/CondensatePeptideSnapshot.jpg)
![MarkowPlot](https://github.com/egorys2000/CondensatesDynamicsML/blob/main/img/MarkovPlot.jpg)

# Preparations
Create the conda environment from the tfenv.yml file with the following commands:

```
conda env create -f tfenv.yml
conda activate tfenv
```

# Usage
To test on test-data, launch Protein_condensate_model_training.ipynb and launch all cells above the **"Run several model iterations saving the best one, to help finding sparcely populated states"** section. To train current model, launch that section. To load weights, run **"Recover the saved model and its training history"** section. To get model output, run
```
temp = model.predict([traj_ord, traj_ord_lag], batch_size=np.shape(X1_vali)[0])
```

# Dataset
The data was simulated with HOOMD-blue package. The pipeline was the following:  
Simulate the formation of ctd condensate ->  
Add peptide chain (E10/P10 in test data) ->  
Run simulation to get binding events ->  
Process .dcd files to .npy and plug into neural network.  

The simulation script, topology file and other model parameters are in *hoomd* folder.

The whole dataset used for training can be found on Zenodo pages:  
- [Ctd condensate + E10 Peptide trajectories](https://zenodo.org/deposit/7474667)
- [Ctd condensate + P10 Peptide trajectories](https://zenodo.org/deposit/7474706)


# Further steps
As for now, the model input is just flattened vector of simulation snapshot, so the network loses geometric meaning of problem and can't get corresponding symmetries. This also forces us to fine-tune our model for new simulations with different peptides because of neurons domination. The next model may be based on PointNet++ architecture with CNN layers to impose rotational and translational equivariance.

# References:
Mardt, A., Pasquali, L., Wu, H., & No??, F. (2018). 
VAMPnets for deep learning of molecular kinetics. 
Nature communications, 9(1), 5.
