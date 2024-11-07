# Molecular Properties Calculation for Electrolyte Research

This project focuses on calculating key properties of molecular components in electrolytes to aid research on electrolyte composition. The following table lists the calculated properties and the tools/methods used for each:

| Property                             | Calculation Method  |
|--------------------------------------|---------------------|
| Maximum Electrostatic Potential      | Psi4               |
| Minimum Electrostatic Potential      | Psi4               |
| HOMO                                 | Psi4               |
| LUMO                                 | Psi4               |
| Polar Surface Area                   | RDKit              |
| Solvent Accessible Surface Area      | RDKit              |
| Gutmann's Donor Number               | ML Prediction      |
| V<sub>red</sub>                      | ML Prediction      |
| V<sub>ox</sub>                       | ML Prediction      |

### Requirements

1. **DFT Calculations**: The DFT calculations rely on the Psi4 quantum chemistry software. You can refer to [espsim repository](https://github.com/hesther/espsim) for additional setup instructions.
  
2. **Machine Learning Model**: The machine learning (ML) portion uses a model trained with MolCLR. Follow the instructions in the [MolCLR repository](https://github.com/yuyangw/MolCLR/tree/master) to set up the environment and train an MLP model for specific property predictions.

### Data Format

The `data.xlsx` file serves as an example, showing how to calculate properties for a large batch of candidate molecules.

### Usage

- **Electrostatic Potential Calculation**: Run `python esp.py` to calculate the electrostatic potential for each molecule.
  
- **HOMO-LUMO Calculation**: Run `python homo_lumo.py` to compute the HOMO and LUMO values.

- **Comprehensive Property Calculation**: Run `python infer.py` to calculate all properties listed in the table and organize the results into a comprehensive table.
