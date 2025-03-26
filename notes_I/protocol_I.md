# prepare (place in robot)

### Information
- This script facilitates the measurement of killcurves.
    - each assay plate represents a timepoint for two antibiotics. 
    - you can perfom at max 2 assays during one run (2x2 antibiotics)
- ADJUST Drugs.xlsx (in plate files folder) to change drugs, stock concentrations, working concentrations

### Prepare
- Configure your experiments using the configuration and setup files.
- Run the main.py file
- Grow O/N strains (3-rep with at least 1ml)
- PBS
- LB
- Agarplates (check setup.agarplates to see how many)
- White 384-plates (check setup.assayplates to see how many)


### Setup
- fill the antibiotic plate as defined in plates.xlsx
- fill the 96-deepwell plate as defined in plates.xlsx
- fill the PBS Trough with cool PBS.
- place everything into the robot as defined in the locations.csv (output folder)
- place agar and assayplates into the incubator as definde in storex_locations.csv (output folder)


### Starting cultures
- Robot will pick the ON cultures and dilute them 1:100
- Subsequently it will incubate for 1h
- Then the cultures are treated and the experiment starts.


### Agarplates (2l)
- Prepare 2l LB Agar
- Add 4ml Nystatin stock


### Nystatin Stock solution
  - Stock solution $5mg/ml$ (DMSO)
  - Working concentration: $10\mu g /ml$
