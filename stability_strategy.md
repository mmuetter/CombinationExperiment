# Antibiotic Combination Experiment Protocol
Setting up the experiment and the procedure of handling and making the drug combinations can have a critical impact on the variance and precision of the conducted experiments as some drugs, most notebly amoxicillin, are very sensitive and tend to degrade very fast.
There are three main causes of variance and inprecision we have to consider: 
- degradation in liquid at four degrees
- degradation due to freezing
- increased variance if the antibiotic is sourced from different original stocks and mixed and diluted at differnent days.
Degradation in liquid may cause a problem as we have to conduct 15 seperate pd-curve experiments and more than one per day was not very feasable. 
Considering that we may not run one experiment every day and sometimes have to repeat an experiment completing all experiments takes approx one month.
This leads to three possible approaches of handling the drugs.
1) Mixing only reservoirs and antibiotic plates when they are needed, and consume sensitve drugs first. 
  - Likely still a notable decrease of drug concentration for some drugs.
  - medium hard to setup (not batch work but we only have to do each stock and dilution row once)
2) Mixing every stock, reservoir and combination plate every day newly.
   - On average quite precies
   - increased variablity between experiments due to inprecisions in handling and weighing
   - Very much effort 
3) premake making all stocks and reservoirs and aliquot subreservoirs. Then freeze all of them at -80 degrees.
   - Likely slight underestimation of the real potency for most drugs due to degradation due to freezing
   - Very low interexperiment variation because  each drug is soureced from the same stocks and reservoirs for all combinations
   - Very easy and flexible to setup all experiments
  

As we are primarily interested in the drug interactions between pairwise combinations we prioritized internal consistency accross experiments over absolute potency.

This experiment investigates all **pairwise combinations** of 6 antibiotics:

- Amoxicillin (AMX)
- Chloramphenicol (CHL)
- Colistin sulfate (COL)
- Fosfomycin (FOS)
- Polymyxin B (POL)
- Tetracycline (TET)

The total number of combinations is:
\[
\binom{6}{2} = 15
\]

The primary goal is **internal consistency across experiments**, not absolute potency.  
Therefore, we prioritize:

- Freezing all drug stock plates **at once** at **−80 °C**
- Using the **same solvent** (20% glycerol in sterile ddH₂O) for all drugs
- Avoiding **day-to-day variation** and **freeze–thaw cycles**

---

## 🧬 Combination Strategy

For each combination, the drugs are assigned a role as **drug A** and **drug B**.

To modularize the experiment, we use the following pre-filled plates:

### Plate Types

| Plate      | Volume per well | Purpose                               |
| ---------- | --------------- | ------------------------------------- |
| A-sub      | 125 µL          | Gradient of drug A                    |
| B-I / B-II | 50 µL           | Gradient of drug B (duplicate plates) |

Each experiment involves:
- One A-sub plate for drug A
- Two B plates for drug B (B-I and B-II)
- Mixing: 50 µL is transferred from A-sub to each B plate before use
- Final volume per well (after mixing): 100 µL

All plates are prepared using **20% glycerol in ddH₂O**, aliquoted, sealed, and **frozen once at −80 °C**.  
Each plate is **used once** and thawed only on the day of the experiment.

---

## 🧊 Antibiotic Stability Considerations

| Antibiotic      | Stability at −80 °C in 20% Glycerol | Notes                                           |
| --------------- | ----------------------------------- | ----------------------------------------------- |
| Chloramphenicol | 🟢 Excellent                         | Chemically stable in aqueous solvents           |
| Polymyxin B     | 🟢 Very good                         | Stable at pH < 6; no known freeze–thaw issues   |
| Colistin        | 🟡 Moderate                          | Slightly less stable than PMB, still acceptable |
| Fosfomycin      | 🟡 Moderate                          | Susceptible to hydrolysis in solution           |
| Tetracycline    | 🔴 Limited                           | Sensitive to pH/light, but stable when frozen   |
| Amoxicillin     | 🔴 Unstable in solution              | Freeze immediately, avoid thawing delays        |

---

## 📦 Plate Count per Antibiotic
## 📦 Plate Count per Antibiotic

Below is the number of plates needed per drug, broken down into:  
- <span style="color:black"><strong>Black</strong></span> = minimum required  
- <span style="color:green"><strong>Green</strong></span> = extra plates for repetitions or errors

| Antibiotic      | A-sub (125 µL)                                                          | B-I / B-II (50 µL each)                                                 |
| --------------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Amoxicillin     | <span style="color:black">5</span> + <span style="color:green">3</span> | <span style="color:black">0</span>                                      |
| Chloramphenicol | <span style="color:black">4</span> + <span style="color:green">2</span> | <span style="color:black">1</span> + <span style="color:green">1</span> |
| Colistin        | <span style="color:black">3</span> + <span style="color:green">2</span> | <span style="color:black">2</span> + <span style="color:green">2</span> |
| Fosfomycin      | <span style="color:black">2</span> + <span style="color:green">2</span> | <span style="color:black">3</span> + <span style="color:green">2</span> |
| Polymyxin B     | <span style="color:black">1</span> + <span style="color:green">1</span> | <span style="color:black">4</span> + <span style="color:green">2</span> |
| Tetracycline    | <span style="color:black">0</span> + <span style="color:green">0</span> | <span style="color:black">5</span> + <span style="color:green">3</span> |

---

## ✅ Summary Guidelines

- Use **20% glycerol in sterile ddH₂O** as the universal solvent
- Aliquot A-sub and B plates in advance and freeze **once** at −80 °C
- Thaw each plate **only once**, on the day of use
- Add antibiotic to LB **on the day of the experiment**
- Final glycerol concentration in assay: ~1.5–2% (non-toxic)

This strategy ensures **stable, reproducible drug exposure** across all 15 pairwise combinations, allowing precise comparative analysis despite limited long-term stability of certain drugs.

