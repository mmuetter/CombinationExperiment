# III pd curve

- prep 12 col columns [1, 4, 7, 10] with 10 ml LB and infect from colonies.
- place all plates in the robot according to
  - locations.csv
  - storex_locations.csv
- This script uses the pintool. Prep the pintool station...
- Theoretically that script works for an unlimited number of drug combinations (as defined in combination_idx = [0], which refers to the keys of the config.combinations dict), however:
  + for each drug combintation (incl. plate I and II) we need one tip block, which needs to be manually replaced during the script (a prompt will show up)
  + measuring each assayplate likely takes 5 minutes, making the cycle for each drug combination 10 minutes. Therefore it is likely not advisable to test more than one drug combination at once.