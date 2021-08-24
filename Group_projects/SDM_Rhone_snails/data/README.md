# Data Description

Data source: This data comes from the monitoring of the Rhône restauration project. 
Responsible: emmanuel.castella@unige.ch

Update: There is now also an .RData file that contains all the modeling data in one data frame.

There are 5 .csv tables, three main tables and two explanatory tables:

- gastero_fauna.csv
  - Contains count data
  - Rows -> Species
  - Columns -> Samples

- gastero_environment.csv
  - Contains environmental data
  - Rows -> environmental variables
  - Columns -> Samples
  
- gastero_samples.csv

  - Contains info about the samples
  - Rows -> samples
  - Columns -> Info 
    - Channel: channel
    - Site: AM: Upstream or AV: Downstream
    - Q: Sample square number within the sit
    - Year: Year
    - Season: E: Summer or P: Spring
 
- gastero_fauna_key.csv
  - Info about Species incl. Family
  
- gastero_environment_key.csv
  - Info about environmental variables
