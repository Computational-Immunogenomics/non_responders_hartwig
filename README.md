# Systematic identification of genomic non-response biomarkers to cancer therapies 


We share code to run analyses described here: https://www.biorxiv.org/content/10.1101/2025.10.06.680777v1)

---

# Overview

The underlying study data comes from the Hartwig Medical Database, and data access was obtained through an approved public data request (DR347). To request access to the data see instructions described on the Hartwig Medical Foundation website (https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/).   

---

# Folders and files

- [0_clinical](0_clinical/)
    - Process clinical data: define tumor type cohorts, organize treatments, derive outcomes 
- [1_database](1_database/)
    - Process data to make easy-to-use tables of various tools output across patients
- [2_biomarkers](2_biomarkers/)
    - Compute genomic biomarkers from processed data (1_database)
- [3_signals](3_signals/)
    - 0_analysis:  Prepare features, define treatment/tumor type cohorts, run systematic testing, clean output
    - 1_figures: Reproduce the figures shown in the study. 

- [helpers](helpers/)
    - Store helper functions used throughout the analyses
- map.r,map.py,run.py
    - Maps of folder locations and run file

---
# Contact
For questions: Joseph Usset - josephusset@vhio.net
