# DA & 5HT Meta Analysis

## Computational environment
The scripts in this folder run off a Docker container. One can run the analysis, 
or edit the scripts by opening an RStudio instance and accessing it through the 
browser. After installing Docker Desktop, navigate into the repository folder in the
Terminal, and then run the following command:

```
docker run -it --rm --name analysis_container -p 8787:8787 -v $(pwd):/home/rstudio/analysis -e PASSWORD=meta anamkr/da5ht_meta_analysis:v1.1

```

Then open your web browser (preferably Chrome) and navigate to 
`http://0.0.0.0:8787/`. The username is `rstudio`, and the password is `meta`.

## Files in this repository

```
.
└── analysis/
    ├── 2024_08_27_final_metaanalysis.csv // Data file listing all assessed studies
    ├── 2024_08_27_aversive_pavlovian_w.csv // Data file listing all 5HT within-subject aversive Pavlovian studies
    ├── 2024_08_27_mb_w.csv // Data file listing all DA within-subject model-based/flexibility studies
    ├── 2024_08_27_risk_w.csv // Data file listing all DA within-subject risk attitude studies
    ├── author_ref_numbers.xlsx // spreadsheet listing all reference numbers from manuscript (make sure Guitart-Masip 2024 review is removed). Used for adding reference numbers to manuscript forest plots
    ├── analysis_script_public.Rmd // Main analysis file
    ├── compute_effect_size.R // Functions to compute effect sizes from reported statistics
    ├── utils.R // Higher level functions used to simplify analysis_script_public.Rmd
    └── Results // Folder where all the figures are saved from running the analysis_script_public.Rmd
```  
