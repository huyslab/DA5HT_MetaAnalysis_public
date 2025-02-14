# This script contains high level functions to simplify reading analysis_script.Rmd

# this function preprocesses the references in the paper to allow us to match correct reference number with correct study in forest plots (JAMA psychiatry policy)
match_ref <- function(reference_file, datALL) {
  # Load the Excel file
  file_path <- reference_file #"author_ref_numbers.xlsx" 
  df_ref <- read_excel(file_path)
  
  # Function to extract the first author and year
  extract_author_and_year <- function(reference) {
    # Extract the first author (before the first space)
    author <- str_extract(reference, "^[^ ]+") 
    
    # Extract the year (four digits after a journal name)
    year <- str_extract(reference, "\\.\\s*([1-2][0-9]{3})[;\\(]") %>%
      str_extract("[1-2][0-9]{3}")  # Extract only the year
    
    # Handle missing values
    if (is.na(author)) author <- "Unknown"
    if (is.na(year)) year <- "Unknown"
    
    # Apply specific modifications - which it does not correctly amend
    if (author == "Erfanian") {
      author <- "Erfanian Abdoust"
    } else if (author == "de") {
      author <- "de Wit"
    } else if (author == "Smith" && year == "Unknown") {
      year <- "2024"
    } else if (author == "Westbrook" && year == "Unknown") {
      year <- "2024"
    }
    
    return(paste(author, year, sep = ", "))
  }
  
  # Apply the function to the second column
  df_ref <- df_ref %>%
    rename(ref_no = 1) %>% 
    mutate(`AuthorY` = sapply(df_ref[[2]], extract_author_and_year)) %>%
    relocate(`AuthorY`, .after = 1)  # Move new column to the correct position
  
  # Save the modified dataframe to a new Excel file
  # write.xlsx(df_ref, "updated_author_ref_numbers.xlsx")
  
  # Matching "AuthorY" in datALL with "Author, Year" in df_ref and creating 'ref_no' column
  # Extract reference number (assumed to be in the first column)
  df_ref_select <- df_ref %>%
    select(`AuthorY`, ref_no)
  
  # Perform the join, ensuring missing matches get NA
  datALL_new <- datAll %>%
    left_join(df_ref_select, by = "AuthorY")
  
  return(datALL_new)  # Return the modified datALL with the ref_no column
}



# This function computes d, vard, and r by design type and available stats,
# and coalesces into columns for d, vard, and r

compare_d_vard <- function(datAll) {
  datAll <- datAll %>%
    mutate(
      # Apply the compare_effectsize_calcs function to each row using pmap
      results = pmap(
        list(N_drug, N_placebo, drug_m, placebo_m, drug_sd, placebo_sd, r_correlation, treat_as, t),
        compare_effectsize_calcs
      ),
      # Extract each of the elements from the list into new columns
      d_rm_orig = map_dbl(results, 1),  
      vard_rm_orig = map_dbl(results, 2), 
      d_rm_stackx = map_dbl(results, 3), 
      vard_rm_stackx = map_dbl(results, 4),  
      d_rm_toster_eq_nct = map_dbl(results, 5),  
      
      vard_rm_toster_eq_nct = map_dbl(results, 6),  
      d_rm_toster_nct= map_dbl(results, 7),  
      vard_rm_toster_nct = map_dbl(results, 8),  
      # d_av_cheung = map_dbl(results, 9), 
      # vard_av_cheung = map_dbl(results, 10), 
      
      d_av_algina = map_dbl(results, 9),  
      vard_av_algina = map_dbl(results, 10), 
      d_rm_algina = map_dbl(results, 11),
      vard_rm_algina = map_dbl(results, 12),
      d_bw = map_dbl(results, 13),
      vard_bw = map_dbl(results, 14),
      d_z = map_dbl(results, 15),
      vard_z = map_dbl(results, 16),
    ) %>%
    select(-results)  # Optionally, remove the intermediate results column

  return(datAll)
}

compute_r <- function(datAll) {
  datAll <- datAll %>%
    mutate(
      r_correlation = if_else(
        treat_as == 'w' & nchar(author_provided_data) == 0 & !is.na(t),
        extract_r_from_t_means_sds(N_drug, t, drug_m, placebo_m, drug_sd, placebo_sd),
        if_else(
          treat_as == 'w' & nchar(author_provided_data) == 0 & is.na(t) & !is.na(SED),
          extract_r_from_sdd_sds(SED, drug_sd, placebo_sd),
          r_correlation  # Keep the existing value if neither condition is met
        )
      )
    )
  
  return(datAll)
}
  
compute_d_vard <- function(datAll) {
  # Compute d, vard, (r when applicable)
  datAll <- datAll %>%
    mutate(
      abs_d_from_means_sds = pmap(
        list(
          treat_as,
          N_drug,
          N_placebo,
          drug_m,
          placebo_m,
          drug_sd,
          placebo_sd
        ),
        abs_d_from_means_sds
      ),
      abs_d_from_t = pmap(
        list(
          treat_as,
          N_drug,
          N_placebo,
          drug_m,
          placebo_m,
          drug_sd,
          placebo_sd,
          t,
          w_method,
          r_correlation,
          t
        ),
        abs_d_from_t
      ),
      abs_d_from_F = pmap(
        list(
          treat_as,
          N_drug,
          N_placebo,
          drug_m,
          placebo_m,
          drug_sd,
          placebo_sd,
          `F`,
          w_method,
          r_correlation,
          t
        ),
        abs_d_from_F
      )
    )
  
  # Split d, vard
  datAll <- datAll %>%
    mutate(
      vard_from_means_sds = map_dbl(abs_d_from_means_sds, "vard"),
      vard_from_t = map_dbl(abs_d_from_t, "vard"),
      vard_from_F = map_dbl(abs_d_from_F, "vard"),
      # r_from_t = map_dbl(abs_d_from_t, "r"),
      # r_from_F = map_dbl(abs_d_from_F, "r"),
      abs_d_from_means_sds = map_dbl(abs_d_from_means_sds, "d"),
      abs_d_from_t = map_dbl(abs_d_from_t, "d"),
      abs_d_from_F = map_dbl(abs_d_from_F, "d")
    )
  
  # Compute for studies reporting d
  datAll <- datAll %>%
    mutate(
      abs_d_from_d_raw = abs(d_raw),
      vard_from_d_raw = ifelse(design_true == "bw", 
                               bw_var_d(N_drug, N_placebo, d_raw),
                               NA)
    )
  
  # Perform checks
  assert("There are between participant papers with insufficient stats",
         with(datAll, sum(design_true == "bw" &
                            is.na(abs_d_from_means_sds) &
                            is.na(abs_d_from_t) &
                            is.na(abs_d_from_F) &
                            is.na(abs_d_from_d_raw)) == 0))
  
  
  # Merge and multiply by direction
  datAll <- datAll %>%
    mutate(
      abs_d = coalesce(
        abs_d_from_means_sds,
        abs_d_from_t,
        abs_d_from_F,
        abs_d_from_d_raw
      ),
      d = abs_d * d.direction.new,
      vard = coalesce(
        vard_from_means_sds,
        vard_from_t,
        vard_from_F,
        vard_from_d_raw
      # ),
      # r = coalesce(
      #   r_from_t,
      #   r_from_F
      )
    )
  
  return(datAll)
}
