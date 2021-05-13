library(oro.nifti)
library(tidyverse)
library(widyr)
library(broom)
library(broom.mixed)
library(lme4)
library(lmerTest)
library(modelr)

subjects_dir <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Data/Subjects_and_exclusions"
data_dir <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/derivatives/cleaning"
predictor_RDM_dir <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Data/Cleaned/Surveys"
predictor_RDMs <- read_csv(file.path(predictor_RDM_dir, "predictor_RDMs_4.csv"))
output_dir <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Results/voxelwise_ISRSA/raw_with_anova"
brain_mask_img <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/voxelwise_analyses/brain_map_consistency/80pct_brain_mask.nii.gz"
args = commandArgs(trailingOnly=TRUE)

# run_lmer_dyad function --------------------------------------------------
run_lmer_dyad <- function(data, formula, dyad_1, dyad_2, drop_Inf = FALSE, lmerTest = FALSE) {
  hier_formula <- update(formula, ~ . + (1 | dyad))
  
  mod_terms <- all.vars(subbars(formula))
  
  #model only uses rows with no NAs so need to select only members of dyads in those rows
  model_df_na_omit <- data %>% dplyr::select(tidyselect::all_of(mod_terms), {{dyad_1}}, {{dyad_2}}) %>% 
    drop_na()
  if(drop_Inf) {
    model_df_na_omit <- model_df_na_omit %>% filter_all(all_vars(. != Inf))
  }
  
  #only use factor levels that actually appear in dataset
  dyad_levels_sub <- fct_c(model_df_na_omit %>% pull({{dyad_1}}) %>% factor, 
                           model_df_na_omit %>% pull({{dyad_2}}) %>% factor) %>% levels()
  #Create new column that includes all the levels (this is used later in creating Zt matrix)
  model_df_na_omit <- model_df_na_omit %>% mutate(
    dyad = rep(dyad_levels_sub, length.out=dim(model_df_na_omit)[1]),
    dyad_1_fct = factor({{dyad_1}}, levels = dyad_levels_sub),
    dyad_2_fct = factor({{dyad_2}}, levels = dyad_levels_sub)
  )
  
  #Fill matrix with dyad 1 and then add together
  dyad_1_mat <- model_df_na_omit %>% #model_df %>%
    modelr::model_matrix(~ 0 + dyad_1_fct) %>%
    dplyr::select(starts_with("dyad_1_fct")) %>%
    rename_all(list(~(str_replace(., fixed("dyad_1_fct"), ""))))
  dyad_2_mat  <- model_df_na_omit %>% #model_df %>%
    modelr::model_matrix(~ 0 + dyad_2_fct) %>%
    dplyr::select(starts_with("dyad_2_fct")) %>%
    rename_all(list(~(str_replace(., fixed("dyad_2_fct"), ""))))
  dyad_mat <- bind_rows(dyad_1_mat %>% rownames_to_column(),
                        dyad_2_mat %>% rownames_to_column()) %>%
    group_by(rowname) %>%
    summarise_all(sum) %>% 
    mutate(rowname_num = as.numeric(rowname)) %>%
    arrange(rowname_num) %>%
    dplyr::select(-rowname, -rowname_num)
  
  lmod <- lFormula(hier_formula, data = model_df_na_omit)
  lmod$reTrms$Ztlist[['1 | dyad']] <- Matrix(t(dyad_mat))
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun)
  res <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  if(lmerTest) {
    res <- lmerTest:::as_lmerModLT(res, devfun)
  }
  return(res)
}


# load data ---------------------------------------------------------------

array_ind <- as.numeric(args[1])
ses <- as.numeric(args[2])
run <- as.numeric(args[3])
n_vox <- as.numeric(args[4])
filter_TR <- as.numeric(args[5])
if (filter_TR == 1){
  TR_start <- as.numeric(args[6])
  TR_end <- as.numeric(args[7])
  print(sprintf('Session %i, run %i, TRs %i through %i',ses,run, TR_start, TR_end))
} else{
  TR_start <- 1
  print(sprintf('Session %i, run %i',ses,run))
}

vox_sel <- (1 + (n_vox * array_ind)):(n_vox * (array_ind + 1))
ses_string <- str_glue("ses", ses, .sep = "-")
run_string <- str_glue("run", run, .sep = "-")

data_tbls <- c() 
ind <- 1

brain_mask <- readNIfTI(brain_mask_img)
imageDim <- dim(brain_mask)[1:3]
mask3D <- array(1, imageDim)
mask3D <- mask3D * (brain_mask[,,] != 0)
mask_arr.ind <- which(mask3D == 1, arr.ind = TRUE)
if (vox_sel[length(vox_sel)] > dim(mask_arr.ind)[1]){
  vox_sel <- vox_sel[vox_sel <= dim(mask_arr.ind)[1]]
}
if (length(vox_sel) == 0){
  stop('No voxels within brain mask for this vox group')
}

# Determine subjects to include
all_subs <- read_csv(str_glue(subjects_dir, "/all_subjects.csv")) %>% pull(sub)
exclusions <- read_csv(str_glue(subjects_dir, sprintf("/exclude_video-watching_aggregate_run-%i.csv",run))) %>% pull(sub)
print('Exclusions: ')
print(exclusions)
subs <- all_subs [! all_subs %in% exclusions]

for (sub in subs) {
  # print(sub)
  sub_string <- str_glue("sub", str_pad(sub, 3, pad = 0), .sep = "-")
  
  fname <- file.path(data_dir, 
                     sub_string,
                     ses_string,
                     "func", 
                     str_glue(sub_string,
                              ses_string,
                              "task-videoWatching",
                              run_string,
                              "space-MNI152NLin2009cAsym_desc-cleaned_bold.nii.gz", .sep = "_"))
  # cat("Loading ", fname)
  print(sprintf("Loading data subject %i: %s", sub, fname))
  img <- readNIfTI(fname)
  noScans <- dim(img)[4]
  if (filter_TR == 0){
    TR_end <- noScans
    print('Including all TRs')
  } else{
    if (TR_end > noScans){
      TR_end <- noScans
    }
    print(sprintf('Including TRs %i through %i', TR_start, TR_end))
  }
  
  for(t in TR_start:TR_end) {
    scan <- img[,,,t]
    
    data_tbls[[ind]] <- tibble(
      sub = sub,
      t = t,
      vox = vox_sel,
      x = mask_arr.ind[vox_sel, 1],
      y = mask_arr.ind[vox_sel, 2],
      z = mask_arr.ind[vox_sel, 3],
      BOLD = scan[mask_arr.ind[vox_sel, ]]
    )
    ind = ind + 1
  }
}

vox_df <- bind_rows(data_tbls)
vox_df <- vox_df %>% 
  filter(t > 4) # Remove the first 4 volumes (countdown timer before video)

ISC_RDM_df <- vox_df %>% 
  group_by(vox, x, y, z) %>%
  nest() %>%
  mutate(
    pairwise_cor = map(data, pairwise_cor, item = sub, feature = t, value = BOLD, upper = FALSE)
  ) %>% unnest(pairwise_cor) %>%
  select(-data) %>%
  rename(SubID1 = item1, SubID2 = item2, ISC = correlation)          

ISRSA_df <- ISC_RDM_df %>% ungroup() %>%
  left_join(predictor_RDMs, by = c("SubID1", "SubID2"))

# Define models to run
control_f <- ~scale(age_distance) + scale(scan_day_distance) + same_gender + same_undergrad + same_community
formulas_forward <- formulas(~scale(ISC),
                ideology = ~scale(ideology_similarity),
                ideology_IUS = ~scale(ideology_similarity) * joint_IUS)
ffs <- formulas_forward

# Check if data file exists
if (filter_TR == 1) {
  out_fname <- str_glue(output_dir, "/ISRSA_results_run-", run, "_TRs-", TR_start, "-", TR_end,"_voxgr-", array_ind, ".RDS")
} else {
  out_fname <- str_glue(output_dir, "/ISRSA_results_run-", run, "_voxgr-", array_ind, ".RDS")
}
print(sprintf('Saving to %s',out_fname))
if (file.exists(out_fname)){
  existing_RDS <- readRDS(out_fname)
  existing_models <- unique(existing_RDS$model_name)
  new_models <- setdiff(names(ffs),existing_models)
  excluded_models <- setdiff(names(ffs), new_models)
  if (length(excluded_models) > 0){
    print(paste("SKIPPING models", paste(excluded_models, collapse = ','), "as they exist in the data on disk.\n"))
  }
  ffs <- ffs[new_models]
}

# Run models here
print("RUNNING the following models:")
print(ffs)

# Test
model_df <- ISRSA_df %>% 
  group_by(vox, x, y, z) %>%
  nest() %>%
  mutate(
    mod_df = map(data, fit_with, run_lmer_dyad, ffs,
                 dyad_1 = quo(SubID1), dyad_2 = quo(SubID2),
                 lmerTest = TRUE),
    enframed_mod_df = map(mod_df, enframe, name = "model_name", value = "mod")
  ) %>% unnest(enframed_mod_df) %>% select(-mod_df) %>% mutate(
    singular = Vectorize(isSingular)(mod),
    rePCA = Vectorize(rePCA)(mod),
    mod_anova = map(mod, anova),
    tidied = map(mod, tidy),
    tidied_anova = map(mod_anova, tidy),
    glanced = map(mod, glance)
  )
model_df_summary <-
  model_df %>%
  unnest(tidied) %>%
  select(-data, -mod, -mod_anova, -tidied_anova, -glanced, -rePCA) %>%
  mutate(test_family = 'summary')
model_df_anova <-
  model_df %>%
  unnest(tidied_anova) %>%
  select(-data, -mod, -mod_anova, -tidied, -glanced, -rePCA) %>%
  mutate(test_family = 'anova')
save_df <-
  bind_rows(model_df_summary, model_df_anova)

# Append to existing data, if there are
if (file.exists(out_fname)){
  save_df <- rbind(existing_RDS, save_df)
}

print('saving...')
saveRDS(save_df, file = out_fname)
print('Done.')