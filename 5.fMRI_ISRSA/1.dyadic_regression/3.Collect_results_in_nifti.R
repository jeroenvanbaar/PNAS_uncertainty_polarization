library(oro.nifti)
library(tidyverse)

# Arguments
args = commandArgs(trailingOnly=TRUE)
run = as.numeric(args[1])
filter_TR = as.numeric(args[2])
test_family = args[3]
if (filter_TR == 1) {
  TR_start = as.numeric(args[3])
  TR_end = as.numeric(args[4])
}
n_vox <- 1000

# Paths
if (test_family == 'summary'){
  data_dir = "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Results/voxelwise_ISRSA/raw"
  out_dir = "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Results/voxelwise_ISRSA/nifti"
} else if (test_family == 'anova'){
  data_dir = "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Results/voxelwise_ISRSA/raw_with_anova"
  out_dir = "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Results/voxelwise_ISRSA/nifti_with_anova"
} else {
  stop(paste("Invalid test family specified:",test_family))
}
brain_mask_img <- "/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/voxelwise_analyses/brain_map_consistency/80pct_brain_mask.nii.gz"
                
# Load dummy nifti to find brain map
print('Loading dummy nifti to get empty brain map...')
brain_mask <- readNIfTI(brain_mask_img)
imageDim <- dim(brain_mask)[1:3]
mask3D <- array(1, imageDim)
mask3D <- mask3D * (brain_mask[,,] != 0)
mask_array <- which(mask3D==1, arr.ind = TRUE)
mask <- which(mask3D == 1)
brain_mask[mask] <- 0
print('...loaded.')

# Get data
print('Extracting data from voxel group results...')
initialized <- FALSE
dats <- list()
vox_grs <- c(0:71)
for (voxgr in vox_grs){
  if (voxgr < 65){
    if ((voxgr %% 10) == 0){
      print(voxgr)
    }
  }
  else{
    print(voxgr)
  }
  # Read results data per voxel group
  if (filter_TR == 1){
    fname <- paste(data_dir,"/ISRSA_results_run-",run,"_TRs-", TR_start, "-", TR_end, "_voxgr-",voxgr,".RDS",sep='')
  } else {
    fname <- paste(data_dir,"/ISRSA_results_run-",run,"_voxgr-",voxgr,".RDS",sep='')
  }
  dats[[voxgr + 1]] <-
    readRDS(fname) %>%
    select(vox, x, y, z, model_name, term, test_family, estimate, statistic, p.value) %>% 
    filter((str_detect(term, "sd_", negate = TRUE)))
}
dat <- bind_rows(dats)
print('...done.')

# Correct p-values with FDR method
dat <-
  dat %>%
  group_by(model_name, term, test_family) %>%
  mutate(q = p.adjust(p.value, method = "fdr")
)

# Generate niftis
print('Writing niftis...')
model_name_list <- dat %>% pull('model_name') %>% unique()
for (model_name_i in model_name_list){
  print(sprintf('...for model %s',model_name_i))
  if (filter_TR == 1){
    fname <- paste(data_dir,"/ISRSA_results_run-",run,"_TRs-", TR_start, "-", TR_end, "_voxgr-",voxgr,".RDS",sep='')
    model_out_dir <- str_glue(out_dir, sprintf('/run-%i_TRs-%i-%i_model-%s', run, TR_start, TR_end, model_name_i))
  } else {
    model_out_dir <- str_glue(out_dir, sprintf('/run-%i_model-%s', run, model_name_i))
  }
  dir.create(model_out_dir)
  model_dat <- dat %>% filter(model_name == model_name_i)
  term_list <- model_dat %>% pull('term') %>% unique()
  for (term_i in term_list){
    # print(term_i)
    term_dat <- dat %>% filter(model_name == model_name_i & term == term_i)
    print(sprintf('......with term %s',term_i))
    out_prefix <- str_glue(model_out_dir,'/dyad_ISRSA_regressor-',term_i)
    out_prefix <- gsub(':','-X-',out_prefix) # Keep colons out of filenames
    # If available, write anova output
    test_families <- term_dat %>% pull(test_family) %>% unique()
    if ('anova' %in% test_families){
      print('Writing anova results')
      anova_dat <- term_dat %>% filter(test_family == 'anova')
      # F map
      dat_mask <- anova_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      fmap <- brain_mask
      fmap[dat_mask] <- anova_dat %>% pull(statistic)
      fname <- str_glue(out_prefix,'_F')
      if (!(file.exists(fname))) {
        print('............writing F map')
        writeNIfTI(fmap, fname, verbose = FALSE)
      } else {
        print('............skipping F map (already exists)')
      }
      # ANOVA P map
      dat_mask <- anova_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      anova_pmap <- brain_mask
      anova_pmap[dat_mask] <- anova_dat %>% pull(p.value)
      fname <- str_glue(out_prefix,'_anova-pval')
      if (!(file.exists(fname))) {
        print('............writing anova p map')
        writeNIfTI(anova_pmap, fname, verbose = FALSE)
      } else {
        print('............skipping anova p map (already exists)')
      }
      # Tresholded F map, p < 0.001 uncorrected
      unc_thresh_dat <- anova_dat %>% filter(p.value <= 0.001)
      dat_mask <- unc_thresh_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      unc_thresh_fmap <- brain_mask
      unc_thresh_fmap[dat_mask] <- unc_thresh_dat %>% pull(statistic)
      fname <- str_glue(out_prefix,'_F-thr-pval-unc-0.001')
      if (!(file.exists(fname))) {
        print('............writing thresholded F map p(unc) < 0.001')
        writeNIfTI(unc_thresh_fmap, fname, verbose = FALSE)
      } else {
        print('............skipping thresholded F map p(unc) < 0.001 (already exists)')
      }
      # FDR
      fdr_thresh_dat <- anova_dat %>% filter(q <= 0.05)
      dat_mask <- fdr_thresh_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      fdr_thresh_fmap <- brain_mask
      fdr_thresh_fmap[dat_mask] <- fdr_thresh_dat %>% pull(statistic)
      fname <- str_glue(out_prefix,'_F-thr-pval-fdr-0.05')
      if (!(file.exists(fname))) {
        print('............writing thresholded F map p(FDR) < 0.05')
        writeNIfTI(fdr_thresh_fmap, fname, verbose = FALSE)
      } else {
        print('............skipping thresholded F map p(FDR) < 0.05 (already exists)')
      }
    }
    ####
    # If available, write summary output
    if ('summary' %in% test_families){
      print('Writing summary results')
      summary_dat <- term_dat %>% filter(test_family == 'summary')
      # Stat map
      dat_mask <- summary_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      tmap <- brain_mask
      tmap[dat_mask] <- summary_dat %>% pull(estimate)
      fname <- str_glue(out_prefix,'_beta')
      if (!(file.exists(fname))) {
        print('............writing beta map')
        writeNIfTI(tmap, fname, verbose = FALSE)
      } else {
        print('............skipping beta map (already exists)')
      }
      # P map
      dat_mask <- summary_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      pmap <- brain_mask
      pmap[dat_mask] <- summary_dat %>% pull(p.value)
      fname <- str_glue(out_prefix,'_summary-pval')
      if (!(file.exists(fname))) {
        print('............writing summary p map')
        writeNIfTI(pmap, fname, verbose = FALSE)
      } else {
        print('............skipping summary p map (already exists)')
      }
      # Tresholded stat map, p < 0.001 uncorrected
      unc_thresh_dat <- summary_dat %>% filter(p.value <= 0.001)
      dat_mask <- unc_thresh_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      unc_thresh_map <- brain_mask
      unc_thresh_map[dat_mask] <- unc_thresh_dat %>% pull(estimate)
      fname <- str_glue(out_prefix,'_beta-thr-pval-unc-0.001')
      if (!(file.exists(fname))) {
        print('............writing thresholded beta map p(unc) < 0.001')
        writeNIfTI(unc_thresh_map, fname, verbose = FALSE)
      } else {
        print('............skipping thresholded beta map p(unc) < 0.001 (already exists)')
      }
      # FDR
      fdr_thresh_dat <- summary_dat %>% filter(q <= 0.05)
      dat_mask <- fdr_thresh_dat %>% ungroup %>% select(x, y, z) %>% data.matrix
      fdr_thresh_map <- brain_mask
      fdr_thresh_map[dat_mask] <- fdr_thresh_dat %>% pull(estimate)
      fname <- str_glue(out_prefix,'_beta-thr-pval-fdr-0.05')
      if (!(file.exists(fname))) {
        print('............writing thresholded beta map p(FDR) < 0.05')
        writeNIfTI(fdr_thresh_map, fname, verbose = FALSE)
      } else {
        print('............skipping thresholded beta map p(FDR) < 0.05 (already exists)')
      }
    }
  }
}
print('...done.')
