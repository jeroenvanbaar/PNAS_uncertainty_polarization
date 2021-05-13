library(tidyverse)
library(mediation)
library(lme4)
library(lmerTest)
library(vroom)
library(plotly)

# GENERAL SETTINGS
base_dir <- 'BASE_DIR_HERE'
library(MuMIn)
source(file.path('lmer_multimemb_reg.R'))

###################################################################################################################
# STEP 1: Load data
###################################################################################################################

## DEFINE MODELS/TERMS
run <- 3
model <- 'ideology_IUS'
term <- 'scale(ideology_similarity)-X-joint_IUS'

# Load data and scale where appropriate,compute interaction term
# This loads ROI data too which is not necessary but I copy-pasted from elsewhere.
var_to_plot <- 'mean_BOLD'
reg_data <- tibble(vroom(sprintf('%s/Results/voxelwise_ISC/nifti/run-%i_model-%s/ROI_regression_data_term-%s_%s.csv',
                                 base_dir, run, model, term, var_to_plot)))
relabeled <- transform(reg_data,
                       Pair_1ID = pmin(SubID1, SubID2),
                       Pair_2ID = pmax(SubID1, SubID2))
reg_data <- subset(relabeled, !duplicated(cbind(Pair_1ID, Pair_2ID)) & (SubID1 != SubID2)) %>% dplyr::select(-Pair_1ID, -Pair_2ID)
roi.cols <- names(reg_data)[10:(length(names(reg_data))-1)]
print(roi.cols)
judg.data <- vroom(sprintf('%s/Results/video_judgment/all_judgment_similarity.csv',base_dir))
DV.cols <- names(judg.data)[3:length(names(judg.data))]
path_data <- tibble(reg_data %>% left_join(judg.data)) # Merge
path_data_scale <- path_data %>% 
  mutate_at(c(3,5,6,10:19,21:27),scale) %>% 
  mutate_at(c(3,5,6,10:19,21:27),funs(as.numeric(.))) %>% 
  mutate(ideo_sim_x_joint_IUS = ideology_similarity*joint_IUS)
print(dim(path_data_scale))
path_data_scale$hi_IUS <- path_data_scale$joint_IUS > median(path_data_scale$joint_IUS)

### Select and reshape data
agreement_sim_df <- path_data_scale %>% dplyr::select(SubID1, SubID2, joint_IUS, ideology_similarity,
                                                      agreement_v2_sim, agreement_v3_sim) %>% 
  pivot_longer(cols = c(agreement_v2_sim, agreement_v3_sim), 
               names_to = "video", names_pattern = "agreement_v(.)_sim", 
               values_to = "agreement_sim") %>%
  mutate(video_fct = factor(video))
# scaling was already done in data prep, means are already 0

###################################################################################################################
# STEP 2: Single-video models
###################################################################################################################

vid.2.model <- run_lmer_dyad(agreement_sim ~ ideology_similarity * joint_IUS, agreement_sim_df %>% filter(video == 2),
                                 dyad_efs = ~(1 | dyad),
                                 dyad_video = FALSE,
                                 lmerTest = TRUE)
summary(vid.2.model)

vid.3.model <- run_lmer_dyad(agreement_sim ~ ideology_similarity * joint_IUS, agreement_sim_df %>% filter(video == 3),
                             dyad_efs = ~(1 | dyad),
                             dyad_video = FALSE,
                             lmerTest = TRUE)
summary(vid.3.model)

###################################################################################################################
# STEP 3: Across-video model
###################################################################################################################

cross.vid.model <- run_lmer_dyad(agreement_sim ~ ideology_similarity * joint_IUS * video_fct, agreement_sim_df,
  dyad_efs = ~(1 | dyad),
  dyad_video = TRUE,
  lmerTest = TRUE)

summary(cross.vid.model)

