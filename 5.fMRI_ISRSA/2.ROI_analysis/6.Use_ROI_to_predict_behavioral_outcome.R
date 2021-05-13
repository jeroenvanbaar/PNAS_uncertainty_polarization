library(tidyverse)
library(mediation)
library(lme4)
library(lmerTest)
library(vroom)
library(plotly)
library(sjPlot)
library(broom.mixed)
library(emmeans)
library(MuMIn)

# GENERAL SETTINGS
base_dir <- 'BASE_DIR_HERE'

# FUNCTION DEFINITIONS

run_lmer_dyad <- function(reg_data, reg_formula, lmerTest = FALSE) {
  # Define terms
  dyad_1 = quo(SubID1)
  dyad_2 = quo(SubID2)
  mod_terms <- all.vars(subbars(reg_formula))
  hier_formula <- update(reg_formula, ~ . + (1 | dyad))
  
  #model only uses rows with no NAs so need to select only members of dyads in those rows. Also drop inf
  model_df_na_omit <- reg_data %>% dplyr::select(tidyselect::all_of(mod_terms), {{dyad_1}}, {{dyad_2}}) %>% 
    drop_na()
  
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
  
  # Fit model
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


###################################################################################################################
# 1. Load data
###################################################################################################################

## DEFINE MODELS/TERMS
run <- 3
model <- 'ideology_IUS'
term <- 'scale(ideology_similarity)-X-joint_IUS'
use.anova.output <- FALSE
if(use.anova.output){
  ROI.data.dir <- sprintf('%s/Results/voxelwise_ISRSA/nifti_with_anova',base_dir)
} else {
  ROI.data.dir <- sprintf('%s/Results/voxelwise_ISRSA/nifti',base_dir)
}
print(sprintf('Using data dir %s',ROI.data.dir))

# Load brain data
var_to_plot <- 'mean_BOLD'
roi.data <- tibble(vroom(sprintf('%s/run-%i_model-%s/ROI_regression_data_term-%s_%s.csv',
                          ROI.data.dir, run, model, term, var_to_plot)))
roi.data <- roi.data %>%
  transform(Pair_1ID = pmin(SubID1, SubID2),
            Pair_2ID = pmax(SubID1, SubID2)) %>%
  subset(!duplicated(cbind(Pair_1ID, Pair_2ID)) & (SubID1 != SubID2)) %>%
  select(-Pair_1ID, -Pair_2ID) %>%
  tibble()
roi.cols <- names(roi.data)[10:(length(names(roi.data))-1)]
print(roi.cols)

# Load judgment data
judg.data <- vroom(sprintf('%s/Results/video_judgment/all_judgment_similarity.csv',base_dir))
DV.cols <- names(judg.data)[3:length(names(judg.data))]

# Merge brain and behavior data
model_data <-
  roi.data %>%
  left_join(judg.data) %>% 
  mutate_at(c(3,5,6,10:19,21:27),scale) %>% 
  mutate_at(c(3,5,6,10:19,21:27),funs(as.numeric(.))) %>% 
  mutate(ideo_sim_x_joint_IUS = ideology_similarity*joint_IUS)
print(dim(model_data))


###################################################################################################################
# 2. Add ROI data to model of judgment and test which has explanatory power.
###################################################################################################################

# Behavior-only model
DV.name <- "agreement_v3_sim"
beh.formula <- as.formula(sprintf('%s ~ ideology_similarity*joint_IUS',DV.name))
beh.model.fit <- run_lmer_dyad(model_data, beh.formula, lmerTest = TRUE)

# This shows results without any brain data:
summary(beh.model.fit)

# Behavior-and-brain models
res <- tibble()
roi.choose <- roi.cols
roi.choose <-sprintf('ISC.ROI.%i',c(2,4,6,7,10))
for (ri in 1:length(roi.choose)) {
  ROI.name <- roi.choose[[ri]]
  print(sprintf('%s',ROI.name))
  
  beh.brain.formula <- as.formula(sprintf('%s ~ ideology_similarity*joint_IUS + %s',DV.name,ROI.name))
  beh.brain.model.fit <- run_lmer_dyad(model_data, beh.brain.formula, lmerTest = TRUE)
  
  model.comp.p <-
    anova(beh.model.fit, beh.brain.model.fit) %>%
    tidy() %>% filter(term == 'beh.brain.model.fit') %>% pull(p.value)
  
  res <- res %>%
    rbind(tidy(beh.brain.model.fit) %>% mutate(ROI.term = term == ROI.name) %>% mutate(model.comp.p = model.comp.p))
}

res %>%
  filter(ROI.term) %>%
  group_by(term) %>%
  mutate(p.one.sided = if(estimate > 0) (p.value / 2) else nan)
  
# This analysis reveals that ROIs 4 (rAI) and possibly 6/10 (precuneus) can both contribute to the model of judgment similarity
