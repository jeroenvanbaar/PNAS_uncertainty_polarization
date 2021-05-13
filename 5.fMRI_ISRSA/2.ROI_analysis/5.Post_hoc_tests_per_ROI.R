library(lme4)
library(lmerTest)
library(vroom)
library(tidyverse)
library(sjPlot)
library(reshape2)
library(broom.mixed)
library(emmeans)

# GENERAL SETTINGS
base_dir <- 'BASE_DIR_HERE'
use.anova.output <- TRUE
if (use.anova.output){
  ROI.data.dir <- sprintf('%s/Results/voxelwise_ISRSA/nifti_with_anova',base_dir)
} else {
  ROI.data.dir <- sprintf('%s/Results/voxelwise_ISRSA/nifti',base_dir)
}

# FUNCTION DEFINITIONS

run_lmer_dyad <- function(reg_data, reg_formula, lmerTest = FALSE) {
  # Define terms
  dyad_1 = quo(SubID1)
  dyad_2 = quo(SubID2)
  mod_terms <- all.vars(subbars(reg_formula))
  hier_formula <- update(reg_formula, ~ . + (1 | dyad))
  
  #model only uses rows with no NAs so need to select only members of dyads in those rows. Also drop inf
  model_df_na_omit <- reg_data %>% 
    dplyr::select(tidyselect::all_of(mod_terms), {{dyad_1}}, {{dyad_2}}) %>% 
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
# STEP 1: Load data
###################################################################################################################

## DEFINE MODELS/TERMS
run <- 3
model <- 'party_IUS'
term <- 'ideology_pair-X-joint_IUS'
var_to_plot <- 'mean_BOLD'

# Load predictor & brain data
predictor.rdms <- as_tibble(vroom(sprintf('%s/Data/Cleaned/Surveys/predictor_RDMs_4.csv',base_dir)))
roi.data <-
  tibble(vroom(sprintf('%s/run-%i_model-%s/ROI_regression_data_term-%s_%s.csv',
                       ROI.data.dir, run, model, term, var_to_plot))) %>%
  dplyr::select(-ideology_similarity, -(age_distance:same_community)) %>%
  left_join(predictor.rdms %>% dplyr::select(SubID1,SubID2,ideology_pair)) %>%
  melt(id.vars = c('SubID1','SubID2','ideology_pair','joint_IUS'), variable.name = 'ROI', value.name = 'synchrony')

roi.names <- roi.data %>% pull(ROI) %>% unique()
print(roi.names)

###################################################################################################################
# STEP 2: Run model in each ROI => post-hoc tests
###################################################################################################################

# Dplyr/tidy way
reg.formula <- as.formula('synchrony ~ ideology_pair * joint_IUS')
reg.stats <-
  roi.data %>% 
  # filter(ROI == "ISC.ROI.14") %>%
  group_by(ROI) %>%
  nest() %>%
  mutate(dyad.reg.stats = purrr::map(data, run_lmer_dyad, reg.formula, lmerTest = TRUE),
         dyad.reg.augment = purrr::map(dyad.reg.stats, broom.mixed::augment),
         dyad.reg.emt = purrr::map(dyad.reg.stats, emtrends, "ideology_pair", var = "joint_IUS"),
         dyad.reg.contrast = purrr::map(dyad.reg.emt, contrast, "trt.vs.ctrl"),
         dyad.reg.stats.tidied = purrr::map(dyad.reg.stats, broom::tidy),
         dyad.reg.contrast.tidied = purrr::map(dyad.reg.contrast, broom::tidy))

reg.predictions <- reg.stats %>%
  unnest(dyad.reg.augment) %>%
  dplyr::select(-data, -dyad.reg.stats)

reg.contrast.out <- reg.stats %>%
  unnest(dyad.reg.contrast.tidied) %>%
  dplyr::select(-data, -dyad.reg.emt, -dyad.reg.contrast, -dyad.reg.stats, -dyad.reg.stats.tidied)

reg.out <- reg.stats %>%
  unnest(dyad.reg.stats.tidied) %>%
  dplyr::select(-data, -dyad.reg.emt, -dyad.reg.contrast, -dyad.reg.stats, -dyad.reg.contrast.tidied, -group)

reg.out %>%
  write_csv(sprintf('%s/run-%i_model-%s/ROI_post_hoc_stats_term-%s.csv',
                    ROI.data.dir, run, model, term))

# Both groups interaction true using dunnett adjustment?
reg.contrast.out %>% 
  filter(adj.p.value < 0.05) %>% 
  count(ROI) %>% mutate(
    both.sig = n == 2
  )

# Same thing but looking at the coeffs
reg.contrast.posandsig <-
  reg.contrast.out %>% 
  group_by(ROI) %>% mutate(
    # n = n(),
    both.sig = sum(adj.p.value < 0.05) == 2,
    both.pos = sum(estimate > 0) == 2,
    both.pos.sig = both.pos & both.sig
  ) %>%
  select(ROI,contrast,estimate, both.pos.sig)
significant.ROIs <- as.character(reg.contrast.posandsig %>% select(ROI,both.pos.sig) %>% filter(both.pos.sig) %>% unique() %>% pull(ROI))
print(significant.ROIs)

# Both groups interaction true?
both.groups.sig.out <-
  reg.out %>%
  filter(term == 'ideology_pairWithin_con:joint_IUS' | term == 'ideology_pairWithin_lib:joint_IUS') %>%
  reshape2::dcast(formula =  ROI ~ term, value.var = 'p.value') %>%
  rename_with(~c('ROI','con.interact','lib.interact')) %>%
  mutate(both.sig = con.interact < 0.05 & lib.interact < 0.05)

# Keep estimates
reg.out %>%
  filter(term == 'ideology_pairWithin_con:joint_IUS' | term == 'ideology_pairWithin_lib:joint_IUS') %>%
  group_by(ROI) %>%
  mutate(sum.sig = sum(p.value < 0.05)) %>% 
  filter(sum.sig == 2) %>%
  reshape2::dcast(formula =  ROI ~ term, value.var = 'estimate') %>%
  rename_with(~c('ROI','con.estimate','lib.estimate'))

both.groups.sig.out
View(both.groups.sig.out)

###################################################################################################################
# STEP 3: Plot raw data with fixed-effects regression lines
###################################################################################################################

plot.data <-
  reg.predictions %>%
  left_join(reg.contrast.posandsig %>% select(ROI,both.pos.sig)) %>%
  filter(both.pos.sig == TRUE)
View(plot.data)

ggplot(data = plot.data) +
  geom_point(aes(x = joint_IUS, y = synchrony, color = ideology_pair), alpha = .2) +
  geom_smooth(aes(x = joint_IUS, y = .fixed, color = ideology_pair), method = 'lm', formula = y ~ x, se = T) +
  # geom_line() +
  facet_wrap(~ROI) +
  scale_color_manual(values=c("purple", "red", "blue"))

select.ROIs <- c("ISC.ROI.10","ISC.ROI.14")
ROI.labels <- c("rOFC","rTPJ")
names(ROI.labels) <- select.ROIs
dyad.labels <- c("CL","CC","LL")
names(dyad.labels) <- c("Between","Within_con","Within_lib")
fname.out <- sprintf('%s/run-%i_model-%s/ROI_raw_and_fixef_plot_term-%s_%s.png',
                     ROI.data.dir, run, model, term, var_to_plot)
print(paste('Saving to ',fname.out))
ggplot(data = plot.data %>% filter(ROI %in% select.ROIs)) +
  geom_point(aes(x = joint_IUS, y = synchrony, color = ideology_pair), alpha = .05) +
  geom_smooth(aes(x = joint_IUS, y = .fixed, color = ideology_pair), method = 'lm', formula = y ~ x) +
  # geom_line() +
  facet_wrap(~ROI, labeller = labeller(ROI = ROI.labels)) +
  scale_color_manual(values=c("purple", "red", "blue"), labels = dyad.labels) +
  labs(color = "Dyad type:") +
  xlab('Joint IUS') + ylab('Neural synchrony') +
  theme_classic() +
  theme(legend.position = 'top',
        strip.background = element_rect(size=0),
        text=element_text(size=16,  family="Arial")) +
  ggsave(fname.out)

# SJPLOT METHOD
plots <- list()
colors <- c('purple','red','blue')
for (i in 1:length(select.ROIs)){
  myroi <- select.ROIs[[i]]
  model <- reg.stats %>% filter(ROI == myroi) %>% pull(dyad.reg.stats)
  plots[[i]] <- sjPlot::plot_model(model[[1]], type = 'emm', terms = c('joint_IUS','ideology_pair'),
                 ci.lvl = .95, SE = FALSE,
                 show.data = TRUE, colors = colors, dot.size = 1, show.legend = TRUE)
}
set_theme(base = theme_classic(), theme.font = 'Arial', axis.ticksize.y = 0,
          title.align = 'center', legend.item.backcol = 'white', legend.pos = 'top', legend.item.bordercol = NA)
library(gridExtra)
labels <- c('CL','CC','LL')
g <- grid.arrange(nrow = 1, 
  plots[[1]] + ylim(-.2,.35) + ggtitle(ROI.labels[[1]]) +
    labs(y = "Neural synchrony", x = "Joint IUS")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14), plot.title = element_text(size = 20, face = "bold")) +
    scale_color_manual(values = colors, name = 'Dyad type', labels = labels) +
    scale_fill_manual(values = colors, name = 'Dyad type', labels = labels),
  plots[[2]] + ylim(-.2,.35) + ggtitle(ROI.labels[[2]]) +
       labs(x = "Joint IUS") + 
       theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
             axis.title=element_text(size=14), plot.title = element_text(size = 20, face = "bold")) +
    scale_color_manual(values = colors, name = 'Dyad type', labels = labels) +
    scale_fill_manual(values = colors, name = 'Dyad type', labels = labels))

ggsave(plot = g, filename = fname.out)

## Formatted regression tables
# Re-run with scaled joint IUS and scaled synchrony
reg.formula <- as.formula('scale(synchrony) ~ ideology_pair * scale(joint_IUS)')
reg.stats <-
  roi.data %>% 
  filter(ROI %in% select.ROIs) %>%
  group_by(ROI) %>%
  nest() %>%
  mutate(dyad.reg.stats = purrr::map(data, run_lmer_dyad, reg.formula, lmerTest = TRUE))

reg.table <- tibble()
for (i in 1:length(select.ROIs)){
  myroi <- select.ROIs[[i]]
  a <- reg.stats %>% filter(ROI == myroi) %>% pull(dyad.reg.stats)
  a <- tidy(a[[1]]) %>%
    mutate(ROI = ROI.labels[[i]])
  reg.table = rbind(reg.table,a)
}
reg.table <- reg.table %>%
  select(ROI, term, estimate, std.error, statistic, df, p.value) %>%
  drop_na()
View(reg.table)
