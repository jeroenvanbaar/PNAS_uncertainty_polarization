library(tidyverse)
library(lme4)
library(lmerTest)
library(vroom)
library(plotly)
library(broom.mixed)
library(sjPlot)
library(gridExtra)

# GENERAL SETTINGS
base_dir <- 'BASE_DIR_HERE'

#####################################################################
# Run regression again on mean ROI data (across voxels)
#####################################################################

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

## DEFINE MODELS/TERMS
run <- 3
model <- 'ideology_IUS'
term <- 'scale(ideology_similarity)-X-joint_IUS'
iv_x <- 'ideology_similarity'
iv_y <- 'joint_IUS'
iv_xy <- sprintf('%s:%s',iv_x,iv_y)
var_to_plot <- 'mean_BOLD'

################################################################
# GRID METHOD (fig 3B)
################################################################

## DEFINE PLOT TYPE
plot_type <- 'surface'
ncol_surface <- 100
color_palette <- viridis::viridis_pal()
azim <- -40
colat <- 20
npoints <- 50
quartzFonts(arial = c("Arial","Arial Black","Arial Narrow", "Arial Rounded MT"))
par(family = 'arial')

# Load data and scale where appropriate
reg_data <- vroom(sprintf('%s/Results/voxelwise_ISRSA/nifti/run-%i_model-%s/ROI_regression_data_term-%s_%s.csv',
                          base_dir, run, model, term, var_to_plot))
relabeled <- transform(reg_data,
                       Pair_1ID = pmin(SubID1, SubID2),
                       Pair_2ID = pmax(SubID1, SubID2))
reg_data <- subset(relabeled, !duplicated(cbind(Pair_1ID, Pair_2ID)) & (SubID1 != SubID2))
reg_data <- reg_data %>%
  dplyr::select(-Pair_1ID, -Pair_2ID)
roi_cols <- names(reg_data)[c(10:(length(names(reg_data))-1))]
roi_select <- c('ISC.ROI.2','ISC.ROI.4','ISC.ROI.6')
roi_names <- c('lTPJ','rOFC','Precuneus')
roi_select

## BUILD FULL REGRESSION MODEL, THEN SIMULATE MARGINAL EFFECTS
print_plot <- TRUE
if (print_plot){
  plot_fname <- sprintf('%s/Results/voxelwise_ISRSA/nifti/run-%i_model-%s/ROI_simulations_%s_select-clusters.pdf',
                        base_dir, run, model, term)
  pdf(plot_fname)
  par(mar=c(2,8,6,6)+.1)
}
z_lower <- -.08
z_upper <- .14
for (dvi in 1:length(roi_select)){
  dv <- roi_select[dvi]
  
  ## RUN DYADIC REGRESSION
  reg_formula <- as.formula(sprintf(
                '%s ~ %s*%s',
                dv, iv_x, iv_y))
  
  res <- run_lmer_dyad(reg_data, reg_formula, lmerTest = TRUE)
  
  ## Simulate data
  # Create IV grid
  iv_x_list <- round(do.call(seq,as.list(c(round(range(reg_data[,iv_x]),2),round(diff(range(reg_data[,iv_x]))/npoints,2)))),2)
  iv_y_list <- round(do.call(seq,as.list(c(round(range(reg_data[,iv_y]),2),round(diff(range(reg_data[,iv_y]))/npoints,2)))),2)
  sim_dat <- as_tibble(expand.grid(iv_x = iv_x_list, iv_y = iv_y_list)) %>%
    mutate(iv_xy = iv_x * iv_y) %>% arrange()
  # Load marginal effects
  beta_int <- fixef(res)['(Intercept)']
  beta_x <- fixef(res)[iv_x][[1]]
  beta_y <- fixef(res)[iv_y][[1]]
  beta_xy <- fixef(res)[iv_xy][[1]]
  # Simulate DV
  ISC <- c()
  for (i in 1:dim(sim_dat)[1]){
    ISC <- c(ISC, c(beta_int
                    + sim_dat$iv_x[i] * beta_x 
                    + sim_dat$iv_y[i] * beta_y 
                    + sim_dat$iv_xy[i] * beta_xy))
  }
  sim_dat$pred_ISC <- ISC
  
  ## PLOT
  if (plot_type == 'scatter'){
    fig <- plot_ly(x=sim_dat$iv_x, y=sim_dat$iv_y, z=sim_dat$pred_ISC,
                   type="scatter3d", mode="markers", color=sim_dat$pred_ISC) %>%
                   layout(scene = list(xaxis = list(title = 'Ideology sim z'),
                                       yaxis = list(title = 'Joint IUS'),
                                       zaxis = list(title = 'ISC z')))
  } else if (plot_type == 'surface'){
    # Compute vertex values and colors
    nrz <- length(iv_x_list)
    ncz <- length(iv_y_list)
    z <- matrix(sim_dat$pred_ISC,nrow=nrz)
    # Start the color scale at the z_min of the current ROI, end at the z_max of the current ROI
    z_range <- range(z)
    col_min <- (z_range[1] - z_lower)/(z_upper - z_lower)
    col_max <- 1 - (z_upper - z_range[2])/(z_upper - z_lower)
    tmp_pal <- viridis::viridis_pal(begin = col_min, end = col_max, option = 'D')
    color <- tmp_pal(ncol_surface)
    zfacet1a <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]) / 4
    facetcol1a <- cut(zfacet1a, ncol_surface)
    # Graph
    mult <- 20
    zlim <- c(floor(min(z*mult))/mult,ceiling(max(z*mult))/mult)
    surf_fig <- persp(iv_x_list, iv_y_list, z, cex.axis = 1.5, cex.lab = 2.5,
                      ticktype="detailed", lwd = 1, border = NA,
                      col=color[facetcol1a], theta=azim, phi = colat, zlim = zlim, axes = T,
                      xlab="\n\n  Ideology\n   similarity", ylab="\n  Joint IUS",
                      zlab="\n\nNeural\nsynchrony", main = sprintf('%s',roi_names[dvi]))
  } else if (plot_type == 'heatmap'){
    # Heatmap
    myplot <- ggplot(sim_dat, aes(iv_x, iv_y)) +
      geom_raster(aes(fill = pred_ISC)) +
      scale_fill_gradientn(colours = viridis::viridis(100)) +
      labs(title = sprintf('ROI %s',strsplit(dv,'ROI.')[[1]][2])) +
      xlab("Ideology similarity (Z)") + ylab('Joint IUS') +
      theme_minimal()
    if (print_plot){
      ggsave(sprintf('%s/Results/voxelwise_ISRSA/nifti/run-%i_model-%s/ROI_simulations_%s_hm_ROI-%s.pdf',
                                           base_dir, run, model, term, strsplit(dv,'ROI.')[[1]][2]), plot = myplot)
    } else{
      show(myplot)
    }
  }
}
if (print_plot){
  dev.off()
}

################################################################
# SJPLOT METHOD (fig 3C)
################################################################

reg_data <- vroom(sprintf('%s/Results/voxelwise_ISRSA/nifti/run-%i_model-%s/ROI_regression_data_term-%s_%s.csv',
                          base_dir, run, model, term, var_to_plot))
relabeled <- transform(reg_data,
                       Pair_1ID = pmin(SubID1, SubID2),
                       Pair_2ID = pmax(SubID1, SubID2))
reg_data <- subset(relabeled, !duplicated(cbind(Pair_1ID, Pair_2ID)) & (SubID1 != SubID2))
reg_data <- reg_data %>%
  # mutate_at(c(3,5,6,10:length(names(reg_data))-3),scale) %>%
  dplyr::select(-Pair_1ID, -Pair_2ID)
roi_cols <- names(reg_data)[c(10:(length(names(reg_data))-1))]
dvs <- roi_cols[c(2,4,6)]
roi_names <- c('TPJ (L)','OFC (R)','Precuneus (L)')

# Build plot
plots <- list()
percentiles_to_plot <- c(.25,.75)
quanta <- reg_data %>% pull(joint_IUS) %>% quantile(probs = percentiles_to_plot)
print(quanta)
marginal_values_string <- paste('joint_IUS [',do.call(sprintf, c(list("%f, %f"), quanta)),']', sep = '')
print(marginal_values_string)
# colors <- c('#8f8f8f', '#696969', '#464646', '#252525', '#000000')
colors <- c('#C8C8C8', '#000000')
# labels <- c('10%','1st quartile','median','3rd quartile','90%')
labels <- c('Q1', 'Q3')
reg.table <- tibble()
for (i in 1:length(dvs)){
  reg_formula <- as.formula(sprintf('%s ~ %s*%s',dvs[[i]], iv_x, iv_y))
  res <- run_lmer_dyad(reg_data, reg_formula, lmerTest = TRUE)
  plots[[i]] <- sjPlot::plot_model(res, type = 'emm', terms = c('ideology_similarity',marginal_values_string),
                                   se = FALSE, ci.lvl = .95, show.data = FALSE, show.legend = TRUE,
                                   colors = colors, legend.title = 'Joint IUS')
  reg.table = rbind(reg.table,tidy(res) %>% mutate(ROI = roi_names[[i]]))
}
# Plot
set_theme(base = theme_classic(), theme.font = 'Arial', title.align = 'center', legend.pos = 'right')
g <- grid.arrange(
  plots[[1]] + ggtitle(roi_names[[1]]) + labs(x = 'Ideology similarity (dyad)', y = "Neural synchrony") +
    scale_color_manual(values = colors, name = 'Joint IUS', labels = labels) +
    scale_fill_manual(values = colors, name = 'Joint IUS', labels = labels),
  plots[[2]] + ggtitle(roi_names[[2]]) + labs(x = 'Ideology similarity (dyad)', y = "Neural synchrony") +
    scale_color_manual(values = colors, name = 'Joint IUS', labels = labels) +
    scale_fill_manual(values = colors, name = 'Joint IUS', labels = labels),
  plots[[3]] + ggtitle(roi_names[[3]]) + labs(x = 'Ideology similarity (dyad)', y = "Neural synchrony") +
    scale_color_manual(values = colors, name = 'Joint IUS', labels = labels) +
    scale_fill_manual(values = colors, name = 'Joint IUS', labels = labels),
  nrow=length(dvs))
ggsave(sprintf('%s/Results/voxelwise_ISC/nifti/run-%i_model-%s/ROI_fixef_marginals_%s_select-clusters.png',
               base_dir, run, model, term), plot = g, width = 70, units = "mm", height = 155, dpi = 300)

################################################################
# REGRESSION TABLE
################################################################
reg.table <- tibble()
table.dat <- reg_data %>%
  select(-age_distance, -scan_day_distance, -same_gender, -same_undergrad, -same_community) %>%
  pivot_longer(!c('SubID1','SubID2','ideology_similarity','joint_IUS','grouping'), names_to = 'ROI', values_to = 'synchrony')
table.dat
reg.table <- tibble()
for (i in 1:length(roi_select)){
  dat <- table.dat %>% filter(ROI == roi_select[[i]])
  res <- run_lmer_dyad(dat, as.formula('scale(synchrony) ~ scale(ideology_similarity) * scale(joint_IUS)'), lmerTest = TRUE)
  reg.table = rbind(reg.table,tidy(res) %>% mutate(ROI = roi_names[[i]]))
}
reg.table <- reg.table %>%
  select(ROI, term, estimate, std.error, statistic, df, p.value) %>%
  drop_na()
View(reg.table)

