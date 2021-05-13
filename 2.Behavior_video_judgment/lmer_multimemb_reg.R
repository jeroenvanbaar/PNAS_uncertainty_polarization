library(tidyverse)
library(lme4)
library(lmerTest)
library(modelr)

get_dyad_memb_mat <- function(formula, data, drop_Inf = FALSE) {
  dyad_1 = quo(SubID1)
  dyad_2 = quo(SubID2)
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
  
  return(list(dyad_mat=as.matrix(dyad_mat), data=model_df_na_omit))
}

#https://bbolker.github.io/mixedmodels-misc/notes/multimember.html
lmer_multimemb <- function(formula, data, memb_mat=list(), lmerTest = FALSE) {
  ## FIXME: pass ... through appropriately
  ## FIXME: test dimensions
  mnms <- names(memb_mat)
  #random effects terms
  fb <- findbars(formula)
  #names of grouping variables in formula
  gvars <- vapply(fb, function(x) deparse(x[[3]]), character(1))
  Ztlist <- list()
  for (i in seq_along(fb)) {
    fbnm <- deparse(fb[[i]])
    m <- gvars[i]
    ## find corresponding random-effects term
    w <- which(mnms==gvars[i])
    if (length(w)>0) {
      M <- Matrix::Matrix(memb_mat[[w]])
      ## extract LHS (effect)
      form <- as.formula(substitute(~z, list(z=fb[[i]][[2]])))
      ## construct model matrix & compute Khatri-Rao product 
      X <- model.matrix(form, data = data)
      Zt <- Matrix::KhatriRao(t(M), t(X), make.dimnames=TRUE)
      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt
      ## if necessary, add factor to data
      if (!m %in% names(data)) {
        ## if the factor has non-trivial ordering, it should be included
        ## in the data.  Do we have to worry about ordering of Z? test!
        data[[gvars[i]]] <- rep(factor(colnames(memb_mat[[w]])), length.out=dim(data)[1])
      }
    } ## if  (length(w)>0)
  } ## for i in seq(fb)
  ## call lFormula  (FIXME: allow glFormula)
  lmod <- lFormula(formula, data = data)
  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
  ## finish fitting
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun)
  res <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  if(lmerTest) {
    res <- lmerTest:::as_lmerModLT(res, devfun)
  }
  return(res)
}

# run_lmer_dyad function --------------------------------------------------
run_lmer_dyad <- function(formula, data, dyad_efs = ~(1 | dyad), dyad_video = FALSE, lmerTest = TRUE) {
  dyad_list <- get_dyad_memb_mat(formula, data)
  dyad_mat <- dyad_list$dyad_mat
  formula <- add_predictors(formula, dyad_efs)
  memb_mat <- list(dyad = dyad_mat)
  if(dyad_video) {
    X <- model.matrix(~0 + video_fct, data = data)
    dyad_video_mat <- t(Matrix::KhatriRao(t(dyad_mat), t(X), make.dimnames=TRUE))
    memb_mat <- list(dyad = dyad_mat, dyad_video = dyad_video_mat)
    formula <- add_predictors(formula, ~(1 | dyad_video))
  }
  m <- lmer_multimemb(formula, dyad_list$data, 
                      memb_mat = memb_mat, lmerTest = lmerTest)
  m
}
