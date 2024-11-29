library(plsmod)
library(tidyverse)
library(tidymodels)
library(ggpubr)

#Local directory
wd <- 'C:/Users/p6exk/Box/UCDavis/'

#Read in and clean data---------------------------------------------------------

#read in spectral data
df_spectra_2021 <- read.csv(paste0(wd, "RMI_grapes/RMI2021/spectroscopy/processed_data/svc_clean_interp.csv"), check.names = F)
df_spectra_2023 <- read.csv(paste0(wd, "RMI_grapes/RMI2023/spectroscopy/svc_clean_interp.csv"), check.names = F)

#merge spectra
df_spectra <- rbind(df_spectra_2021, df_spectra_2023)

# extract leaf spectra and clean
df_spectra <- df_spectra %>%
  group_by(Block, Row, Vine, Variety, date) %>%
  filter(grepl(pattern ="berry", rep)) %>%
  select(!c(rep, ID, scan)) %>%
  filter(`1500` <= 25) %>% #screen out bad spectra
  summarise_all(mean, na.rm = T) %>%
  select(!(`967`:`1021`)) %>% #screen out hot pixels
  select(!c(`1340`:`1445`, `1790`:`1955`, `2401`:`2500`)) %>% #screen out water absorption pixels
  mutate(NDVI = (`750` - `670`) / (`750` + `670`),
         PRI = (`531` - `570`) / (`531` + `570`),
         TVI = 0.5 * (120 * (`750` - `550`) - 200 * (`670` - `550`))) %>%
  mutate(Row = parse_number(Row),
         Vine = parse_number(Vine),
         date = as.Date(date, format = "%m/%d/%Y"))

#read in origin info
df_origin <- read.csv(paste0(wd, "RMI_grapes/RMI2021/VarietyOrigins.csv"))

#read in berry composition data
df_berry <- read.csv(paste0(wd, "RMI_grapes/RMI2023/berry_chemistry/berrychem_clean_2021_23.csv"))

#clean berry df
df_berry <- df_berry %>%
  mutate(Row = parse_number(Row),
         Vine = parse_number(Vine),
         brix = as.numeric(brix),
         ta = as.numeric(ta),
         date = as.Date(date, format = "%m/%d/%Y"))

#Merging dfs ----------------------------------------

#merge spectra with origins
df_data <- merge(df_spectra, df_origin, by = "Variety")

#merge spectra with berry composition
df_data <- list(df_data, df_berry) %>%
  reduce(left_join, by = c("date", "Block", "Row", "Vine")) %>%
  mutate(date = as.Date(date))

#quick plots------------------------------------------------
#plot spectra
plt_spectra <- df_data %>%
  pivot_longer(cols = `400`:`2400`, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = interaction(Block, Row, Vine), color = Variety))+
  geom_line(show.legend = T)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0,100))+
  scale_color_viridis_d(option = "H")+
  facet_wrap(~date, nrow = 3)+
  theme_bw()

#plot berry composition
plt_phys <- df_data %>%
  pivot_longer(cols = c(NDVI, PRI, TVI, brix, ph, ta), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = as.factor(date), y = value, group = interaction(Variety, date), color = Variety))+
  geom_boxplot(show.legend = F)+
  scale_x_discrete("date")+
  scale_y_continuous("value")+
  facet_wrap(~param, scales = "free_y")+
  theme_bw()

#plot ta zoomed in
plt_ta <- df_data %>%
  filter(date>= as.Date("2023-07-01")) %>%
  ggplot(aes(date, ta, color=Variety, group = interaction(Block, Row, Vine)))+
  geom_point()+
  geom_line()

#plot VI to berry scatterplots
plt_reg <- df_data %>%
  #ggplot(aes(x = brix, y = ((`531`-`570`)/(`531`+`570`)), color = as.factor(date)))+    #PRI
  #ggplot(aes(x = brix, y = NDVI, color = as.factor(date)))+     #NDVI
  ggplot(aes(x = brix, y = TVI, color = as.factor(date)))+     #TVI
  geom_point()+
  stat_cor(aes(color = NULL))+
  theme_bw()

#pls data cleaning and splitting--------------------------------

#Keep only certain columns
df_pls <- df_data#[,-c(2:5, 6:8, 2111, 2114:2115)]
df_pls$date <- as.character(df_pls$date)

#data split into train and test
set.seed(234)
perm_split <- initial_split(df_pls, prop = 0.7, strata = date)
perm_train <- training(perm_split)
perm_test <- testing(perm_split)

#make only post veraison
perm_train_veraison <- perm_train %>%
  filter(!date == as.Date("2021-06-24"))
perm_test_veraison <- perm_test %>%
  filter(!date == as.Date("2021-06-24"))

#V-fold cross validation
set.seed(345)
perm_folds <-  vfold_cv(perm_train)
perm_folds_veraison <-  vfold_cv(perm_train_veraison)

#Preparing pls model ---------------------------

#recipe brix
pls_rec_brix <- recipe(brix ~ ., data = perm_train) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(brix, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

pls_rec_brix_veraison <- recipe(brix ~ ., data = perm_train_veraison) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(brix, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

#recipe ta
pls_rec_ta <- recipe(ta ~ ., data = perm_train) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(ta, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

pls_rec_ta_veraison <- recipe(ta ~ ., data = perm_train_veraison) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(ta, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

#recipe ph
pls_rec_ph <- recipe(ph ~ ., data = perm_train) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(ph, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

pls_rec_ph_veraison <- recipe(ph ~ ., data = perm_train_veraison) %>%
  update_role(NDVI, PRI, TVI, Variety, Block, Row, Vine, date, Region_of_Origin, variety, Variety_f, brix, ph, ta, new_role = 'id') %>%
  update_role(ph, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_naomit(all_outcomes())

#create empty data frames to save into
rec_list <- list(pls_rec_brix, pls_rec_ta, pls_rec_ph)
pred_list <- c("brix", "ta", "ph")
df_tune <- data.frame()
df_vip <- data.frame()
df_train <- data.frame()
df_test <- data.frame()

#PLSR----------------------------------------------
#loop to run PLSR over the recipes
for (i in 1:length(rec_list)){
  print(paste("In progress:", i, "out of", length(rec_list)))
  #identify recipe
  pls_rec <- rec_list[[i]]
  
  #pls tuning specifications
  pls_tuning <- parsnip::pls(num_comp = tune()) %>% 
    set_mode("regression") %>% 
    set_engine("mixOmics")
  
  #create tuning workflow
  pls_workflow <- workflow() %>% 
    add_recipe(pls_rec) %>% 
    add_model(pls_tuning)
  
  #enable parallel processing to speed up
  doParallel::registerDoParallel()
  
  #number of pc to test
  comp_grid <- tibble(num_comp = seq(from = 1, to = 20, by = 1))
  
  #choosing resample folds
  if (grepl("veraison", pred_list[[i]], fixed = TRUE)) {
    folds = perm_folds_veraison
  } else {
    folds = perm_folds
  }
  
  #pls model tuning
  tuned_pls_results <- pls_workflow %>% 
    tune_grid(resamples = folds,
              grid = comp_grid,
              metrics = metric_set(mae, rmse, rsq))
  
  #select best model
  tuned_best <- tuned_pls_results %>% 
    #select_best("rsq")    #absolute best based on either "mae", "rmse", or "rsq"
    select_by_pct_loss(    #selecting least complex model with no more than limit % loss of rmse
      metric = "rmse",
      limit = 3, 1
    )
  
  #extract tuning summary
  tune <- tuned_pls_results %>% 
    collect_metrics() %>%
    mutate(param = pred_list[i], 
           best_comp = tuned_best$num_comp)
  
  #output to dataframe
  df_tune <- rbind(df_tune, tune)
  
  #update model and workflow with tuned parameters
  updated_pls_model <-  parsnip::pls(num_comp = tuned_best$num_comp) %>% 
    set_mode("regression") %>% 
    set_engine("mixOmics")
  
  updated_workflow <- pls_workflow %>% 
    update_model(updated_pls_model)
  
  #create pls model
  if (grepl("veraison", pred_list[[i]], fixed = TRUE)) {
    pls_model <- updated_workflow %>% 
      fit(data = perm_train_veraison)
  } else {
    pls_model <- updated_workflow %>% 
      fit(data = perm_train)
  }
  

  
  #assessing pls model -
  #extract weights
  pls_weights <- pls_model %>% 
    extract_fit_parsnip() %>% 
    tidy()
  
  #extract weights
  weights <- pls_weights %>% 
    filter(term != "Y", # outcome variable col name
           component == tuned_best$num_comp) %>% 
    select(value) %>%
    rename(weight = value)
  
  #extract vip
  library(pls)
  var_list <- pls_rec$term_info %>% filter(role == "predictor" | role == "outcome") %>% select(variable)
  perm_test2 <- subset(perm_test, select = c(var_list$variable))
  plsr_mod <- plsr(as.formula(paste(gsub("_veraison", "", pred_list[i]), "~", ".")), ncomp = tuned_best$num_comp, data = perm_test2, validation = "LOO")
  #get regression coefficient
  regcoef <- as.data.frame(coef(plsr_mod, intercept = F, ncomp = tuned_best$num_comp))
  colnames(regcoef) <- 'regcoef'
  #calculate vip
  SS <- c(plsr_mod$Yloadings)^2 * colSums(plsr_mod$scores^2)
  Wnorm2 <- colSums(plsr_mod$loading.weights^2)
  SSW <- sweep(plsr_mod$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
  vip <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
  
  vip <- as_tibble(vip[tuned_best$num_comp,], rownames = 'wl') %>%
    slice(-n()) %>%
    unnest(-wl) %>% 
    mutate(wl = parse_number(wl),
           regcoef%>%
             slice(-n()),
           weights,
           ncomp = tuned_best$num_comp,
           param = pred_list[i]) %>%
    rename(vip = value)
  
  #output vip to dataframe
  df_vip <- rbind(df_vip, vip)
  
  #extract pls training predictions
  if (grepl("veraison", pred_list[[i]], fixed = TRUE)) {
    pls_train <- pls_model %>% 
      predict(new_data = perm_train_veraison) %>% 
      mutate(truth = perm_train_veraison[,gsub("_veraison", "", pred_list[i])], date = perm_train_veraison$date, Variety = perm_train_veraison$Variety, Variety_f = perm_train_veraison$Variety_f, Vine = perm_train_veraison$Vine, Row = perm_train_veraison$Row, Block = perm_train_veraison$Block, param = pred_list[i], model = "train", ndvi = perm_train_veraison$NDVI, pri = perm_train_veraison$PRI) 
  } else {
    pls_train <- pls_model %>% 
      predict(new_data = perm_train) %>% 
      mutate(truth = perm_train[,gsub("_veraison", "", pred_list[i])], date = perm_train$date, Variety = perm_train$Variety, Variety_f = perm_train$Variety_f, Vine = perm_train$Vine, Row = perm_train$Row, Block = perm_train$Block, param = pred_list[i], model = "train", ndvi = perm_train$NDVI, pri = perm_train$PRI) 
  }
  
  #output training results to dataframe
  df_train <- rbind(df_train, pls_train)
  
  #extract pls testing predictions
  if (grepl("veraison", pred_list[[i]], fixed = TRUE)) {
    pls_test <- pls_model %>% 
      predict(new_data = perm_test_veraison) %>% 
      mutate(truth = perm_test_veraison[,gsub("_veraison", "", pred_list[i])], date = perm_test_veraison$date, Variety = perm_test_veraison$Variety, Variety_f = perm_test_veraison$Variety_f, Vine = perm_test_veraison$Vine, Row = perm_test_veraison$Row, Block = perm_test_veraison$Block, param = pred_list[i], model = "test", ndvi = perm_test_veraison$NDVI, pri = perm_test_veraison$PRI) 
    
  } else {
    pls_test <- pls_model %>% 
      predict(new_data = perm_test) %>% 
      mutate(truth = perm_test[,gsub("_veraison", "", pred_list[i])], date = perm_test$date, Variety = perm_test$Variety, Variety_f = perm_test$Variety_f, Vine = perm_test$Vine, Row = perm_test$Row, Block = perm_test$Block, param = pred_list[i], model = "test", ndvi = perm_test$NDVI, pri = perm_test$PRI) 
  }  
  
  #output testing data results to dataframe
  df_test <- rbind(df_test, pls_test)

}


#plots--------------------------------------------------

#plot tuning summary
df_tune %>% 
  ggplot(aes(num_comp, mean, col = param)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(n.breaks = 20) +
  labs(x = "Number of components",
       y = "Indicator",
       title = "Plot of MAE, RMSE and R-SQ vs number of components for TRAINING dataset, with 10-fold repeated cross validation") +
  facet_grid(.metric ~., scales = "free_y") +
  theme_bw()


#plot vip
df_vip %>% 
  ggplot(aes(x = wl, y = vip, color = interaction(param, ncomp)))+
  geom_line(show.legend = T, linewidth = 1)+
  scale_x_continuous("wavelength (nm)")+
  scale_y_continuous("Variable Importance in Projection (VIP)", limits = c(0, 4))+
  scale_color_manual(values = c("#264653", "#2a9d8f", "#8ab17d", "#e76f51", "#f4a261", "#e9c46a"))+
  theme_bw()


#plot scatterplot training data
df_train %>% 
  ggplot(aes(truth, .pred, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), method = 'lm', formula = y~x) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  labs(title = "Actual vs Predicted for TRAINING dataset",
       x = "Observed value",
       y = "Predicted value") +
  stat_cor(aes(color = NULL, label = paste(rr.label, ..p.label.., sep = "~`,`~")))+
  theme_bw()+theme(aspect.ratio = 1)+
  facet_wrap(~param, scales = "free")

#plot scatterplot testing data
df_test %>% 
  ggplot(aes(truth, .pred, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), method = 'lm', formula = y~x) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  labs(title = "Actual vs Predicted for TESTING dataset",
       x = "Observed value",
       y = "Predicted value") +
  stat_cor(aes(color = NULL, label = paste(rr.label, ..p.label.., sep = "~`,`~")))+
  theme_bw()+theme(aspect.ratio = 1)+
  facet_wrap(~param, scales = "free")

#PLSR spectral strength results------------------------------------------------ 
#VIPS
plt_vip <- df_vip %>% 
  ggplot(aes(x = wl, y = vip, color = interaction(param, ncomp)))+
  geom_line(show.legend = T, linewidth = 1)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous(expression(atop("Variable Importance", "in Projection (VIP)")), limits = c(0, 2.6))+
  scale_color_manual(NULL,values = c("#264653", "#8ab17d", "#f4a261"), labels = c("Tartaric Acid", "Brix", "pH"))+
  theme_bw()

#Regression coefficients
plt_reg <- df_vip %>% 
  ggplot(aes(x = wl, y = regcoef, color = interaction(param, ncomp)))+
  geom_line(show.legend = T, linewidth = 1)+
  scale_x_continuous("wavelength (nm)")+
  scale_y_continuous(expression(atop("Regression", "coefficient")), limits = c(-.05, .05))+
  scale_color_manual(NULL,values = c("#264653", "#8ab17d", "#f4a261"), labels = c("Tartaric Acid", "Brix", "pH"))+
  theme_bw()

#Model weights
plt_weight <- df_vip %>% 
  ggplot(aes(x = wl, y = weight, color = interaction(param, ncomp)))+
  geom_line(show.legend = T, linewidth = 1)+
  scale_x_continuous("wavelength (nm)")+
  scale_y_continuous("Model weighing", limits = c(-.1, .1))+
  scale_color_manual(NULL,values = c("#264653", "#8ab17d", "#f4a261"), labels = c("Tartaric Acid", "Brix", "pH"))+
  theme_bw()

#Saving spec figure
png("figures/fig_vip_berry.png", units="in", height = 7, width = 8, res=500)

ggarrange(plt_weight, 
          plt_reg,
          plt_vip,
          ncol = 1, nrow = 3, align = "hv", common.legend = T, legend = "right",
          labels = c("(a)", "(b)", "(c)", "(d)"))

dev.off()   # Stop writing to the file

#scatter plots-----------------------------------------------

#get stats
df_stats <- df_test %>%
  group_by(param) %>%
  na.omit() %>%
  reframe(rsq = cor.test(x = .pred, y = truth)$estimate^2,
          p = cor.test(x = .pred, y = truth)$p.value,
          mse = median((.pred - truth)^2),
          rmse = sqrt(mse),
          rmspe = sqrt(median(((truth - .pred) / truth)^2))*100,
          mae = mean(abs(.pred - truth)))

#get stats NDVI
df_stats_ndvi <- df_test %>%
  group_by(param) %>%
  na.omit() %>%
  reframe(rsq = cor.test(x = ndvi, y = truth)$estimate^2,
          p = cor.test(x = ndvi, y = truth)$p.value,
          mse = median((ndvi - truth)^2),
          rmse = sqrt(mse),
          rmspe = sqrt(median(((truth - ndvi) / truth)^2))*100,
          mae = mean(abs(ndvi - truth)))
#get stats PRI
df_stats_pri <- df_test %>%
  group_by(param) %>%
  na.omit() %>%
  reframe(rsq = cor.test(x = pri, y = truth)$estimate^2,
          p = cor.test(x = pri, y = truth)$p.value,
          mse = median((pri - truth)^2),
          rmse = sqrt(mse),
          rmspe = sqrt(median(((truth - pri) / truth)^2))*100,
          mae = mean(abs(pri - truth)))

col_all <- c('#e03524', '#f07c12', '#ffc200', '#90bc1a', '#21b534', '#0095ac', '#1f64ad', '#4040a0', '#903498')
col_ver <- c('#f07c12', '#ffc200', '#90bc1a', '#21b534', '#0095ac', '#1f64ad', '#4040a0', '#903498')

plt_sp_brix_plsr <- df_test %>% 
  filter(param == "brix") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous("Predicted Brix", limits = c(0,35))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_all)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "brix", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "brix", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "brix", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_brix_plsr_v <- df_test %>% 
  filter(param == "brix_veraison") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous("Predicted Brix", limits = c(0,35))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "brix_veraison", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "brix_veraison", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "brix_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_brix_ndvi <- df_test %>% 
  filter(param == "brix") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_all)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "brix", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "brix", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_brix_ndvi_v <- df_test %>% 
  filter(param == "brix_veraison") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "brix_veraison", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "brix_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_brix_pri <- df_test %>% 
  filter(param == "brix") %>%
  ggplot(aes(x = pri, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("PRI", limits = c(-.1,.1))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_all)+
  annotate("text", x = 0, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_pri[df_stats_pri$param == "brix", "rsq"]$rsq, 2),
                         cut(df_stats_pri[df_stats_pri$param == "brix", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_brix_pri_v <- df_test %>% 
  filter(param == "brix_veraison") %>%
  ggplot(aes(x = pri, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("PRI", limits = c(-.1,.1))+
  scale_y_continuous("Observed Brix", limits = c(0,35))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = 0, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_pri[df_stats_pri$param == "brix_veraison", "rsq"]$rsq, 2),
                         cut(df_stats_pri[df_stats_pri$param == "brix_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_plsr <- df_test %>% 
  filter(param == "ta") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous(expression(atop("Predicted", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_color_manual(values = col_all)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "ta", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "ta", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "ta", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_plsr_v <- df_test %>% 
  filter(param == "ta_veraison") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous(expression(atop("Predicted", "Tartaric acid (g/L)")), limits = c(0,15))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,15))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "ta_veraison", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "ta_veraison", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "ta_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_ndvi <- df_test %>% 
  filter(param == "ta") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_color_manual(values = col_all)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "ta", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "ta", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_ndvi_v <- df_test %>% 
  filter(param == "ta_veraison") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "ta_veraison", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "ta_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_pri <- df_test %>% 
  filter(param == "ta") %>%
  ggplot(aes(x = pri, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("PRI", limits = c(-.1,.1))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_color_manual(values = col_all)+
  annotate("text", x = 0, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_pri[df_stats_pri$param == "ta", "rsq"]$rsq, 2),
                         cut(df_stats_pri[df_stats_pri$param == "ta", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ta_pri_v <- df_test %>% 
  filter(param == "ta_veraison") %>%
  ggplot(aes(x = pri, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("PRI", limits = c(-.1,.1))+
  scale_y_continuous(expression(atop("Observed", "Tartaric acid (g/L)")), limits = c(0,20))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = 0, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_pri[df_stats_pri$param == "ta_veraison", "rsq"]$rsq, 2),
                         cut(df_stats_pri[df_stats_pri$param == "ta_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())


plt_sp_ph_plsr <- df_test %>% 
  filter(param == "ph") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous("Predicted pH", limits = c(0,8))+
  scale_y_continuous("Observed pH", limits = c(0,8))+
  scale_color_manual(values = col_all)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "ph", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "ph", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "ph", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ph_plsr_v <- df_test %>% 
  filter(param == "ph_veraison") %>%
  ggplot(aes(x = .pred, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')+
  scale_x_continuous("Predicted pH", limits = c(0,8))+
  scale_y_continuous("Observed pH", limits = c(0,8))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = 0, y = Inf, vjust = 1.35, hjust = 0, 
           label = paste("RMSPE =", round(df_stats[df_stats$param == "ph_veraison", "rmspe"]$rmspe, 2), "%",
                         "\nR² =", round(df_stats[df_stats$param == "ph_veraison", "rsq"]$rsq, 2),
                         cut(df_stats[df_stats$param == "ph_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ph_ndvi <- df_test %>% 
  filter(param == "ph") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous("Observed pH", limits = c(0,8))+
  scale_color_manual(values = col_all)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "ph", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "ph", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())

plt_sp_ph_ndvi_v <- df_test %>% 
  filter(param == "ph_veraison") %>%
  ggplot(aes(x = ndvi, y = truth, color = as.factor(date))) +
  geom_point() +
  geom_smooth(aes(color = NULL), color = 'black', method = 'lm', formula = y~x, se = F) +
  scale_x_continuous("NDVI", limits = c(0.5,1))+
  scale_y_continuous("Observed pH", limits = c(0,8))+
  scale_color_manual(values = col_ver)+
  annotate("text", x = .5, y = Inf, vjust = 2, hjust = 0, 
           label = paste("R² =", round(df_stats_ndvi[df_stats_ndvi$param == "ph_veraison", "rsq"]$rsq, 2),
                         cut(df_stats_ndvi[df_stats_ndvi$param == "ph_veraison", "p"]$p, 
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))))+
  theme_bw()+theme(aspect.ratio = 1, legend.title = element_blank())


#Saving sp figure
png("figures/fig_scatter_berry.png", units="in", height = 10, width = 7, res=500)

ggarrange(as_ggplot(text_grob('NDVI', vjust = 3, hjust = .5, size = 14)),
          as_ggplot(text_grob('PLSR', vjust = 3, hjust = .5, size = 14)),
          plt_sp_brix_ndvi, plt_sp_brix_plsr, plt_sp_ta_ndvi, plt_sp_ta_plsr, plt_sp_ph_ndvi, plt_sp_ph_plsr,
          ncol = 2, nrow = 4, common.legend = T, legend = "right", align = "hv", 
          heights = c(.1,1,1,1,.1,1,1,1,1),
          labels = c("", "", "(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "", "", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))

dev.off()   # Stop writing to the file

#Plot spectra---------------------------------------------------------------

plt_spec_aug230814 <- df_data %>%
  filter(date == as.Date("2023-08-14")) %>%
  pivot_longer(cols = `400`:`2400`, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = interaction(Block, Row, Vine), color = Variety_f))+
  geom_line(show.legend = T)+
  annotate("text",  x=Inf, y = Inf, label = "2023-08-14", vjust = 2, hjust = 1.5)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0, 65))+
  scale_color_viridis_d("Variety", option = "H")+
  theme_bw()

plt_spec_aug230828 <- df_data %>%
  filter(date == as.Date("2023-08-28")) %>%
  pivot_longer(cols = `400`:`2400`, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = interaction(Block, Row, Vine), color = Variety_f))+
  geom_line(show.legend = T)+
  annotate("text",  x=Inf, y = Inf, label = "2023-08-28", vjust = 2, hjust = 1.5)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0, 65))+
  scale_color_viridis_d("Variety", option = "H")+
  theme_bw()

plt_spec_sept230907 <- df_data %>%
  filter(date == as.Date("2023-09-07")) %>%
  pivot_longer(cols = `400`:`2400`, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = interaction(Block, Row, Vine), color = Variety_f))+
  geom_line(show.legend = T)+
  annotate("text",  x=Inf, y = Inf, label = "2023-09-07", vjust = 2, hjust = 1.5)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0, 65))+
  scale_color_viridis_d("Variety", option = "H")+
  theme_bw()

plt_spec_sept230918 <- df_data %>%
  filter(date == as.Date("2023-09-18")) %>%
  pivot_longer(cols = `400`:`2400`, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = interaction(Block, Row, Vine), color = Variety_f))+
  geom_line(show.legend = T)+
  annotate("text",  x=Inf, y = Inf, label = "2023-09-18", vjust = 2, hjust = 1.5)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0, 65))+
  scale_color_viridis_d("Variety", option = "H")+
  theme_bw()

#Saving spec figure
png("figures/fig_spectra_berry.png", units="in", height = 10, width = 12, res=500)

ggarrange(plt_spec_aug230814,
          plt_spec_aug230828,
          plt_spec_sept230907,
          plt_spec_sept230918,
          ncol = 2, nrow = 2, align = "hv", common.legend = T, legend = "right",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"))

dev.off()   # Stop writing to the file

