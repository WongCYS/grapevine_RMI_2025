library(tidyverse)
library(viridis)

#Local directory
wd <- 'C:/Users/cyswong/Box/UCDavis/'

#read in data
df_svc <- read.csv(paste0(wd, "RMI_grapes/RMI2021/spectroscopy/processed_data/svc_clean_interp.csv"), check.names = F)

#read in origin info
df_origin <- read.csv(paste0(wd, "RMI_grapes/RMI2021/VarietyOrigins.csv"))

#merge spectra with origins
df_data <- merge(df_svc, df_origin, by = "Variety")

#quick spectra plot
df_svc %>%
  pivot_longer(cols = 9:2109, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = ID, color = Variety))+
  geom_line(show.legend = T)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)")+
  facet_wrap(~date)+
  theme_bw()


#PCA via tidymodels
library(tidymodels)


#Keep only certain columns
df_pca <- df_data[,-c(2:7)]

#screen for na
df_pca <- df_pca[complete.cases(df_pca), ] 

df_pca_date <- df_pca %>%
  #filter(date == as.Date("2021-06-24"))
  #filter(date == as.Date("2021-07-15"))
  #filter(date == as.Date("2021-08-12"))
  filter(date == as.Date("2021-09-09"))

#recipe for pca
pca_rec <- recipe(~., data = df_pca) %>%
  update_role(Variety, date, Region_of_Origin, new_role = 'id') %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors())

#prep pca
pca_prep <- prep(pca_rec)

#pca and eigenvalues
tidied_pca <- tidy(pca_prep, 2)
eig_pca <- tidy(pca_prep, 2, type = 'variance')

#plot PC breakdown
tidied_pca %>%
  filter(component %in% paste0("PC", 1:6)) %>%
  mutate(component = fct_inorder(component),
         terms = parse_number(terms)) %>%
  ggplot(aes(as.numeric(terms), value)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~component) +
  labs(y = NULL)

#plot biplot
bake(pca_prep, new_data = NULL) %>%
  ggplot(aes(PC1, PC2, label = Variety)) +
  geom_point(aes(color = Region_of_Origin, shape = date), alpha = 0.7, size = 2, show.legend = F) +
  geom_text(check_overlap = TRUE, hjust = "inward") +
  labs(color = NULL)+
  scale_color_viridis(discrete = T)




#plot variance explained
eig_pca %>%
  filter(terms == 'percent variance') %>%
  ggplot(aes(x = component, y = value))+
  geom_line()+
  scale_y_continuous('Percent Variance')+
  scale_x_continuous(limits = c(0,5))
