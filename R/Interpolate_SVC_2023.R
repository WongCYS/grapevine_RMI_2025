require(tidyverse)

#Local directory
wd <- 'C:/Users/p6exk/Box/UCDavis/'

#list of all svc files
file_list <- list.files(paste0(wd, "RMI_grapes/RMI2023/spectroscopy/"), 
                        pattern = ".csv",
                        recursive = F, full.names = T)
#list of svc data
data_list <- list()
for (i in 1:length(file_list)){
  data_list[[i]] <- read.csv(file_list[[i]], check.names = F)
}

#define function to linear interpolate spectra to 1 nm from 400 to 2500 nm
interpolaterfl <- function(df){
  #range to keep spectra
  wl_range <- seq(400, 2500, by = 1)
  
  #get index location for first and last wavebands
  wl_start <- 330
  wl_end <- 2550
  
  #get wavelengths
  wl <- as.numeric(colnames(df))
  
  #get boundaries
  wl_start  <- which.min( abs(wl_start  - wl))
  wl_end <- which.min( abs(wl_end - wl))
  
  #keep only meta data
  df_met <- df[1:wl_start-1]
  
  #convert df format wide to long
  df <- df %>%
    pivot_longer(cols = wl_start:wl_end, names_to = "wl", values_to = "rfl")
  
  #interpolate data to 1 nm
  df_interpolate <- approx(df$wl, df$rfl, xout=wl_range, method = "linear")
  
  #extract interpolated data
  df_interp <- data.frame(wl = df_interpolate$x, rfl = df_interpolate$y)
  
  #convert df format long to wide
  df_interp <- df_interp %>%
    pivot_wider(names_from = wl, values_from = rfl)
  
  #merge with meta data
  df_interp <- merge(df_met, df_interp)
  
  return(df_interp)
}


#Loop the function to interpolate data
df_svc <- data.frame()
for (i in 1:length(data_list)){
  df <- data_list[[i]]
  for (j in 1:nrow(df)){
    print(paste("list", i, "out of", length(data_list), "row", j, "out of", nrow(df)))
    temp_data <- interpolaterfl(df[j,])
    df_svc <- rbind(df_svc, temp_data)
  }
}



#quick spectra plot
df_svc %>%
  pivot_longer(cols = 9:2109, names_to = "wl", values_to = "rfl") %>%
  ggplot(aes(x = as.numeric(wl), y = rfl, group = ID, color = Variety))+
  geom_line(show.legend = T)+
  scale_x_continuous("Wavelength (nm)")+
  scale_y_continuous("Reflectance (%)", limits = c(0,100))+
  facet_wrap(~date)+
  theme_bw()

#save df
write.csv(df_svc, paste0(wd, "RMI_grapes/RMI2023/spectroscopy/svc_clean_interp.csv"), row.names = F)
