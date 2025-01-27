# Author: Murray SA Thompson
# Contact: murray.thompson@cefas.gov.uk
# Version: 1
# January 2025
# prepared for EMODnet


# rm(list=ls()); gc()

pkgs = c("tidyverse", "data.table", "vroom", "mapplots", "patchwork", "ggpubr") #
for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
} 


### Here is how to load and pre-process all of the fish survey data downloaded from https://data.cefas.co.uk/view/21421
# 
# setwd('Y:/C8424_Pyramids_Of_Life_SMMR/Data_Storage/Biological/fish survey data/Lynam dataproduct')
# load('survey_data_taxonomy.Rdata') #txdf
# 
# # read all survey data txt files:
# list_of_files = list.files(
#   pattern = "\\.txt$", 
#   full.names = TRUE)
# 
# # Read all the files and create a FileName column to store filenames
# DT = rbindlist(sapply(list_of_files, fread, simplify = FALSE),
#                use.names = TRUE, fill=TRUE, idcol = "FileName")
# 
# alldf = DT %>% as_tibble() %>%
#   left_join(txdf, c('SciName'='taxa')) %>%
#   rename(taxa=new_taxa, year=Year, latitude=ShootLat_degdec, longitude=ShootLong_degdec) %>%
#   mutate(ices=ices.rect2(longitude, latitude),
#          # to iron out missing biomass info
#          ind_weight_g = exp(log(LWRa)+LWRb*log(FishLength_cm)),# \1000 to convert to kg
#          DensBiom_kg_Sqkm = case_when(is.na(DensBiom_kg_Sqkm) ~ (ind_weight_g/1000)*DensAbund_N_Sqkm,
#                                       TRUE ~ DensBiom_kg_Sqkm),
#          # add body mass bins
#          bin_num = case_when(ind_weight_g < 6.73e-8 ~ 1,
#                              ind_weight_g < 4.38e-7 ~ 2,
#                              ind_weight_g < 2.85e-6 ~ 3,
#                              ind_weight_g < 1.85e-5 ~ 4,
#                              ind_weight_g < 1.20e-4 ~ 5,
#                              ind_weight_g < 7.83e-4 ~ 6,
#                              ind_weight_g < 5.09e-3 ~ 7,
#                              ind_weight_g < 3.31e-2 ~ 8,
#                              ind_weight_g < 2.16e-1 ~ 9,
#                              ind_weight_g < 1.40e+0 ~ 10,
#                              ind_weight_g < 9.12e+0 ~ 11,
#                              ind_weight_g < 5.93e+1 ~ 12,
#                              ind_weight_g < 3.86e+2 ~ 13,
#                              ind_weight_g < 2.51e+3 ~ 14,
#                              ind_weight_g < 1.63e+4 ~ 15,
#                              ind_weight_g < 1.06e+5 ~ 16,
#                              ind_weight_g < 6.90e+5 ~ 17,
#                              ind_weight_g < 4.49e+6 ~ 18,
#                              ind_weight_g < 2.92e+7 ~ 19,
#                              ind_weight_g < 1.90e+8 ~ 20))

# save(alldf, file='Y:/C8424_Pyramids_Of_Life_SMMR/Data_Storage/Biological/EMODnet/combined_survey_data_for_EMODnet.R')

#################################################
#setwd('Y:/C8424_Pyramids_Of_Life_SMMR/Data_Storage/Biological/EMODnet')
setwd("C:\Users\KC05\CEFAS\Murray Thompson (Cefas) - data")

# load processed survey data
load('combined_survey_data_for_EMODnet.R')

# NMDS scores generated for this study: https://doi.org/10.5194/essd-2024-102
# Note to self - Feeding guilds_table_08_23.csv needs to be published here: https://data.cefas.co.uk/view/21771
load('guilds_cluster_NMDS_08_23.Rdata')

glds=sub_glds
glds$bin_num = glds$bin_number

pred_scrs = as.data.frame(nmds.sor$points)
pred_scrs$pred_id = rownames(pred_scrs)
rownames(pred_scrs) = NULL

tbl = sub_glds %>% 
  left_join(pred_scrs, 'pred_id') %>%
  #left_join(guild_hrachy, c('pred_taxa', 'bin_number')) %>%
  dplyr::select(pred_taxa, bin_number, MDS1, MDS2) %>%
  rename(taxa=pred_taxa) %>%
  arrange(taxa)

# Gobiidae and Ammodytidae taxa resolved to species across different data:
go_am_tx = c("Pomatoschistus minutus",
             "Hyperoplus lanceolatus", 
             "Gobius niger", 
             "Gobiusculus flavescens", 
             "Gobius paganellus",
             "Ammodytes dubius",
             "Aphia minuta")

# select data to use - Otter trawls between 1997-2020
# GNSIntOT3 is subset because it overlaps spatially with GNSIntOT1, the latter is more comparable in timing with other surveys
otterdf = alldf %>%
  filter(Gear == 'GOV',
         year %in% c(1997:2020),
         !Survey_Acronym == 'GNSIntOT3',
         DensBiom_kg_Sqkm >0,
         !is.na(bin_num)) %>%
  distinct() %>%
  # resolve taxa, ready for integration with feeding guild information
  mutate(taxa = case_when(species %in% go_am_tx ~ species,
                          family == 'Gobiidae' ~ 'Gobiidae',
                          family == 'Ammodytidae' ~ 'Ammodytes',
                          TRUE ~ species)) %>% 
  left_join(tbl, c('bin_num'='bin_number', 'taxa')) 

# calculate biomass weighted mean axis scores
surv_wgts = otterdf %>%
  group_by(HaulID, Survey_Acronym, ices, year) %>%
  summarise(b_mds1 = weighted.mean(MDS1, DensBiom_kg_Sqkm, na.rm=T),
            b_mds2 = weighted.mean(MDS2, DensBiom_kg_Sqkm, na.rm=T)) %>%
  pivot_longer(cols = -(HaulID:year),
               names_to = 'scrs',
               values_to = 'value') %>%
  separate(scrs, c('response', 'axes'), sep='_', remove=T) %>%
  pivot_wider(names_from = response,
              values_from = value)

##############################################################

# ices-level means for converting to NetCDF 
dat_for_netcdf = surv_wgts %>%
  group_by(ices, year, axes) %>%
  summarise(av_bw_axis_scores=mean(b))


##############################################################
### spatial and temporal plots

# estimating temporal correlations

tempdf_axes = data.frame(ices=NA, axes=NA,
                    b_cor=NA, n_cor=NA,
                    b_p=NA, n_p=NA)

for(ss in unique(dat_for_netcdf$ices)){ #ss=dat_for_netcdf$ices[1]

  dat = subset(dat_for_netcdf, ices %in% ss)

  for(axs in unique(dat_for_netcdf$axes)){ # axs='mds1'

    datg = subset(dat, axes %in% axs)

    b_cor = try(cor.test(datg$year, datg$av_bw_axis_scores, use='complete.obs'), silent=T)
    n_cor = try(cor.test(datg$year, datg$av_n, use='complete.obs'), silent=T)

    ss_df = data.frame(ices=ss, axes=axs,
                        b_cor=ifelse(is(b_cor,"try-error"), NA, b_cor$estimate),
                        n_cor=ifelse(is(n_cor,"try-error"), NA, n_cor$estimate),

                        b_p=ifelse(is(b_cor,"try-error"), NA, b_cor$p.value),
                        n_p=ifelse(is(n_cor,"try-error"), NA, n_cor$p.value))

    tempdf_axes = rbind(tempdf_axes, ss_df)
  }
}
tempdf_axes = tempdf_axes[-1,]


##############################################################
#### map plots

# extract ices.rect coordinates
coords = ices.rect(tempdf_axes$ices)
tempdf_axes = cbind(ungroup(tempdf_axes), coords)

coords = ices.rect(dat_for_netcdf$ices)
dat_for_netcdf = cbind(ungroup(dat_for_netcdf), coords)

# set map limits with full dataset
xlims = range(dat_for_netcdf$lon)
ylims = range(dat_for_netcdf$lat)


##########
# spatial 

# subset No guild to plot
toplot = unique(dat_for_netcdf$axes)

world_shp = sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

##########
# spatial means

# consider subsetting years to check each 
plt_df = dat_for_netcdf #%>%
  #filter(year==2010)

plot_list_b = list()
for(axs in unique(plt_df$axes)){ #axs = unique(plt_df$axes)[1]

  pltdat = subset(plt_df, axes==axs & !is.na(av_bw_axis_scores)) %>%
    group_by(ices, lon, lat) %>%
    summarise(av_bw_axis_scores=mean(av_bw_axis_scores))

  lims = range(0, max(pltdat$av_bw_axis_scores))

  p = pltdat %>% ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = av_bw_axis_scores)) +
    viridis::scale_fill_viridis(name = '',
      option = 'rocket') +
    geom_sf(data = world_shp,
            fill = 'black',
            color = 'black',
            size = 0.1) +
    coord_sf(xlim = xlims, ylim = ylims)+
    labs(x='', y='', title=axs) +
    guides(fill = guide_colourbar(barwidth = 1)) +
    theme(panel.background = element_rect(fill = 'grey80'),
          panel.border = element_rect(colour='black', fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot_list_b[[axs]] = p
}

av_plts_b = do.call(ggarrange, c(plot_list_b[toplot], align='hv', ncol=2, nrow=1)) %>%
  annotate_figure(fig.lab = 'Spatial biomass', fig.lab.size = 14)

##########
# temporal
plot_list_b_temp = list()
for(axs in unique(tempdf_axes$axes)){ #axs ='Benthivore'

  pltdat = subset(tempdf_axes, axes==axs & !is.na(b_cor))
  lims = range(pltdat$b_cor)

  p = pltdat %>% ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = b_cor)) +
    scale_fill_gradient2("", low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = lims) +
    geom_tile(aes(x = lon, y = lat), fill=NA, size=0.2, color = 'black', data = subset(pltdat, b_p<0.05)) +
    geom_sf(data = world_shp,
            fill = 'black',
            color = 'black',
            size = 0.1) +
    coord_sf(xlim = xlims, ylim = ylims)+
    labs(x='', y='', title=axs) +
    guides(fill = guide_colourbar(barwidth = 1)) +
    theme(panel.background = element_rect(fill = 'grey80'),
          panel.border = element_rect(colour='black', fill=NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot_list_b_temp[[axs]] = p
}

temp_plts_b = do.call(ggarrange, c(plot_list_b_temp[toplot], align='hv', ncol=2, nrow=1)) %>%
  annotate_figure(fig.lab = 'Temporal biomass', fig.lab.size = 14)

av_plts_b + temp_plts_b +
    plot_layout(nrow=2)
