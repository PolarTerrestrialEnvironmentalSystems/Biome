
###################################
# R version 4.4.1                 #
# Operating System: Windows 10    #
# Code for data Visualization     #
# Supplement to: Li, C., Dallmeyer, A., Ni, J., Chevalier, M., Willeit, M., Andreev, A. A., Cao, X., Schild, L., Heim, B., and Herzschuh, U.: Global biome changes over the last 21,000 years inferred from model-data comparisons, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-1862, 2024. #
# Contact: Chenzhi Li (chenzhi.li@awi.de)[Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Potsdam, Germany 2024] #
###################################

# Figure List:
## 1. Figure 1: Spatial patterns of megabiome distributions at 0 cal. ka BP and their agreement with modern potential natural megabiomes
## 2. Figure 2. Spatial distributions of megabiomes at 21, 16, 13, 9, 6, and 3 cal. ka BP
## 3. Figure 3：Temporal changes in the latitudinal location of megabiome and ice-sheet
## 4. Figure 4：Spatiotemporal patterns of EMD between the pollen-based reconstructions and ESM-based simulation ensemble
## 5. Appendix Figure A1：Spatial distribution and sources of fossil pollen records in the LegacyPollen 2.0 dataset
## 6. Appendix Figure A2: Spatial patterns of megabiome distributions at 0 cal. ka BP and their agreement with modern potential natural megabiomes
## 7. Appendix Figure A3：Differences in bioclimatic variables between ESM-based simulations at 0ka and observations

# Note: Biomes and their abbreviations and codes 
## 1 - Tropical forest (TRFO); 2 - Subtropical forest (WTFO); 3 - Temperate forest (TEFO); 4 - Boreal forest (BOFO); 
## 5 - (Warm) savanna and dry woodland (SAVA); 6 - Grassland and dry shrubland (STEP); 7 - (Warm) desert (DESE); 8 - Tundra and polar desert (TUND)

# Note: Bioclimatic variables
## BIO10 = Mean Temperature of Warmest Quarter; BIO11 = Mean Temperature of Coldest Quarter; BIO18 = Precipitation of Warmest Quarter; BIO19 = Precipitation of Coldest Quarter

# choose directory:
#setwd("~/Supplementary code/Visualization")  # select the folder "Supplementary code" from your directory

setwd("N:/bioing/user/cli/Biomization_Globally/data-model verification/Figure in the manuscript/Second submission")

# install packages if not installed
#install.packages(c("ggplot2", "rnaturalearth", "rnaturalearthdata", "dplyr", "ggpubr", "ggnewscale", "tidyr"))

# loading packages
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggpubr)
library(ggnewscale)
library(tidyr)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 1. Figure 1: Spatial patterns of megabiome distributions at 0 cal. ka BP and their agreement with modern potential natural megabiomes  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_pollen_ESM_site_0ka                   <- read.csv2("data/Fig. 1 data.1-Reconstructed and simulated megabiome at 0ka for record.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_megabiomes_site    <- read.csv2("data/Fig. 1 data.2-Modern potential natural megabiomes for record.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_megabiomes_grid    <- read.csv2("data/Fig. 1 data.3-Modern potential natural biomes for grid (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_pollen_ESM_site_0ka                   <- type.convert(Biome_pollen_ESM_site_0ka, as.is = TRUE) 
Modern_potential_natural_megabiomes_site    <- type.convert(Modern_potential_natural_megabiomes_site, as.is = TRUE) 
Modern_potential_natural_megabiomes_grid    <- type.convert(Modern_potential_natural_megabiomes_grid, as.is = TRUE) 

# Join datasets
Modern_potential_natural_pollen_ESM_megabiomes_site <- left_join(Biome_pollen_ESM_site_0ka, Modern_potential_natural_megabiomes_site[ ,c(1,7)], by = "Dataset_ID")
names(Modern_potential_natural_pollen_ESM_megabiomes_site)[c(7,14:15)] <- c("Reconstruction", "Simulation_ensemble", "PNV") ## Rename columns for clarity

# Calculate the agreement
{
  # Initialize an empty result matrix
  Modern_potential_natural_pollen_ESM_megabiomes_site_agreement <- data.frame(matrix(NA, nrow=nrow(Modern_potential_natural_pollen_ESM_megabiomes_site), ncol= ncol(Modern_potential_natural_pollen_ESM_megabiomes_site)-1, byrow=TRUE))
  colnames(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement) <- names(Modern_potential_natural_pollen_ESM_megabiomes_site)[1:14]
  
  # Copy metadata from the original dataset to the agreement matrix.
  Modern_potential_natural_pollen_ESM_megabiomes_site_agreement[ ,1:6] <- Modern_potential_natural_pollen_ESM_megabiomes_site[ ,1:6]
  
  # Loop through columns 7 to 14 to check agreement
  for (i in 7:14) {
    
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement[ ,i] <- ifelse(Modern_potential_natural_pollen_ESM_megabiomes_site[ ,i] == Modern_potential_natural_pollen_ESM_megabiomes_site[ ,15], "Yes", "No")
    
  }
  
}

# Map settings and customization
{
  # Convert the Biome column to a factor with specified levels
  Modern_potential_natural_megabiomes_grid$Biome <- factor(Modern_potential_natural_megabiomes_grid$Biome,
                                                           levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
  
  # Define the colors corresponding to each megabiome category
  cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
            "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
  
  # Define the breaks (categories) for megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
  
  # Define the labels for the megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
          "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
  
  # Load coastline data
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  
}


# 1.1 Spatial patterns of megabiome distributions at 0 cal. ka BP
# -------------------------------------------------------------------------------------------------

{
  # (1) Pollen-based reconstruction
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$Reconstruction <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$Reconstruction, 
                                                                                 levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_Reconstruction_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= '(a) Megabiome distribution of pollen-based reconstruction', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= Reconstruction), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_Reconstruction_site_map
      
    }
    
  }
  
  
  # (2) ESM-based simulation ensemble
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$Simulation_ensemble <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$Simulation_ensemble, 
                                                                                      levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_Simulation_ensemble_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of ESM-based simulation ensemble', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= Simulation_ensemble), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_Simulation_ensemble_site_map
      
    }
    
  }
  
  
}


# 1.2 Megabiome agreement with modern potential natural biomes
# -------------------------------------------------------------------------------------------------

{
  # (1) Pollen-based reconstruction
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$Reconstruction <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$Reconstruction, 
                                                                                           levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_Reconstruction_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= '(b) Megabiome agreement of pollen-based reconstruction', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= Reconstruction), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_Reconstruction_site_map
      
    }
    
  }
  
  
  # (2) ESM-based simulation ensemble
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$Simulation_ensemble <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$Simulation_ensemble, 
                                                                                                levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_Simulation_ensemble_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of ESM-based simulation ensemble', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= Simulation_ensemble), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_Simulation_ensemble_site_map
      
    }
    
  }
  
  
}


# -------------------------------------------------------------------------------------------------
#Arranging plots
Biome_0ka_agreement_site_pollen_models_map   <- ggarrange(ggarrange(Biome_0ka_Reconstruction_site_map, Biome_0ka_Simulation_ensemble_site_map, 
                                                                    nrow = 1, ncol = 2, common.legend = T, legend = "bottom"),
                                                          ggarrange(Biome_agreement_0ka_Reconstruction_site_map, Biome_agreement_0ka_Simulation_ensemble_site_map,
                                                                    nrow = 1, ncol = 2, common.legend = T, legend = "bottom"),
                                                          nrow = 2, ncol = 1, heights = c(0.515, 0.485), common.legend = F)

# Save the final figure as a PNG file
ggsave("result/Figure 1-Spatial patterns of megabiome distributions at 0 cal. ka BP and their agreement with modern potential natural megabiomes.png", width = 14, height = 8, units = "in", dpi = 300)



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 2. Figure 2. Spatial distributions of megabiomes at 21, 16, 13, 9, 6, and 3 cal. ka BP  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_timeslice_pollen_ensemble_simulation_grid    <- read.csv2("data/Fig. 2 data.1-Reconstructed and simulated megabiome at selected timeslices for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid        <- read.csv2("data/Fig. 2 data.2-Ice sheet of ICE5G ICE6G and GLAC1D since 21 ka BP for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_timeslice_pollen_ensemble_simulation_grid    <- type.convert(Biome_timeslice_pollen_ensemble_simulation_grid, as.is = TRUE) 
ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid        <- type.convert(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, as.is = TRUE) 

# Map settings and customization
{
  # Define the colors corresponding to each megabiome category
  cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
            "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
  
  # Define the breaks (categories) for megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
  
  # Define the labels for the megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
          "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
  
  # Load coastline data
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  
}


# 2.1 Pollen-based reconstruction
# -------------------------------------------------------------------------------------------------
{
  # Convert biome to a factor
  Biome_timeslice_pollen_ensemble_simulation_grid$Biome_pollen  <- factor(Biome_timeslice_pollen_ensemble_simulation_grid$Biome_pollen,
                                                                          levels = c("1", "2", "3", "4", "5", "6", "7", "8"))  
  
  # (1) 21 cal. ka BP
  {
    Biome_timeslice_pollen_grid_21ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "(a) Pollen-based reconstruction at 21 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 21), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==21), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_21ka_map
  }
  
  # (2) 16 cal. ka BP
  {
    Biome_timeslice_pollen_grid_16ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "Pollen-based reconstruction at 16 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 16), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==16), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_16ka_map
  }
  
  
  # (3) 13 cal. ka BP
  {
    Biome_timeslice_pollen_grid_13ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "Pollen-based reconstruction at 13 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 13), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==13), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_13ka_map
  }
  
  # (4) 9 cal. ka BP
  {
    Biome_timeslice_pollen_grid_9ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "Pollen-based reconstruction at 9 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 9), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==9), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_9ka_map
  }
  
  
  # (5) 6 cal. ka BP
  {
    Biome_timeslice_pollen_grid_6ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "Pollen-based reconstruction at 6 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 6), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==6), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_6ka_map
  }
  
  # (6) 3 cal. ka BP
  {
    Biome_timeslice_pollen_grid_3ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "Pollen-based reconstruction at 3 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 3), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==3), aes(x=Longitude, y=Latitude, fill= Biome_pollen), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_pollen_grid_3ka_map
  }
  
}


# 2.2 ESM-based simulation ensemble
# -------------------------------------------------------------------------------------------------
{
  # Convert biome to a factor
  Biome_timeslice_pollen_ensemble_simulation_grid$Biome_ensemble_simulation  <- factor(Biome_timeslice_pollen_ensemble_simulation_grid$Biome_ensemble_simulation,
                                                                                       levels = c("1", "2", "3", "4", "5", "6", "7", "8"))  
  
  # (1) 21 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_21ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "(b) ESM-based simulation ensemble at 21 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 21), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==21), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_21ka_map
  }
  
  # (2) 16 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_16ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "ESM-based simulation ensemble at 16 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 16), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==16), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_16ka_map
  }
  
  
  # (3) 13 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_13ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "ESM-based simulation ensemble at 13 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 13), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==13), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_13ka_map
  }
  
  # (4) 9 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_9ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "ESM-based simulation ensemble at 9 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 9), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==9), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_9ka_map
  }
  
  
  # (5) 6 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_6ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "ESM-based simulation ensemble at 6 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 6), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==6), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))+ 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_6ka_map
  }
  
  # (6) 3 cal. ka BP
  {
    Biome_timeslice_ensemble_simulation_grid_3ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.3) +
      labs(title= "ESM-based simulation ensemble at 3 cal. ka BP", x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_tile(data= subset(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, Timeslice == 3), aes(x=Longitude, y=Latitude, fill = Ice_ensemble), alpha = 0.8) +
      scale_fill_manual(name="Ice-sheet", values= c("Ice-sheet"  = "lightblue"), labels = "") +
      guides(fill = guide_legend(override.aes = list(), title.position = "bottom", title.theme = element_text(size=14, hjust = 0.5, vjust = 0.5), order = 2)) +
      new_scale_fill() +
      geom_tile(data= subset(Biome_timeslice_pollen_ensemble_simulation_grid, Timeslice ==3), aes(x=Longitude, y=Latitude, fill= Biome_ensemble_simulation), alpha = 0.8) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm")) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) + 
      guides(fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE, order = 1)) 
    
    Biome_timeslice_ensemble_simulation_grid_3ka_map
  }
  
}


# -------------------------------------------------------------------------------------------------
# Arranging plots
Biome_timeslice_pollen_ensemble_simulation_grid_map   <- ggarrange(Biome_timeslice_pollen_grid_21ka_map, Biome_timeslice_ensemble_simulation_grid_21ka_map,
                                                                   Biome_timeslice_pollen_grid_16ka_map, Biome_timeslice_ensemble_simulation_grid_16ka_map,
                                                                   Biome_timeslice_pollen_grid_13ka_map, Biome_timeslice_ensemble_simulation_grid_13ka_map,
                                                                   Biome_timeslice_pollen_grid_9ka_map,  Biome_timeslice_ensemble_simulation_grid_9ka_map,
                                                                   Biome_timeslice_pollen_grid_6ka_map,  Biome_timeslice_ensemble_simulation_grid_6ka_map, 
                                                                   Biome_timeslice_pollen_grid_3ka_map,  Biome_timeslice_ensemble_simulation_grid_3ka_map,
                                                                   nrow = 6, ncol = 2, common.legend = T, legend = "bottom")

# Save the final figure as a PNG file
ggsave("result/Figure 2-Spatial distributions of megabiomes at 21, 16, 13, 9, 6, and 3 cal. ka BP.png", width = 13, height = 21, units = "in", dpi = 300)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 3. Figure 3：Temporal changes in the latitudinal location of megabiome and ice-sheet ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_timeslice_pollen_grid                 <- read.csv2("data/Fig. 3 data.1-Pollen-based reconstructed biome at each timeslice for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Biome_timeslice_ensemble_simulation_grid    <- read.csv2("data/Fig. 3 data.2-Ensemble simulated biome at each timeslice for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid <- read.csv2("data/Fig. 3 data.3-Ice sheet of ICE5G ICE6G and GLAC1D since 21 ka BP for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_PNV_grid                             <- read.csv2("data/Fig. 3 data.4-Modern potential megabiomes for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion 
Biome_timeslice_pollen_grid                 <- type.convert(Biome_timeslice_pollen_grid, as.is = TRUE) 
Biome_timeslice_ensemble_simulation_grid    <- type.convert(Biome_timeslice_ensemble_simulation_grid, as.is = TRUE)
ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid <- type.convert(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, as.is = TRUE) 
Modern_PNV_grid                             <- type.convert(Modern_PNV_grid, as.is = TRUE) 

# Filter datasets
## Join pollen and ensemble simulation data by coordinates and timeslice, keeping relevant columns
Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka   <- left_join(Biome_timeslice_pollen_grid[ ,c(1:3,6)], Biome_timeslice_ensemble_simulation_grid[Biome_timeslice_ensemble_simulation_grid$Source == "Balanced ensemble", c(2:4,7)], by= c("Longitude", "Latitude", "Timeslice"))
names(Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka)[4:5] <- c("Biome_pollen", "Biome_ensemble_simulation") ## Rename columns for clarity
Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka   <- Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka[complete.cases(Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka), ] ## Remove rows with missing data

ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid  <- na.omit(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[ , c(1:3,7)]) ## Filter ice sheet dataset to keep relevant columns and assign "Ice-sheet" as a fixed value
names(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid)[4]  <- "Value"
ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid$Value <- "Ice-sheet"

# Combine datasets
Biome_timeslice_pollen_grid_0_21ka <- Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka[ ,1:4] ## Create a unified dataset for pollen and ice sheet data
names(Biome_timeslice_pollen_grid_0_21ka)[4] <- "Value"
Biome_Ice_timeslice_pollen_grid_0_21ka <- rbind(Biome_timeslice_pollen_grid_0_21ka, ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid)

# Convert biome and timeslice values into factors with specific levels and labels
Biome_Ice_timeslice_pollen_grid_0_21ka$Value <- factor(Biome_Ice_timeslice_pollen_grid_0_21ka$Value,
                                                       levels = c("Ice-sheet", "1", "2", "3", "4", "5", "6", "7", "8"),
                                                       labels = c("Ice-sheet", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))                    
Biome_Ice_timeslice_pollen_grid_0_21ka$Timeslice <- factor(Biome_Ice_timeslice_pollen_grid_0_21ka$Timeslice ,
                                                           levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", 
                                                                      "12", "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5", "21"))                    

# Repeat similar steps for ensemble simulation data
Biome_timeslice_ensemble_simulation_grid_0_21ka <- Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka[ ,c(1:3,5)]
names(Biome_timeslice_ensemble_simulation_grid_0_21ka)[4] <- "Value"
Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka <- rbind(Biome_timeslice_ensemble_simulation_grid_0_21ka, ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid)
Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka$Value <- factor(Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka$Value,
                                                                    levels = c("Ice-sheet", "1", "2", "3", "4", "5", "6", "7", "8"),
                                                                    labels = c("Ice-sheet", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))                    
Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka$Timeslice <- factor(Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka$Timeslice ,
                                                                        levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", 
                                                                                   "12", "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5", "21"))                    

# Statistical modern latitudinal position
PNV_median_latitude_site_grid <- na.omit(left_join(Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka[Biome_timeslice_pollen_ensemble_simulation_grid_0_21ka$Timeslice == 0,1:2], Modern_PNV_grid, by = c("Longitude", "Latitude"))) %>%
  group_by(PNV_MPI_ESM_V3_GLAC1D) %>%
  summarize(Latitude_median = median(abs(Latitude))) ## Calculate the median latitude for each modern potential vegetation (PNV) type, based on sites at timeslice = 0
names(PNV_median_latitude_site_grid)[1] <- "Value" ## Rename column for clarity
PNV_median_latitude_site_grid$Value <- factor(PNV_median_latitude_site_grid$Value,
                                              levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                                              labels = c("TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))   ## Convert PNV values into descriptive factors

# -------------------------------------------------------------------------------------------------

# plot
{
  # set the y limiit manually
  {
    # Calculate Q1 (25th percentile) and Q3 (75th percentile) for pollen data grouped by Value (biome type) and Timeslice
    biome_pollen_latitude_Q1_Q3_grid <- Biome_Ice_timeslice_pollen_grid_0_21ka %>%
      group_by(Value, Timeslice) %>%
      summarize(Q1 = quantile(abs(Latitude), c(0.25)),
                Q3 = quantile(abs(Latitude), c(0.75))) %>%
      as.data.frame()
    
    # Calculate Q1 and Q3 for ensemble simulation data grouped by Value and Timeslice
    biome_ensemble_simulation_latitude_Q1_Q3_grid <- Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka %>%
      group_by(Value, Timeslice) %>%
      summarize(Q1 = quantile(abs(Latitude), c(0.25)),
                Q3 = quantile(abs(Latitude), c(0.75))) %>%
      as.data.frame()
    
    # Combine pollen and ensemble simulation Q1/Q3 data into a single dataset
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_grid   <- rbind(biome_pollen_latitude_Q1_Q3_grid, biome_ensemble_simulation_latitude_Q1_Q3_grid)
    
    # Calculate the minimum Q1 and maximum Q3 across all Timeslices for each Value (biome type)
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid <- biome_pollen_ensemble_simulation_latitude_Q1_Q3_grid %>%
      group_by(Value) %>%
      summarize(Q1_min = min(Q1),
                Q3_max  = max(Q3)) %>%
      as.data.frame()
    
    # Round Q1_min and Q3_max to the nearest multiple of 5 for ymin and ymax limits
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid$Q1_range <- as.integer(biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid$Q1_min / 5) *5
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid$Q3_range <- (as.integer(biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid$Q3_max/ 5) + 1) *5
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid$n        <- 5
    
    # Rename columns for compatibility with ggplot2 scaling
    names(biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid)[c(1,4,5)] <- c("Panel", "ymin", "ymax")
    
    # Copy the dataset for manual adjustments
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final <- biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid
    
    # Manually set ymin and ymax values for each biome type
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final$ymin <- c(55, 0, 5, 25, 35, 0, 20, 15, 40)
    biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final$ymax <- c(90, 35, 40, 60, 70, 35, 55, 50, 75)
    
    # Create a list of data frames for each biome type, including ymin, ymax, and n
    biome_pollen_ensemble_simulation_latitude_scales_grid_list <- list('Ice-sheet' = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[1, c(1,4:6)],
                                                                       TRFO = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[2, c(1,4:6)],
                                                                       WTFO = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[3, c(1,4:6)],
                                                                       TEFO = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[4, c(1,4:6)],
                                                                       BOFO = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[5, c(1,4:6)],
                                                                       SAVA = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[6, c(1,4:6)],
                                                                       STEP = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[7, c(1,4:6)],
                                                                       DESE = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[8, c(1,4:6)],
                                                                       TUND = biome_pollen_ensemble_simulation_latitude_Q1_Q3_range_grid_final[9, c(1,4:6)])
    
    # Create ggplot2 scale_y_continuous objects for each biome type with specific ymin, ymax, and breaks
    biome_pollen_ensemble_simulation_latitude_scales_grid_scales <- lapply(biome_pollen_ensemble_simulation_latitude_scales_grid_list, function(x) {
      scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n, expand = c(0,0))
    })
    
  }
  
  
  {
    # Pollen-based reconstruction 
    Latitude_pollen_ice_grid_boxplot <- ggplot(data=Biome_Ice_timeslice_pollen_grid_0_21ka, aes(x = Timeslice, y = abs(Latitude))) +
      geom_boxplot(aes(fill = Value), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", alpha = 0.5) +
      stat_summary(fun=median, geom="line", aes(group=Value), position=position_dodge(width=0.9), color="black", size=1) +
      geom_hline(data = PNV_median_latitude_site_grid, aes(yintercept = Latitude_median), color = "red", linetype = "dashed") +
      labs(x="Age (cal. ka BP)", y = "", title = "(a) Pollen-based reconstruction") +
      scale_x_discrete(breaks = seq(0,21,by=1),
                       labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_fill_manual(name="", values=c("Ice-sheet" = "lightblue", "TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#ffeda0", "DESE"="#FFFF00", "TUND"="#6a3d9a"), guide = "none") +
      facet_grid(Value ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
                 
      ) +
      ggh4x::facetted_pos_scales(y = biome_pollen_ensemble_simulation_latitude_scales_grid_scales) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(10, 0, 0, 0), "pt"), size=18),
            axis.text.x=element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 10, 0, 0), "pt")),
            axis.text.y=element_text(size=16)) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(1.2, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=16, face="bold"),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 15, 0, 0), "pt")),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) 
  }
  
  
  {
    # Ensemble simulations 
    
    Latitude_ensemble_simulation_ice_grid_boxplot <- ggplot(data=Biome_Ice_timeslice_ensemble_simulation_grid_0_21ka, aes(x = Timeslice, y = abs(Latitude))) +
      geom_boxplot(aes(fill = Value), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", alpha = 0.5) +
      stat_summary(fun=median, geom="line", aes(group=Value), position=position_dodge(width=0.9), color="black", size=1) +
      geom_hline(data = PNV_median_latitude_site_grid, aes(yintercept = Latitude_median), color = "red", linetype = "dashed") +
      labs(x="Age (cal. ka BP)", y = "", title = "(b) ESM-based simulation ensemble") +
      scale_x_discrete(breaks = seq(0,21,by=1),
                       labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_fill_manual(name="", values=c("Ice-sheet" = "lightblue", "TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#ffeda0", "DESE"="#FFFF00", "TUND"="#6a3d9a"), guide = "none") +
      facet_grid(Value ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
                 
      ) +
      ggh4x::facetted_pos_scales(y = biome_pollen_ensemble_simulation_latitude_scales_grid_scales) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(10, 0, 0, 0), "pt"), size=18),
            axis.text.x=element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 10, 0, 0), "pt")),
            axis.text.y=element_text(size=16)) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(1.2, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=16, face="bold"),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 10, 0, 0), "pt")),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) 
  }
  
}


# -------------------------------------------------------------------------------------------------
# Arranging plots
Latitude_pollen_ensemble_simulation_ice_grid_boxplot   <- ggarrange(Latitude_pollen_ice_grid_boxplot, Latitude_ensemble_simulation_ice_grid_boxplot, 
                                                                    nrow = 1, ncol = 2, common.legend = F)

# Save the final figure as a PNG file
ggsave("result/Figure 3-Temporal changes in the latitudinal location of megabiome and ice-sheet.png", width = 13, height = 21, units = "in", dpi = 300)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 4. Figure 4：Spatiotemporal patterns of EMD between the pollen-based reconstructions and ESM-based simulation ensemble ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Map settings and customization
{
  # Define custom colors for clusters
  cols <- c("10"  = "#e31a1c",  "5"  = "#fb9a99",  "1" = "#33a02c",   "4" = "#1f78b4",  "7"  = "#ff7f00", 
            "6"   = "#fdbf6f",  "9"  = "#ffffb3",  "2" = "#6a3d9a",   "3" = "#cab2d6",  "8" = "#b15928")
  
  # Define the cluster order for consistent mapping
  brks = c("1", "2",  "3",  "4", "5",  "6",  "7",  "8",  "9",  "10")
  
  # Define labels for clusters
  labs = c("Cluster-1", "Cluster-2",  "Cluster-3",  "Cluster-4", "Cluster-5",  "Cluster-6",  "Cluster-7",  "Cluster-8",  "Cluster-9",  "Cluster-10")
  
  # Load coastline data using the Natural Earth package
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  
}


# 4.1. Spatial pattern of data-model EMD (Fig. 4a)
# -------------------------------------------------------------------------------------------------

# Load datasets
EMD_spatial_pollen_simulations_grid   <- read.csv2("data/Fig. 4 data.1-EMD between reconstructed biome score and each ensemble simulation matrix for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
EMD_spatial_pollen_simulations_grid   <- type.convert(EMD_spatial_pollen_simulations_grid, as.is = TRUE) 


# Plotting 
{
  # Compute 5th and 95th percentiles of EMD values
  Q_0.05 <- quantile(EMD_spatial_pollen_simulations_grid$EMD, probs = 0.05)
  Q_0.95 <- quantile(EMD_spatial_pollen_simulations_grid$EMD, probs = 0.95)
  
  # Cap values below the 5th percentile and above the 95th percentile
  EMD_spatial_pollen_simulations_grid$EMD[EMD_spatial_pollen_simulations_grid$EMD <= Q_0.05] <- Q_0.05
  EMD_spatial_pollen_simulations_grid$EMD[EMD_spatial_pollen_simulations_grid$EMD >= Q_0.95] <- Q_0.95
  
  # Determine the range of EMD values for color scale
  min <- range(EMD_spatial_pollen_simulations_grid$EMD)[1]
  max <- range(EMD_spatial_pollen_simulations_grid$EMD)[2]
  
  EMD_spatial_pollen_simulations_grid_map <- ggplot() +
    geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
    geom_tile(data= EMD_spatial_pollen_simulations_grid, aes(x=Longitude, y=Latitude, fill = EMD), alpha = 0.6) +  
    scale_fill_gradient(name = "EMD", low = "yellow", high = "red", limits = c(min, max), breaks = seq(0.2,0.4, by =0.05), labels = seq(0.2,0.4, by =0.05)) +  
    labs(title= "(a)", x="", y="") +
    scale_x_continuous(limits=c(-180,180), breaks= NULL, expand = c(0,0)) +
    coord_sf(ylim = c(-83.5, 83.5), expand = FALSE) +
    theme_bw() +  # Reset the plot background
    theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
    theme(axis.title.x = element_blank(), # Remove x-axis title
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.title.y = element_blank(), # Remove y-axis title
          axis.text.y = element_blank(),  # Remove y-axis text
          axis.ticks.y = element_blank()) + # Remove y-axis ticks
    theme(legend.position="bottom",
          legend.text=element_text(size=14, hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
          legend.key.size = unit(1.2, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.key.width = unit(2, "lines"),
          legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
          legend.box = "horizontal", # Adjusted to horizontal placement
          legend.justification = c(0.5, 0.5),
          legend.direction = "horizontal", # Placing legends side by side
          legend.box.spacing = unit(0.5, "cm")) + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_rect(linetype="solid", fill=NA),
          panel.background=element_rect(fill="white")) +
    theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm")) +
    guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
  
  EMD_spatial_pollen_simulations_grid_map
  
}


# 4.2 Spatial pattern of clustering (Fig. 4b)
# -------------------------------------------------------------------------------------------------

# Load datasets
EMD_clustering_spatial_pollen_simulations_grid   <- read.csv2("data/Fig. 4 data.2-Region clustering based on EMD of data-model at each timeslice since 21 ka for grid.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
EMD_clustering_spatial_pollen_simulations_grid   <- type.convert(EMD_clustering_spatial_pollen_simulations_grid, as.is = TRUE) 


# Plotting 
{
  # Convert the `Region` column to a factor with defined levels for consistent ordering
  EMD_clustering_spatial_pollen_simulations_grid$Region  <- factor(EMD_clustering_spatial_pollen_simulations_grid$Region,
                                                                   levels = c("Europe_1", "Europe_2", "Asia_1", "Asia_2", "Asia_3", "Asia_4", "North America_1", "North America_2", "North America_3",
                                                                              "Africa_1", "Africa_2", "Indopacific_1", "Indopacific_2", "South America_1", "South America_2"))           
  {
    EMD_clustering_spatial_pollen_simulations_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      labs(title= "(b)", x="", y="") +
      geom_tile(data= EMD_clustering_spatial_pollen_simulations_grid, aes(x=Longitude, y=Latitude, fill= Region), alpha = 0.6) +
      scale_fill_manual(name="Region \n (code)", 
                        labels=c("Europe 1 \n (EU 1)", "Europe 2 \n (EU 2)", "Asia 1 \n (AS 1)", "Asia 2 \n (AS 2)", "Asia 3 \n (AS 3)", "Asia 4 \n (AS 4)", 
                                 "North America 1 \n (NA 1)", "North America 2 \n (NA 2)", "North America 3 \n (NA 3)",
                                 "Africa 1 \n (AF 1)", "Africa 2 \n (AF 2)", "Indo-Pacific 1 \n (OA 1)", "Indo-Pacific 2 \n (OA 2)", "South America 1 \n (SA 1)", "South America 2 \n (SA 2)"), 
                        breaks=c("Europe_1", "Europe_2", "Asia_1", "Asia_2", "Asia_3", "Asia_4", "North America_1", "North America_2", "North America_3",
                                 "Africa_1", "Africa_2", "Indopacific_1", "Indopacific_2", "South America_1", "South America_2"),
                        values=c("Europe_1" = "#90EE90", "Europe_2" = "#33a02c", "Asia_1" = "#6a3d9a", "Asia_2" = "#fdbf6f", "Asia_3" = "#fb9a99", 
                                 "Asia_4" = "#807dba", "North America_1" = "#1f78b4", "North America_2" = "#87CEEB", "North America_3" = "#74a9cf", "Africa_1" = "#ff7f00", 
                                 "Africa_2" = "#FFBF00", "Indopacific_1" = "#C19A6B", "Indopacific_2" = "#A52A2A", "South America_2" = "#e31a1c", "South America_1" = "#fb6a4a")) +
      scale_x_continuous(limits=c(-180,180), breaks=NULL, expand = c(0,0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14, hjust = 0.5, vjust = 0.5),
            legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "horizontal", # Adjusted to horizontal placement
            legend.justification = c(0.5, 0.5),
            legend.direction = "horizontal", # Placing legends side by side
            legend.box.spacing = unit(0.5, "cm"),
            legend.margin=margin(t = -10, r = 0, b = 0, l = 0)) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm")) +
      guides(fill=guide_legend(override.aes=list(), nrow=3, byrow=FALSE))
    
    EMD_clustering_spatial_pollen_simulations_grid_map
  }  
  
}


# 4.3 Temporal pattern of data-model EMD_anomaly (Fig. 4c-f)
# -------------------------------------------------------------------------------------------------

# Load datasets
EMD_anomaly_pollen_simulations_grid          <- read.csv2("data/Fig. 4 data.3-EMD anomaly of data-model at each timeslice since 21 ka for grid by region.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
EMD_anomaly_median_pollen_simulations_grid   <- read.csv2("data/Fig. 4 data.4-Median EMD anomaly of data-model at each timeslice since 21 ka for grid by region.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
EMD_anomaly_pollen_simulations_grid          <- type.convert(EMD_anomaly_pollen_simulations_grid, as.is = TRUE) 
EMD_anomaly_median_pollen_simulations_grid   <- type.convert(EMD_anomaly_median_pollen_simulations_grid, as.is = TRUE) 

# Rename columns of the dataset
names(EMD_anomaly_median_pollen_simulations_grid) <- c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", "12", 
                                                       "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5", "21", "Region", "Region_code")
# Assign row names to the dataset
row.names(EMD_anomaly_median_pollen_simulations_grid) <- c("Global", "EU 1", "AS 1", "NA 1", "NA 2", "AS 2", "AS 3", "NA 3", "AS 4", "EU 2", "AF 1", "AF 2", "SA 1", "SA 2", "OA 1", "OA 2")


# Subset dataset by region and timeslice
{
  {
    # Convert `Timeslice` to a factor with ordered levels for consistent sorting
    EMD_anomaly_pollen_simulations_grid$Timeslice  <- factor(EMD_anomaly_pollen_simulations_grid$Timeslice,
                                                             levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", 
                                                                        "12", "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5", "21"))  
    # Convert `Region` to a factor with specified levels for grouping regions
    EMD_anomaly_pollen_simulations_grid$Region    <- factor(EMD_anomaly_pollen_simulations_grid$Region,
                                                            levels = c("Global", "EU 1", "EU 2", "AS 1", "AS 2", "AS 3", "AS 4", "NA 1", "NA 2", "NA 3", "AF 1", "AF 2", "OA 1", "OA 2", "SA 1", "SA 2"))
    
    # Subset data for specific regions
    EMD_anomaly_pollen_simulations_grid_global   <- subset(EMD_anomaly_pollen_simulations_grid, Region_code == "Global")
    EMD_anomaly_pollen_simulations_grid_cluster1 <- subset(EMD_anomaly_pollen_simulations_grid, Region_code %in% c("AS 1", "NA 1", "AS 2", "AS 4"))
    EMD_anomaly_pollen_simulations_grid_cluster2 <- subset(EMD_anomaly_pollen_simulations_grid, Region_code %in% c("AF 1", "EU 2", "NA 3", "SA 2", "OA 1"))
    EMD_anomaly_pollen_simulations_grid_cluster3 <- subset(EMD_anomaly_pollen_simulations_grid, Region_code %in% c("AS 3", "NA 2", "EU 1", "OA 2", "AF 2", "SA 1"))
    
    # Set region order for Cluster 1
    EMD_anomaly_pollen_simulations_grid_cluster1$Region_code <- factor(EMD_anomaly_pollen_simulations_grid_cluster1$Region_code,
                                                                       levels = c("AS 1", "NA 1", "AS 2", "AS 4"))
    # Set region order for Cluster 2
    EMD_anomaly_pollen_simulations_grid_cluster2$Region_code <- factor(EMD_anomaly_pollen_simulations_grid_cluster2$Region_code,
                                                                       levels = c("AF 1", "EU 2", "NA 3", "SA 2", "OA 1"))
    
    # Set region order for Cluster 3
    EMD_anomaly_pollen_simulations_grid_cluster3$Region_code <- factor(EMD_anomaly_pollen_simulations_grid_cluster3$Region_code,
                                                                       levels = c("AS 3", "NA 2", "EU 1", "OA 2", "AF 2", "SA 1"))
  }
  
  # Subset median EMD values for each region and timeslice
  {
    # Reshape `EMD_anomaly_median_pollen_simulations_grid` to long format
    EMD_timeslice_box_region_anomaly_median_long  <- pivot_longer(EMD_anomaly_median_pollen_simulations_grid[ ,-44], 
                                                                  cols = 1:43,  
                                                                  names_to = "Timeslice",  
                                                                  values_to = "EMD") 
    
    # Convert `Region_code` and `Timeslice` to factors for consistent ordering
    EMD_timeslice_box_region_anomaly_median_long$Region_code <- as.factor(EMD_timeslice_box_region_anomaly_median_long$Region_code)
    EMD_timeslice_box_region_anomaly_median_long$Timeslice   <- as.factor(EMD_timeslice_box_region_anomaly_median_long$Timeslice)
    
    # Subset median EMD data for specific regions
    EMD_timeslice_box_region_anomaly_median_long_global <-  subset(EMD_timeslice_box_region_anomaly_median_long, Region_code %in% c("Global"))
    EMD_timeslice_box_region_anomaly_median_long_cluster1 <-  subset(EMD_timeslice_box_region_anomaly_median_long, Region_code %in% c("AS 1", "NA 1", "AS 2", "AS 4"))
    EMD_timeslice_box_region_anomaly_median_long_cluster2 <-  subset(EMD_timeslice_box_region_anomaly_median_long, Region_code %in% c("AF 1", "EU 2", "NA 3", "SA 2", "OA 1"))
    EMD_timeslice_box_region_anomaly_median_long_cluster3 <-  subset(EMD_timeslice_box_region_anomaly_median_long, Region_code %in% c("AS 3", "NA 2", "EU 1", "OA 2", "AF 2", "SA 1"))
    
    # Set region order for Cluster 1 in the median dataset
    EMD_timeslice_box_region_anomaly_median_long_cluster1$Region_code <- factor(EMD_timeslice_box_region_anomaly_median_long_cluster1$Region_code,
                                                                                levels = c("AS 1", "NA 1", "AS 2", "AS 4"))
    # Set region order for Cluster 2 in the median dataset
    EMD_timeslice_box_region_anomaly_median_long_cluster2$Region_code <- factor(EMD_timeslice_box_region_anomaly_median_long_cluster2$Region_code,
                                                                                levels = c("AF 1", "EU 2", "NA 3", "SA 2", "OA 1"))
    # Set region order for Cluster 3 in the median dataset
    EMD_timeslice_box_region_anomaly_median_long_cluster3$Region_code <- factor(EMD_timeslice_box_region_anomaly_median_long_cluster3$Region_code,
                                                                                levels = c("AS 3", "NA 2", "EU 1", "OA 2", "AF 2", "SA 1"))
  }
  
}


## plot
{
  ### Global
  {
    Q05 <- quantile(EMD_anomaly_pollen_simulations_grid_global$EMD, probs = 0.05)
    Q95 <- quantile(EMD_anomaly_pollen_simulations_grid_global$EMD, probs = 0.95)
    
    EMD_anomaly_pollen_simulations_grid_global_boxplot <- ggplot(data=EMD_anomaly_pollen_simulations_grid_global, aes(x = Timeslice, y = EMD)) +
      geom_boxplot(width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", fill = "grey", alpha = 0.5, fatten = NULL) +
      geom_line(data = EMD_timeslice_box_region_anomaly_median_long_global, aes(x = Timeslice, y = EMD, group=Region_code), color="black", linewidth=0.9) +
      labs(x="Age (cal. ka BP)", y = "", title = "(c)") +
      scale_x_discrete(breaks = seq(0,21,by=1), labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_y_continuous(limits=c(Q05, Q95), breaks=seq(0.2,0.4, by= 0.1), labels=seq(0.2,0.4, by= 0.1), expand = c(0,0)) +
      facet_grid(Region_code ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
      ) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=16),
            axis.text.x=element_text(size=14, angle = 0, vjust = 0.5, hjust = 0.5, face = 'bold')) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt"), face = 'bold'),
            axis.text.y=element_text(size=14, face = 'bold')) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.8, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=14),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm"))  
    
    EMD_anomaly_pollen_simulations_grid_global_boxplot
    
  }
  
  ### Region 1
  {
    EMD_anomaly_pollen_simulations_grid_cluster1_boxplot <- ggplot(data=EMD_anomaly_pollen_simulations_grid_cluster1, aes(x = Timeslice, y = EMD)) +
      geom_boxplot(aes(fill = Region_code), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", alpha = 0.5) +
      geom_line(data =EMD_timeslice_box_region_anomaly_median_long_cluster1, aes(x = Timeslice, y = EMD, group=Region_code), color="black", linewidth=0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth=0.5, alpha = 0.9) +
      labs(x="Age (cal. ka BP)", y = "", title = "(d)") +
      scale_x_discrete(breaks = seq(0,21,by=1), labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_y_continuous(limits=c(-0.2,0.2), breaks=seq(-0.2,0.2, by= 0.1), labels=seq(-0.2,0.2, by= 0.1)) +
      scale_fill_manual(values=c("EU 1" = "#90EE90", "EU 2" = "#33a02c", "AS 1" = "#6a3d9a", "AS 2" = "#fdbf6f",  "AS 3" = "#fb9a99", 
                                 "AS 4" = "#807dba", "NA 1" = "#1f78b4", "NA 2" = "#87CEEB", "NA 3" = "#74a9cf",  "AF 1" = "#ff7f00", 
                                 "AF 2" = "#FFBF00", "OA 1" = "#C19A6B", "OA 2" = "#A52A2A", "SA 2" = "#e31a1c", "SA 1" = "#fb6a4a"), guide = "none") +
      facet_grid(Region_code ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
      ) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=16),
            axis.text.x=element_text(size=14, angle = 0, vjust = 0.5, hjust = 0.5, face = 'bold')) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt"), face = 'bold'),
            axis.text.y=element_text(size=14, face = 'bold')) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.8, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=14),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm")) 
    
    EMD_anomaly_pollen_simulations_grid_cluster1_boxplot
  }
  
  ### Region 2
  {
    EMD_anomaly_pollen_simulations_grid_cluster2_boxplot <- ggplot(data=EMD_anomaly_pollen_simulations_grid_cluster2, aes(x = Timeslice, y = EMD)) +
      geom_boxplot(aes(fill = Region_code), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", alpha = 0.5) +
      geom_line(data =EMD_timeslice_box_region_anomaly_median_long_cluster2, aes(x = Timeslice, y = EMD, group=Region_code), color="black", linewidth=0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth=0.5, alpha = 0.9) +
      labs(x="Age (cal. ka BP)", y = "", title = "(e)") +
      scale_x_discrete(breaks = seq(0,21,by=1), labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_y_continuous(limits=c(-0.2,0.2), breaks=seq(-0.2,0.2, by= 0.1), labels=seq(-0.2,0.2, by= 0.1)) +
      scale_fill_manual(values=c("EU 1" = "#90EE90", "EU 2" = "#33a02c", "AS 1" = "#6a3d9a", "AS 2" = "#fdbf6f",  "AS 3" = "#fb9a99", 
                                 "AS 4" = "#807dba", "NA 1" = "#1f78b4", "NA 2" = "#87CEEB", "NA 3" = "#74a9cf",  "AF 1" = "#ff7f00", 
                                 "AF 2" = "#FFBF00", "OA 1" = "#C19A6B", "OA 2" = "#A52A2A", "SA 2" = "#e31a1c", "SA 1" = "#fb6a4a"), guide = "none") +
      facet_grid(Region_code ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
      ) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=16),
            axis.text.x=element_text(size=14, angle = 0, vjust = 0.5, hjust = 0.5, face = 'bold')) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt"), face = 'bold'),
            axis.text.y=element_text(size=14, face = 'bold')) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.8, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=14),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm")) 
    
    EMD_anomaly_pollen_simulations_grid_cluster2_boxplot
  }
  
  ### Region 3
  {
    EMD_anomaly_pollen_simulations_grid_cluster3_boxplot <- ggplot(data=EMD_anomaly_pollen_simulations_grid_cluster3, aes(x = Timeslice, y = EMD)) +
      geom_boxplot(aes(fill = Region_code), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, coef = 0, color = "grey70", alpha = 0.5) +
      geom_line(data =EMD_timeslice_box_region_anomaly_median_long_cluster3, aes(x = Timeslice, y = EMD, group=Region_code), color="black", linewidth=0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth=0.5, alpha = 0.9) +
      labs(x="Age (cal. ka BP)", y = "", title = "(f)") +
      scale_x_discrete(breaks = seq(0,21,by=1), labels = c("0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20", "")) +
      scale_y_continuous(limits=c(-0.2,0.2), breaks=seq(-0.2,0.2, by= 0.1), labels=seq(-0.2,0.2, by= 0.1)) +
      scale_fill_manual(values=c("EU 1" = "#90EE90", "EU 2" = "#33a02c", "AS 1" = "#6a3d9a", "AS 2" = "#fdbf6f",  "AS 3" = "#fb9a99", 
                                 "AS 4" = "#807dba", "NA 1" = "#1f78b4", "NA 2" = "#87CEEB", "NA 3" = "#74a9cf",  "AF 1" = "#ff7f00", 
                                 "AF 2" = "#FFBF00", "OA 1" = "#C19A6B", "OA 2" = "#A52A2A", "SA 2" = "#e31a1c", "SA 1" = "#fb6a4a"), guide = "none") +
      facet_grid(Region_code ~ ., scales = "free_y",
                 switch = "y", # flip the facet labels along the y axis from the right side to the left
      ) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=16)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=16),
            axis.text.x=element_text(size=14, angle = 0, vjust = 0.5, hjust = 0.5, face = 'bold')) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt"), face = 'bold'),
            axis.text.y=element_text(size=14, face = 'bold')) + 
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.8, "lines")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=14),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.1,0.1,0.1),"cm")) 
    
    EMD_anomaly_pollen_simulations_grid_cluster3_boxplot
  }
  
  
}


# -------------------------------------------------------------------------------------------------

# Arrange plots
blank_plot <- ggplot() + theme_void()

EMD_anomaly_spatial_pollen_simulations_grid_final <- ggarrange(
  ggarrange(EMD_spatial_pollen_simulations_grid_map, EMD_clustering_spatial_pollen_simulations_grid_map,
            ncol = 2, nrow = 1, common.legend = FALSE, widths = c(0.48, 0.52)),
  ggarrange(blank_plot, ggarrange(EMD_anomaly_pollen_simulations_grid_global_boxplot, EMD_anomaly_pollen_simulations_grid_cluster1_boxplot,
                                  ncol = 1, nrow = 2, common.legend = FALSE, heights = c(0.275, 0.725)), 
            EMD_anomaly_pollen_simulations_grid_cluster2_boxplot, EMD_anomaly_pollen_simulations_grid_cluster3_boxplot, blank_plot,
            ncol = 5, nrow = 1, common.legend = FALSE, widths = c(0.05, 0.3, 0.3, 0.3, 0.05)),
  ncol = 1, nrow = 2, common.legend = FALSE, heights = c(0.4, 0.6))

# Save the final figure as a PNG file
ggsave("result/Figure 4-Spatiotemporal patterns of EMD between the pollen-based reconstructions and ESM-based simulation ensemble.png", width = 21, height = 15, units = "in", dpi = 300)



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 5. Appendix Figure A1：Spatial distribution and sources of fossil pollen records in the LegacyPollen 2.0 dataset ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
LegacyPollen_2.0_Overview_global    <- read.csv2("data/Appendix Fig. A1 data.1-Overview of LegacyPollen2.0 dataset.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
LegacyPollen_2.0_Overview_global    <- type.convert(LegacyPollen_2.0_Overview_global, as.is = TRUE) 


# Plotting
{
  # Load coastline data
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  
  # Reorder and set levels for the Data_Source factor
  LegacyPollen_2.0_Overview_global$Data_Source  <- factor(LegacyPollen_2.0_Overview_global$Data_Source,
                                                          levels = c("LegacyPollen 1.0",  "Neotoma",  "ACER 1.0", "Zhou et al. (2022)", "Cao et al. (2022)",  "AWI"))  
  
  {
    Legacypollen_2_source_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_point(data= LegacyPollen_2.0_Overview_global, aes(x=Longitude, y=Latitude, color= Data_Source), alpha = 0.8) +
      scale_colour_manual(name="Source", values = c("LegacyPollen 1.0"  = "black",   "Neotoma"  = "red",     "Zhou et al. (2022)" = "orange", 
                                                    "Cao et al. (2022)" = "green",   "ACER 1.0" = "#00FFFF",    "AWI"             = "blue")) +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      labs(title= "", x="", y="") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=20, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=18)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=18)) +
      theme(legend.position="bottom",
            legend.text=element_text(size=18),
            legend.title = element_text(size=20, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(color = guide_legend(override.aes=list(size=6), nrow=1, byrow=TRUE))
    
    Legacypollen_2_source_map
  }
  
  # Save the plot as a high-resolution PNG file.
  ggsave("result/Appendix Figure A1-Spatial distribution and sources of fossil pollen records in the LegacyPollen 2.0 dataset.png", width = 18, height = 10, units = "in", dpi = 300)
}  


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 6. Appendix Figure 2: Spatial patterns of megabiome distributions at 0 cal. ka BP and their agreement with modern potential natural megabiomes  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_pollen_ESM_site_0ka                   <- read.csv2("data/Appendix Fig. A2 data.1-Reconstructed and simulated megabiome at 0ka for record.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_megabiomes_site    <- read.csv2("data/Appendix Fig. A2 data.2-Modern potential natural megabiomes for record.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_megabiomes_grid    <- read.csv2("data/Appendix Fig. A2 data.3-Modern potential natural biomes for grid (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_pollen_ESM_site_0ka                   <- type.convert(Biome_pollen_ESM_site_0ka, as.is = TRUE) 
Modern_potential_natural_megabiomes_site    <- type.convert(Modern_potential_natural_megabiomes_site, as.is = TRUE) 
Modern_potential_natural_megabiomes_grid    <- type.convert(Modern_potential_natural_megabiomes_grid, as.is = TRUE) 

# Join datasets
Modern_potential_natural_pollen_ESM_megabiomes_site <- left_join(Biome_pollen_ESM_site_0ka, Modern_potential_natural_megabiomes_site[ ,c(1,7)], by = "Dataset_ID")
names(Modern_potential_natural_pollen_ESM_megabiomes_site)[c(7,14)] <- c("Reconstruction", "PNV") ## Rename columns for clarity

# Calculate the agreement
{
  # Initialize an empty result matrix
  Modern_potential_natural_pollen_ESM_megabiomes_site_agreement <- data.frame(matrix(NA, nrow=nrow(Modern_potential_natural_pollen_ESM_megabiomes_site), ncol= ncol(Modern_potential_natural_pollen_ESM_megabiomes_site)-1, byrow=TRUE))
  colnames(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement) <- names(Modern_potential_natural_pollen_ESM_megabiomes_site)[1:13]
  
  # Copy metadata from the original dataset to the agreement matrix.
  Modern_potential_natural_pollen_ESM_megabiomes_site_agreement[ ,1:6] <- Modern_potential_natural_pollen_ESM_megabiomes_site[ ,1:6]
  
  # Loop through columns 7 to 13 to check agreement
  for (i in 7:13) {
    
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement[ ,i] <- ifelse(Modern_potential_natural_pollen_ESM_megabiomes_site[ ,i] == Modern_potential_natural_pollen_ESM_megabiomes_site[ ,14], "Yes", "No")
    
  }
  
}

# Map settings and customization
{
  # Convert the Biome column to a factor with specified levels
  Modern_potential_natural_megabiomes_grid$Biome <- factor(Modern_potential_natural_megabiomes_grid$Biome,
                                                           levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
  
  # Define the colors corresponding to each megabiome category
  cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
            "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
  
  # Define the breaks (categories) for megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
  
  # Define the labels for the megabiomes
  # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
  labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
          "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
  
  # Load coastline data
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  
}


# 6.1 Spatial patterns of biome distributions at 0 cal. ka BP
# -------------------------------------------------------------------------------------------------

{
  # (1) MPI-ESM_GLAC1D
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$MPI.ESM_GLAC1D <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$MPI.ESM_GLAC1D, 
                                                                                 levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_MPI_GLAC1D_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= '(a) Megabiome distribution of MPI-ESM_GLAC1D', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= MPI.ESM_GLAC1D), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_MPI_GLAC1D_site_map
      
    }
    
  }
  
  
  # (2) MPI-ESM_ICE6G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$MPI.ESM_ICE6G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$MPI.ESM_ICE6G, 
                                                                                levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_MPI_ICE6G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of MPI-ESM_ICE6G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= MPI.ESM_ICE6G), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_MPI_ICE6G_site_map
      
    }
    
  }
  
  
  # (3) CLIMBER-X_GLAC1D
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$CLIMBER.X_GLAC1D <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$CLIMBER.X_GLAC1D, 
                                                                                   levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_CLIMBER_GLAC1D_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of CLIMBER-X_GLAC1D', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= CLIMBER.X_GLAC1D), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_CLIMBER_GLAC1D_site_map
      
    }
    
  }
  
  
  # (4) CLIMBER-X_ICE6G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$CLIMBER.X_ICE6G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$CLIMBER.X_ICE6G, 
                                                                                  levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_CLIMBER_ICE6G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of CLIMBER-X_ICE6G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= CLIMBER.X_ICE6G), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_CLIMBER_ICE6G_site_map
      
    }
    
  }
  
  
  # (5) TRACE-21K-I_ICE5G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$TRACE.21K.I_ICE5G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$TRACE.21K.I_ICE5G, 
                                                                                    levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_TRACE_I_ICE5G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of TRACE-21K-I_ICE5G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= TRACE.21K.I_ICE5G), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_TRACE_I_ICE5G_site_map
      
    }
    
  }
  
  
  # (6) TRACE-21K-II_ICE5G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site$TRACE.21K.II_ICE5G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site$TRACE.21K.II_ICE5G, 
                                                                                     levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    
    {
      Biome_0ka_TRACE_II_ICE5G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        labs(title= 'Megabiome distribution of TRACE-21K-II_ICE5G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site, aes(x=Longitude, y=Latitude, color= TRACE.21K.II_ICE5G), size = 0.6) +
        scale_colour_manual(name="Megabiome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=4), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_TRACE_II_ICE5G_site_map
      
    }
    
  }
  
}


# 6.2 Simulated megabiome agreement with modern potential natural biomes
# -------------------------------------------------------------------------------------------------

{
  # (1) MPI-ESM_GLAC1D
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$MPI.ESM_GLAC1D <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$MPI.ESM_GLAC1D, 
                                                                                           levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_MPI_GLAC1D_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= '(b) Megabiome agreement of MPI-ESM_GLAC1D', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= MPI.ESM_GLAC1D), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_MPI_GLAC1D_site_map
      
    }
    
  }
  
  
  # (2) MPI-ESM_ICE6G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$MPI.ESM_ICE6G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$MPI.ESM_ICE6G, 
                                                                                          levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_MPI_ICE6G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of MPI-ESM_ICE6G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= MPI.ESM_ICE6G), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_MPI_ICE6G_site_map
      
    }
    
  }
  
  
  # (3) CLIMBER-X_GLAC1D
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$CLIMBER.X_GLAC1D <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$CLIMBER.X_GLAC1D, 
                                                                                             levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_CLIMBER_GLAC1D_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of CLIMBER-X_GLAC1D', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= CLIMBER.X_GLAC1D), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE))  
      
      Biome_agreement_0ka_CLIMBER_GLAC1D_site_map
      
    }
    
  }
  
  
  # (4) CLIMBER-X_ICE6G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$CLIMBER.X_ICE6G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$CLIMBER.X_ICE6G, 
                                                                                            levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_CLIMBER_ICE6G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of CLIMBER-X_ICE6G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= CLIMBER.X_ICE6G), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_CLIMBER_ICE6G_site_map
      
    }
    
  }
  
  
  # (5) TRACE-21K-I_ICE5G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$TRACE.21K.I_ICE5G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$TRACE.21K.I_ICE5G, 
                                                                                              levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_TRACE_I_ICE5G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of TRACE-21K-I_ICE5G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= TRACE.21K.I_ICE5G), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE))  
      
      Biome_agreement_0ka_TRACE_I_ICE5G_site_map
      
    }
    
  }
  
  
  # (6) TRACE-21K-II_ICE5G
  {
    # Convert biome to a factor
    Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$TRACE.21K.II_ICE5G <- factor(Modern_potential_natural_pollen_ESM_megabiomes_site_agreement$TRACE.21K.II_ICE5G, 
                                                                                               levels = c("Yes", "No"))
    
    {
      Biome_agreement_0ka_TRACE_II_ICE5G_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_megabiomes_grid, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Megabiome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= 'Megabiome agreement of TRACE-21K-II_ICE5G', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Modern_potential_natural_pollen_ESM_megabiomes_site_agreement, aes(x=Longitude, y=Latitude, color= TRACE.21K.II_ICE5G), size = 0.6) +
        scale_colour_manual(name="Status", values=c("Yes" = "yellow", "No" = "red"), labels= c("Yes" = "Agree", "No" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
        theme(axis.title.x = element_blank(), # Remove x-axis title
              axis.text.x = element_blank(),  # Remove x-axis text
              axis.ticks.x = element_blank(), # Remove x-axis ticks
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text.y = element_blank(),  # Remove y-axis text
              axis.ticks.y = element_blank()) + # Remove y-axis ticks
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
        guides(color=guide_legend(override.aes=list(size=4), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_TRACE_II_ICE5G_site_map
      
    }
    
  }
  
}


# -------------------------------------------------------------------------------------------------
#Arranging plots
Biome_0ka_agreement_site_all_models_map   <- ggarrange(ggarrange(Biome_0ka_MPI_GLAC1D_site_map, Biome_0ka_MPI_ICE6G_site_map, Biome_0ka_CLIMBER_GLAC1D_site_map, 
                                                                 Biome_0ka_CLIMBER_ICE6G_site_map, Biome_0ka_TRACE_I_ICE5G_site_map, Biome_0ka_TRACE_II_ICE5G_site_map, 
                                                                 nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                       ggarrange(Biome_agreement_0ka_MPI_GLAC1D_site_map, Biome_agreement_0ka_MPI_ICE6G_site_map, Biome_agreement_0ka_CLIMBER_GLAC1D_site_map, 
                                                                 Biome_agreement_0ka_CLIMBER_ICE6G_site_map, Biome_agreement_0ka_TRACE_I_ICE5G_site_map, Biome_agreement_0ka_TRACE_II_ICE5G_site_map,
                                                                 nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                       nrow = 2, ncol = 1, common.legend = F, heights = c(0.51,0.49))

# Save the final figure as a PNG file
ggsave("result/Appendix Figure A2-Spatial patterns of megabiome distributions at 0 ka and their agreement with modern potential natural megabiomes.png", width = 13, height = 21, units = "in", dpi = 300)



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 7. Appendix Figure 3：Differences in bioclimatic variables between ESM-based simulations at 0ka and observations ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Modern_bioclimatic_variables_CRU_model_difference_grid     <- read.csv2("data/Appendix Fig. A3 data.1-Differences in bioclimatic variables between ESM-based simulations at 0 cal. ka BP and observations.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Modern_bioclimatic_variables_CRU_model_difference_grid     <- type.convert(Modern_bioclimatic_variables_CRU_model_difference_grid, as.is = TRUE) 

# Load coastline data
coastline <- ne_coastline(scale = "medium", returnclass = "sf")


# 7.1 BIO10 = Mean Temperature of Warmest Quarter
# -------------------------------------------------------------------------------------------------
{
  # Calculate the 5th and 95th percentiles for variable
  Q5  <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10, probs = 0.05)
  Q95 <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10, probs = 0.95)
  
  # Cap values outside the 5th and 95th percentiles to those thresholds
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10 <= Q5] <- Q5
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO10 >= Q95] <- Q95
  
  ## (1) MPI-ESM_GLAC1D
  {
    BIO10_difference_MPI_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold("(a) ") * bold(T[Warm]) * ": MPI-ESM_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_MPI_GLAC1D_grid_map
    
  }
  
  ## (2) MPI-ESM_ICE6G
  {
    BIO10_difference_MPI_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +   
      labs(title = expression(bold(T[Warm]) * ": MPI-ESM_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_MPI_ICE6G_grid_map
    
  }
  
  ## (3) CLIMBER-X_GLAC1D
  {
    BIO10_difference_CLIMBER_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +   
      labs(title = expression(bold(T[Warm]) * ": CLIMBER-X_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_CLIMBER_GLAC1D_grid_map
    
  }
  
  ## (4) CLIMBER-X_ICE6G
  {
    BIO10_difference_CLIMBER_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Warm]) * ": CLIMBER-X_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_CLIMBER_ICE6G_grid_map
    
  }
  
  ## (5) TRACE-21K-I_ICE5G
  {
    BIO10_difference_TRACE_I_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_I - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Warm]) * ": TRACE-21K-I_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_TRACE_I_grid_map
    
  }
  
  ## (6) TRACE-21K-II_ICE5G
  {
    BIO10_difference_TRACE_II_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_II - CRU"), aes(x=Longitude, y=Latitude, fill = BIO10)) +
      scale_fill_gradient2(name = "BIO 10", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,4, by =4), Q95), labels = c(seq(-8,4, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Warm]) * ": TRACE-21K-II_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO10_difference_TRACE_II_grid_map
    
  }
  
}


# 7.2 BIO11 = Mean Temperature of Coldest Quarter
# -------------------------------------------------------------------------------------------------
{
  # Calculate the 5th and 95th percentiles for variable
  Q5  <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11, probs = 0.05)
  Q95 <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11, probs = 0.95)
  
  # Cap values outside the 5th and 95th percentiles to those thresholds
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11 <= Q5] <- Q5
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO11 >= Q95] <- Q95
  
  ## (1) MPI-ESM_GLAC1D
  {
    BIO11_difference_MPI_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold("(b) ") * bold(T[Cold]) * ": MPI-ESM_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_MPI_GLAC1D_grid_map
    
  }
  
  ## (2) MPI-ESM_ICE6G
  {
    BIO11_difference_MPI_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +   
      labs(title = expression(bold(T[Cold]) * ": MPI-ESM_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_MPI_ICE6G_grid_map
    
  }
  
  ## (3) CLIMBER-X_GLAC1D
  {
    BIO11_difference_CLIMBER_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +   
      labs(title = expression(bold(T[Cold]) * ": CLIMBER-X_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_CLIMBER_GLAC1D_grid_map
    
  }
  
  ## (4) CLIMBER-X_ICE6G
  {
    BIO11_difference_CLIMBER_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Cold]) * ": CLIMBER-X_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_CLIMBER_ICE6G_grid_map
    
  }
  
  ## (5) TRACE-21K-I_ICE5G
  {
    BIO11_difference_TRACE_I_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_I - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Cold]) * ": TRACE-21K-I_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_TRACE_I_grid_map
    
  }
  
  ## (6) TRACE-21K-II_ICE5G
  {
    BIO11_difference_TRACE_II_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_II - CRU"), aes(x=Longitude, y=Latitude, fill = BIO11)) +
      scale_fill_gradient2(name = "BIO 11", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-8,8, by =4), Q95), labels = c(seq(-8,8, by =4), expression("        [°C]"))) +  
      labs(title = expression(bold(T[Cold]) * ": TRACE-21K-II_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO11_difference_TRACE_II_grid_map
    
  }
  
}


# 7.3 BIO18 = Precipitation of Warmest Quarter
# -------------------------------------------------------------------------------------------------
{
  # Calculate the 5th and 95th percentiles for variable
  Q5  <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18, probs = 0.05)
  Q95 <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18, probs = 0.95)
  
  # Cap values outside the 5th and 95th percentiles to those thresholds
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18 <= Q5] <- Q5
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO18 >= Q95] <- Q95
  
  ## (1) MPI-ESM_GLAC1D
  {
    BIO18_difference_MPI_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold("(c) ") * bold(P[Warm]) * ": MPI-ESM_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_MPI_GLAC1D_grid_map
    
  }
  
  ## (2) MPI-ESM_ICE6G
  {
    BIO18_difference_MPI_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +   
      labs(title = expression(bold(P[Warm]) * ": MPI-ESM_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_MPI_ICE6G_grid_map
    
  }
  
  ## (3) CLIMBER-X_GLAC1D
  {
    BIO18_difference_CLIMBER_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +   
      labs(title = expression(bold(P[Warm]) * ": CLIMBER-X_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_CLIMBER_GLAC1D_grid_map
    
  }
  
  ## (4) CLIMBER-X_ICE6G
  {
    BIO18_difference_CLIMBER_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Warm]) * ": CLIMBER-X_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_CLIMBER_ICE6G_grid_map
    
  }
  
  ## (5) TRACE-21K-I_ICE5G
  {
    BIO18_difference_TRACE_I_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_I - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Warm]) * ": TRACE-21K-I_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_TRACE_I_grid_map
    
  }
  
  ## (6) TRACE-21K-II_ICE5G
  {
    BIO18_difference_TRACE_II_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_II - CRU"), aes(x=Longitude, y=Latitude, fill = BIO18)) +
      scale_fill_gradient2(name = "BIO 18", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,250, by =100), Q95), labels = c(seq(-200,250, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Warm]) * ": TRACE-21K-II_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO18_difference_TRACE_II_grid_map
    
  }
  
}


# 7.4 BIO19 = Precipitation of Coldest Quarter
# -------------------------------------------------------------------------------------------------
{
  # Calculate the 5th and 95th percentiles for variable
  Q5  <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19, probs = 0.05)
  Q95 <- quantile(Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19, probs = 0.95)
  
  # Cap values outside the 5th and 95th percentiles to those thresholds
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19 <= Q5] <- Q5
  Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19[Modern_bioclimatic_variables_CRU_model_difference_grid$BIO19 >= Q95] <- Q95
  
  ## (1) MPI-ESM_GLAC1D
  {
    BIO19_difference_MPI_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold("(d) ") * bold(P[Cold]) * ": MPI-ESM_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_MPI_GLAC1D_grid_map
    
  }
  
  ## (2) MPI-ESM_ICE6G
  {
    BIO19_difference_MPI_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "MPI-ESM_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +   
      labs(title = expression(bold(P[Cold]) * ": MPI-ESM_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_MPI_ICE6G_grid_map
    
  }
  
  ## (3) CLIMBER-X_GLAC1D
  {
    BIO19_difference_CLIMBER_GLAC1D_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_GLAC1D - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +   
      labs(title = expression(bold(P[Cold]) * ": CLIMBER-X_GLAC1D - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_CLIMBER_GLAC1D_grid_map
    
  }
  
  ## (4) CLIMBER-X_ICE6G
  {
    BIO19_difference_CLIMBER_ICE6G_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "CLIMBER_ICE6G - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Cold]) * ": CLIMBER-X_ICE6G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_CLIMBER_ICE6G_grid_map
    
  }
  
  ## (5) TRACE-21K-I_ICE5G
  {
    BIO19_difference_TRACE_I_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_I - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Cold]) * ": TRACE-21K-I_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_TRACE_I_grid_map
    
  }
  
  ## (6) TRACE-21K-II_ICE5G
  {
    BIO19_difference_TRACE_II_grid_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= subset(Modern_bioclimatic_variables_CRU_model_difference_grid, Source == "TRACE_II - CRU"), aes(x=Longitude, y=Latitude, fill = BIO19)) +
      scale_fill_gradient2(name = "BIO 19", low = "#1E90FF", mid = "white", high = "red", midpoint = 0, limits = c(Q5, Q95), breaks = c(seq(-200,150, by =100), Q95), labels = c(seq(-200,150, by =100), expression("        [mm]"))) +  
      labs(title = expression(bold(P[Cold]) * ": TRACE-21K-II_ICE5G - Observation"), x = "", y = "") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0)) +
      theme(axis.title.x = element_blank(), # Remove x-axis title
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(), # Remove x-axis ticks
            axis.title.y = element_blank(), # Remove y-axis title
            axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank()) + # Remove y-axis ticks
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill = guide_colorbar(direction = "horizontal", barwidth = 21, barheight = 0.8, title.position = "left"))
    
    BIO19_difference_TRACE_II_grid_map
    
  }
  
}


# -------------------------------------------------------------------------------------------------
#Arranging plots
BIO10_BIO11_difference_all_models_grid_map   <- ggarrange(ggarrange(BIO10_difference_MPI_GLAC1D_grid_map, BIO10_difference_MPI_ICE6G_grid_map, BIO10_difference_CLIMBER_GLAC1D_grid_map, 
                                                                    BIO10_difference_CLIMBER_ICE6G_grid_map, BIO10_difference_TRACE_I_grid_map, BIO10_difference_TRACE_II_grid_map, 
                                                                    nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                          ggarrange(BIO11_difference_MPI_GLAC1D_grid_map, BIO11_difference_MPI_ICE6G_grid_map, BIO11_difference_CLIMBER_GLAC1D_grid_map, 
                                                                    BIO11_difference_CLIMBER_ICE6G_grid_map, BIO11_difference_TRACE_I_grid_map, BIO11_difference_TRACE_II_grid_map,
                                                                    nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                          nrow = 2, ncol = 1, common.legend = F)

# Save the final figure as a PNG file
ggsave("result/Appendix Figure A3-1-Differences in bioclimatic variables between ESM-based simulations at 0 ka and observations_BIO 10 and BIO 11.png", width = 13, height = 20, units = "in", dpi = 300)


BIO18_BIO19_difference_all_models_grid_map   <- ggarrange(ggarrange(BIO18_difference_MPI_GLAC1D_grid_map, BIO18_difference_MPI_ICE6G_grid_map, BIO18_difference_CLIMBER_GLAC1D_grid_map, 
                                                                    BIO18_difference_CLIMBER_ICE6G_grid_map, BIO18_difference_TRACE_I_grid_map, BIO18_difference_TRACE_II_grid_map,
                                                                    nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                          ggarrange(BIO19_difference_MPI_GLAC1D_grid_map, BIO19_difference_MPI_ICE6G_grid_map, BIO19_difference_CLIMBER_GLAC1D_grid_map, 
                                                                    BIO19_difference_CLIMBER_ICE6G_grid_map, BIO19_difference_TRACE_I_grid_map, BIO19_difference_TRACE_II_grid_map,
                                                                    nrow = 3, ncol = 2, common.legend = T, legend = "bottom"),
                                                          nrow = 2, ncol = 1, common.legend = F)

# Save the final figure as a PNG file
ggsave("result/Appendix Figure A3-2-Differences in bioclimatic variables between ESM-based simulations at 0 ka and observations_BIO 18 and BIO 19.png", width = 13, height = 20, units = "in", dpi = 300)


