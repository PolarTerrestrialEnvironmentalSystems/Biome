
###################################
# R version 4.4.1                 #
# Operating System: Windows 10    #
# Code for statistical analysis   #
# Supplement to: Li, C., Dallmeyer, A., Ni, J., Chevalier, M., Willeit, M., Andreev, A. A., Cao, X., Schild, L., Heim, B., and Herzschuh, U.: Global biome changes over the last 21,000 years inferred from model-data comparisons, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-1862, 2024. #
# Contact: Chenzhi Li (chenzhi.li@awi.de)[Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Potsdam, Germany 2024] #
###################################

# Statistical analysis list:
## 1. Statistical analysis 1: Select the pollen sample of the target timeslice from the time-nearest sample (Asia as example)
## 2. Statistical analysis 2: Calculate the biome affinity score of pollen-based reconstructions for site (Asia as example)
## 3. Statistical analysis 3: Assign ice-sheet and biome of site to grid-cell 
## 4. Statistical analysis 4: Validation of model climate bias
## 5. Statistical analysis 5: Space Constrained Clustering
## 6. Statistical analysis 6: EMD between pollen-based reconstructions and ESM-based simulations

# Note: Biomes and their abbreviations and codes 
## 1 - Tropical forest (TRFO); 2 - Subtropical forest (WTFO); 3 - Temperate forest (TEFO); 4 - Boreal forest (BOFO); 
## 5 - (Warm) savanna and dry woodland (SAVA); 6 - Grassland and dry shrubland (STEP); 7 - (Warm) desert (DESE); 8 - Tundra and polar desert (TUND)

# Note: Bioclimatic variables
## BIO10 = Mean Temperature of Warmest Quarter; BIO11 = Mean Temperature of Coldest Quarter; 
## BIO18 = Precipitation of Warmest Quarter; BIO19 = Precipitation of Coldest Quarter

# choose directory:
setwd("~/Supplementary code/Statistical analyses")  # select the folder "Supplementary code" from your directory

# install packages if not installed
# install.packages(c("paleotools", "dplyr", "data.table", "ggpubr", "TSclust", "spdep", "adespatial", "factoextra"))

# loading packages
library(paleotools)
library(dplyr)
library(data.table)
library(ggpubr)
library(TSclust)
library(spdep)
library(adespatial)
library(factoextra)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 1. Statistical analysis 1: Select the pollen sample of the target timeslice from the time-nearest sample (Asia as example) ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Pollen_percentages_sample_site_Asia_0_21ka   <- read.csv2("data/Statistical analysis 1_data.1-Pollen percentages of sample for site in Aisa_0-21ka.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Pollen_percentages_sample_site_Asia_0_21ka   <- type.convert(Pollen_percentages_sample_site_Asia_0_21ka, as.is = TRUE) 

# Calculation
{
  # Define a sequence of timeslices from 0 to 21 (in steps of 0.5)
  timeslice_list <- seq(0,21, by = 0.5)
  
  # Extract unique Dataset_IDs from the dataset
  DatasetID_list <- unique(Pollen_percentages_sample_site_Asia_0_21ka$Dataset_ID)
  
  # Initialize an empty data container to store the final results
  Pollen_percentages_timeslice_site_Asia_0_21ka <- NULL
  
  # Loop through each unique Dataset_ID
  for (ID in DatasetID_list) {
    
    # Print progress message for the current Dataset_ID
    print(paste0("+++++ Dataset_ID ",ID," | (",which(DatasetID_list == ID),"/",length(DatasetID_list),") +++++"))
    
    # Subset the data for the current Dataset_ID
    subset_site <- subset(Pollen_percentages_sample_site_Asia_0_21ka, Dataset_ID == ID) 
    
    # Create a data frame to store intermediate results for the current site
    Pollen_percentages_timeslice_site_0_21ka <- data.frame(matrix(NA, nrow=length(timeslice_list), ncol=ncol(Pollen_percentages_sample_site_Asia_0_21ka)-5, byrow=TRUE))
    colnames(Pollen_percentages_timeslice_site_0_21ka) <- c(names(Pollen_percentages_sample_site_Asia_0_21ka)[1:14], "Timeslice_ka", names(Pollen_percentages_sample_site_Asia_0_21ka)[21:ncol(Pollen_percentages_sample_site_Asia_0_21ka)])
    
    # Replicate metadata for each timeslice and assign timeslices
    Pollen_percentages_timeslice_site_0_21ka[ ,1:14] <- unique(subset_site[ ,1:14])[rep(seq_len(nrow(unique(subset_site[ ,1:14]))), each=length(timeslice_list)), ]
    Pollen_percentages_timeslice_site_0_21ka[ ,15]   <- timeslice_list
    
    # Loop through each timeslice
    for(i in 1:length(timeslice_list)) {
      
      # Subset the data to include rows within ±0.25 ka of the current timeslice
      subset_site_sample <- subset(subset_site, meanAgeBP >= timeslice_list[i] - 0.25 & meanAgeBP <= timeslice_list[i] + 0.25)
      
      if(nrow(subset_site_sample) == 0) {
        
        # If no data matches the current timeslice, set the values to NA
        Pollen_percentages_timeslice_site_0_21ka[i, 16:ncol(Pollen_percentages_timeslice_site_0_21ka)] <- NA
        
      }else{
        
        if(ncol(subset_site_sample) == 1) {
          
          # If only one row matches, assign the corresponding values
          Pollen_percentages_timeslice_site_0_21ka[i, 16:ncol(Pollen_percentages_timeslice_site_0_21ka)] <- subset_site_sample[1,21:ncol(subset_site_sample)]
          
        }else{
          
          # If multiple rows match, find the one closest to the current timeslice and assign its values
          Pollen_percentages_timeslice_site_0_21ka[i, 16:ncol(Pollen_percentages_timeslice_site_0_21ka)] <- subset_site_sample[which.min(abs(subset_site_sample$meanAgeBP - timeslice_list[i])),21:ncol(subset_site_sample)]
          
        }
        
      }
      
      
    }
    
    
    # Append the processed results for the current site to the main result container
    Pollen_percentages_timeslice_site_Asia_0_21ka <- rbind(Pollen_percentages_timeslice_site_Asia_0_21ka, Pollen_percentages_timeslice_site_0_21ka)
    
    # Remove rows with missing values (NA) from the main result container
    Pollen_percentages_timeslice_site_Asia_0_21ka <- Pollen_percentages_timeslice_site_Asia_0_21ka[complete.cases(Pollen_percentages_timeslice_site_Asia_0_21ka[, 15:ncol(Pollen_percentages_timeslice_site_Asia_0_21ka)]), ]
  }
  
  # Save the final results as a CSV file
  write.csv(Pollen_percentages_timeslice_site_Asia_0_21ka, file="result/Statistical analysis 1_result.1-Pollen percentages of timeslice for site in Aisa_0-21ka.csv", row.names=FALSE)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 2. Statistical analysis 2: Calculate the biome affinity score of pollen-based reconstructions for site (Asia as example) ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load pollen data
pollen   <- read.csv2("data/Statistical analysis 2_data.1-Pollen percentages of timeslice for site in Aisa_0-21ka.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
pollen   <- type.convert(pollen, as.is = TRUE) # Convert data types
pollen   <- pollen[, which(names(pollen) != "Indeterminable")]  # Remove the "Indeterminable" taxa column

# Load Taxa-PFTs data
pft              <- read.csv2("data/Statistical analysis 2_data.2-Taxa-PFTs_Asia.csv", sep=";", stringsAsFactors=FALSE, check.names=FALSE)     # change for other regions 
pft              <- type.convert(pft, as.is = TRUE) # Convert data types
colnames(pft)[1] <- "taxa" # Rename the first column to "taxa"

# Load PFTs-Biomes data
bio              <- read.csv2("data/Statistical analysis 2_data.3-PFTs-MegaBiomes_Asia.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep=";")     # change for other regions 
bio              <- type.convert(bio, as.is = TRUE) # Convert data types
colnames(bio)[1] <- "biome" # Rename the first column to "biome"

# Load Taxa-variable data
variable              <- read.csv2("data/Statistical analysis 2_data.4-Taxa-variable_Asia.csv", dec=".", sep=";", stringsAsFactors=FALSE)   # change for other regions
variable              <- type.convert(variable, as.is = TRUE) # Convert data types
colnames(variable)[1] <- "taxa" # Rename the first column to "taxa"
variable              <- variable[-which(variable$taxa == "Indeterminable"), ]  # Remove the "Indeterminable" taxa row

# Helper function to find the nth highest value in a vector: https://stackoverflow.com/questions/2453326/fastest-way-to-find-second-third-highest-lowest-value-in-vector-or-column
maxN <- function(x, N=1){
  len <- length(x)
  if(N > len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
} # This function identifies the nth largest value in a numeric vector

# Define a custom negation operator for "not in"
`%notin%` <- Negate(`%in%`)


# -------------------------------------------------------------------------------------------------
# Calculate biome scores
{
  # Get the list of unique Dataset_IDs from the pollen data
  list_ID <- unique(pollen$Dataset_ID)
  
  # Initialize a result container for all sites
  resulttable_complete <- NULL
  
  # Loop through each Dataset_ID to process individual sites
  for (ID in list_ID) { #add loop for lots of sites
    
    print(paste0("+++++ SITE ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
    
    # Subset pollen data for the current site
    subset_site <- subset(pollen, Dataset_ID == ID) 
    
    # Get the number of pollen taxa and their names
    pollen.no   <- nrow(variable)
    pollennames <- variable$taxa
    
    # Split site data into metadata (first 15 columns) and taxa-specific data
    subset_site_info <- subset_site[ ,c(1:15)]
    subset_site_taxa <- subset_site[ ,names(subset_site) %in% pollennames]
    
    # Calibrate pollen percentages based on thresholds and weights
    taxa_corrected <- data.frame(matrix(NA, nrow=dim(subset_site_taxa)[1], ncol=dim(subset_site_taxa)[2], byrow=TRUE))
    colnames(taxa_corrected) <- names(subset_site_taxa)
    
    # Loop for threshold and weight calibration
    for(k in 1:pollen.no){
      
      # Extract pollen percentage for the current taxon
      a <- as.vector(unlist(subset_site_taxa[pollennames[k]])) # subtract pollen percentage 
      
      # Set pollen percentages below the threshold to 0
      a[(a < variable[k,"threshold"])] <- 0  
      
      # Multiply the adjusted percentage by the weight
      taxa_corrected[ ,pollennames[k]] <- a * variable[k,"weight"] 
      
    } 
    
    # Calculate biome scores sample by sample
    sample.no <- nrow(subset_site)  # Number of samples
    info.no   <- ncol(subset_site_info) # Number of metadata columns
    
    biome.no   <- nrow(bio) # Number of biomes
    biomenames <- bio$biome # Biome names
    
    # Initialize result matrix for the current site
    res             <- data.frame(matrix(NA, nrow=sample.no, ncol=info.no+biome.no+4, byrow=TRUE))
    res[,1:info.no] <- subset_site_info
    colnames(res)   <- c(names(subset_site_info), biomenames, "Best_Biome", "Best_AffinityScore", "Second_Best_Biome", "Second_Best_AffinityScore")  ##add "Best_Column", "Second_Best_Column", 
    
    # Loop through each sample
    for(i in 1:sample.no){
      
      # Loop through each biome
      for(j in 1:biome.no){
        
        # Select the current biome and identify its PFTs
        biome1 <- biomenames[j]  # Current biome name                       
        pft1   <- colnames(bio)[which(bio[j, ] == 1)]  # Identify PFTs present in this biome
        
        # Identify taxa that belong to the selected PFTs
        if(length(pft1) > 1){
          
          # If there are multiple PFTs, find all taxa contributing to these PFTs
          pollen_pft <- pollennames[rowSums(pft[ ,which(names(pft) %in% pft1)]) > 0]
          
        }else{
          
          # If there is only one PFT, find taxa contributing to it
          pollen_pft <- pollennames[pft[ ,which(names(pft) %in% pft1)] > 0]
          
        }  
        
        # Calculate the affinity score for the current biome
        res[i,biome1] <- sum(sqrt(taxa_corrected[i,pollen_pft])) 
        
      }
      
      
      ### Check if any biomes have the same affinity score
      
      if(any(duplicated(as.numeric(res[i,biomenames])))){
        
        # Extract duplicated affinity scores
        dups <- res[i,biomenames][duplicated(as.numeric(res[i,biomenames])) | duplicated(as.numeric(res[i,biomenames]), fromLast=TRUE)]
        
        # Case 1: Multiple pairs of duplicates with different affinity scores
        
        if(length(unique(as.numeric(dups))) > 1){
          
          # Identify the duplicate with the highest affinity score
          higher_dups <- dups[ ,dups == sort(unique(as.numeric(dups)), decreasing=TRUE)[1]]
          
          # Find PFT counts for these duplicates
          pft_dups <- bio[which(bio$biome %in% names(higher_dups)), ]  
          pft_dups_counts <- cbind.data.frame(pft_dups[,1],rowSums(pft_dups[,-1]))
          colnames(pft_dups_counts) <- c("biome","count")
          pft_dups_counts <- pft_dups_counts[order(pft_dups_counts$count), ]
          
          # Case 1.1: The duplicate biome is the one with the highest affinity score
          
          if(unique(pft_dups_counts$biome %in% names(res[i,biomenames][which(res[i,biomenames] == sort(unique(as.numeric(res[i,biomenames])), decreasing=TRUE)[1])]))){
            
            res[i,"Best_Biome"]                <- names(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-1]])], N=1))   # duplicates with the higher number of PFT are excluded from the list (biomes with less PFT gets priority)
            res[i,"Best_AffinityScore"]        <- as.vector(unlist(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-1]])], N=1))) 
            res[i,"Second_Best_Biome"]         <- names(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-2]])], N=1))           # duplicate with the lower number of PFT is excluded from the list (as it is already chosen in the best biome)
            res[i,"Second_Best_AffinityScore"] <- as.vector(unlist(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-2]])], N=1))) 
            
          }else{
            
            # Case 1.2: The duplicate biome is not the one with the highest affinity score
            
            res[i,"Best_Biome"]                <- names(maxN(res[i,biomenames], N=1)) 
            res[i,"Best_AffinityScore"]        <- as.vector(unlist(maxN(res[i,biomenames], N=1))) 
            res[i,"Second_Best_Biome"]         <- names(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-1]])], N=2))  # duplicate with the higher number of PFT is excluded from the list (biomes with less PFT gets priority)
            res[i,"Second_Best_AffinityScore"] <- as.vector(unlist(maxN(res[i,(biomenames[biomenames %notin% pft_dups_counts$biome[-1]])], N=2))) 
            
          }
          
          
        }else{
          
          # Case 2: Only one pair of duplicates
          pft_dups <- bio[which(bio$biome %in% names(dups)), ]  
          pft_dups_counts <- cbind.data.frame(pft_dups[,1],rowSums(pft_dups[,-1]))
          colnames(pft_dups_counts) <- c("biome","count")
          pft_dups_counts <- pft_dups_counts[order(pft_dups_counts$count), ]
          
          # Case 2.1: The duplicate biome is the one with the highest affinity score
          
          if(unique(pft_dups_counts$biome %in% names(res[i,biomenames][which(res[i,biomenames] == sort(unique(as.numeric(res[i,biomenames])), decreasing=TRUE)[1])]))){
            
            res[i,"Best_Biome"]                <- biomenames[biomenames %in% pft_dups_counts$biome[1]]  # duplicate with the lowest number of PFT is chosen
            res[i,"Best_AffinityScore"]        <- as.vector(res[i,(biomenames[biomenames %in% pft_dups_counts$biome[1]])]) 
            res[i,"Second_Best_Biome"]         <- biomenames[biomenames %in% pft_dups_counts$biome[2]]  # duplicate with the second lowest number of PFT is chosen
            res[i,"Second_Best_AffinityScore"] <- as.vector(res[i,(biomenames[biomenames %in% pft_dups_counts$biome[2]])]) 
            
          }else{
            
            # Case 2.2: The duplicate biome is not the one with the highest affinity score
            
            res[i,"Best_Biome"]                <- names(maxN(res[i,biomenames], N=1)) 
            res[i,"Best_AffinityScore"]        <- as.vector(unlist(maxN(res[i,biomenames], N=1))) 
            res[i,"Second_Best_Biome"]         <- ifelse(names(maxN(res[i,biomenames], N=2)) %in% pft_dups_counts$biome, biomenames[biomenames %in% pft_dups_counts$biome[1]], names(maxN(res[i,biomenames], N=2)))  # duplicate with the lowest number of PFT is chosen
            res[i,"Second_Best_AffinityScore"] <- ifelse(names(maxN(res[i,biomenames], N=2)) %in% pft_dups_counts$biome, as.vector(res[i,(biomenames[biomenames %in% pft_dups_counts$biome[1]])]), as.vector(unlist(maxN(res[i,biomenames], N=2))))
            
          }
          
        }
        
        
      }else{
        
        # If no duplicate affinity scores
        
        res[i,"Best_Biome"]                <- names(maxN(res[i,biomenames], N=1)) 
        res[i,"Best_AffinityScore"]        <- as.vector(unlist(maxN(res[i,biomenames], N=1))) 
        res[i,"Second_Best_Biome"]         <- names(maxN(res[i,biomenames], N=2))
        res[i,"Second_Best_AffinityScore"] <- as.vector(unlist(maxN(res[i,biomenames], N=2))) 
        
      }
      
    } 
    
    # Append the results for the current site to the overall result container
    resulttable_complete <- rbind(resulttable_complete, res)                
    
  }
  
  # Save the final results as a CSV file
  write.csv(resulttable_complete, file="result/Statistical analysis 2_result.1-Biome affinity score of pollen-based reconstructions in Asia for site_0-21ka.csv", row.names=FALSE)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 3. Statistical analysis 3: Assign ice-sheet and biome of site to grid-cell ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Coordinate_list_T31                                   <- read.csv2("data/Statistical analysis 3_data.1-Coordinate list of T31 resolution.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
ICE5G_ice_sheet_data                                  <- read.csv2("data/Statistical analysis 3_data.2-Ice sheet of ICE5G since 21 ka BP.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
ICE6G_GLAC1D_ice_sheet_grid                           <- read.csv2("data/Statistical analysis 3_data.3-Ice sheet of ICE6G and GLAC1D since 26 ka BP.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Reconstructed_simulated_biome_timeslice_site_0_21ka   <- read.csv2("data/Statistical analysis 3_data.4-Reconstructed and simulated biome at timeslice for site since 21 ka BP.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Coordinate_list_T31                                   <- type.convert(Coordinate_list_T31, as.is = TRUE) 
ICE5G_ice_sheet_data                                  <- type.convert(ICE5G_ice_sheet_data, as.is = TRUE) 
ICE6G_GLAC1D_ice_sheet_grid                           <- type.convert(ICE6G_GLAC1D_ice_sheet_grid, as.is = TRUE) 
Reconstructed_simulated_biome_timeslice_site_0_21ka   <- type.convert(Reconstructed_simulated_biome_timeslice_site_0_21ka, as.is = TRUE) 

# 3.1 Ice-sheet ensemble of grid-cell
# -------------------------------------------------------------------------------------------------
{
  # (1) Assign ICE5G to grid-cell (3.75° resolution)
  {
    # Determine the number of unique longitudes and latitudes
    Longitude_number <- length(unique(Coordinate_list_T31$Longitude))
    Latitude_number  <- length(unique(Coordinate_list_T31$Latitude))
    
    # Define the timeslices for analysis
    timeslice_list <- c(seq(0,17,by=0.5), seq(18,21,by=1))
    
    # Initialize an empty data frame to store ICE5G grid results
    ICE5G_ice_sheet_grid <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset ICE5G data for the current timeslice
      timeslice_subset <- subset(ICE5G_ice_sheet_data, Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      Ice_ICE5G_box <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(Ice_ICE5G_box) <- c("Longitude","Latitude", "Timeslice", "Grid_number", "ICE5G")
      
      # Populate longitude, latitude, and timeslice data
      Ice_ICE5G_box[ ,1:2] <- Coordinate_list_T31
      Ice_ICE5G_box[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude bounds for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset ICE5G data that falls within the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update grid count and ICE5G value based on the subset
          Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If only one matching record, assign directly
            Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$ICE5G
            
          } else {
            
            # Handle multiple matching records
            if(sum(table(subset_site$ICE5G) == max(table(subset_site$ICE5G))) == 1) {
              
              # No tie, assign the most frequent value
              Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$ICE5G))[which.max(table(subset_site$ICE5G))]
              
            } else {
              
              Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 5]     <- 100
              
            }
            
          }
          
        }else{
          
          # No matching records
          Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
          Ice_ICE5G_box[Ice_ICE5G_box$Longitude == Coordinate_list_T31[i, 1] & Ice_ICE5G_box$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
          
        }
        
      }
      
      # Append the result for the current timeslice
      ICE5G_ice_sheet_grid <- rbind(ICE5G_ice_sheet_grid, Ice_ICE5G_box)
      
    }
    
  }  
  
  
  # (2) Create Ice-sheet ensemble
  {
    # Process ICE5G data
    ICE5G_ice_sheet_grid[ ,5]       <- ifelse(ICE5G_ice_sheet_grid[ ,5] == 100, -2, NA)
    
    # Filter and adjust ICE6G and GLAC1D datasets
    ICE6G_GLAC1D_ice_sheet_timeslice           <- ICE6G_GLAC1D_ice_sheet_grid[ICE6G_GLAC1D_ice_sheet_grid$Timeslice %in% seq(-21000, 0, by = 500), ]
    ICE6G_GLAC1D_ice_sheet_timeslice$Timeslice <- ICE6G_GLAC1D_ice_sheet_timeslice$Timeslice /-1000
    ICE6G_GLAC1D_ice_sheet_timeslice[, 4:5][ICE6G_GLAC1D_ice_sheet_timeslice[, 4:5] == 0] <- NA
    
    # Combine ICE5G and ICE6G/GLAC1D datasets
    ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid <- full_join(ICE6G_GLAC1D_ice_sheet_timeslice, ICE5G_ice_sheet_grid[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice"))
    
    # Calculate the ensemble ice-sheet
    ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid$Ice_ensemble <- ifelse(!is.na(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 4]), ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 4], ifelse(!is.na(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 5]), ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 5], ifelse(!is.na(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 6]), ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[, 6], NA)))  
    
    # Save the final ensemble results
    write.csv(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid, file="result/Statistical analysis 3_result.1-Ice sheet ensemble for grid (T31)_0-21ka.csv", row.names=FALSE)
    
  }
  
}


# 3.2 Exclude the site covered by ice-sheet ensemble
# -------------------------------------------------------------------------------------------------
{
  ICE5G_ICE6G_GLAC1D_ice_sheet_ensemble_timeslice_grid  <- na.omit(ICE5G_ICE6G_GLAC1D_ice_sheet_timeslice_grid[ ,c(1:3,7)])
  
  {
    timeslice_list <- seq(0,21,by=0.5)
    
    resulttable_complete <- NULL
    
    for(t in 1:length(timeslice_list)){
      
      Ice_sheet_subset <- subset(ICE5G_ICE6G_GLAC1D_ice_sheet_ensemble_timeslice_grid, Timeslice == timeslice_list[t])
      
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        if(nrow(subset(Ice_sheet_subset, Longitude == Coordinate_list_T31[i, 1] & Latitude == Coordinate_list_T31[i, 2])) > 0){
          
          Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
          Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
          
          Biome_subset     <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka[ ,1:4], Timeslice == timeslice_list[t] & Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
          
        }else{
          
          Biome_subset <- NULL
          
        }
        
        resulttable_complete <- rbind(resulttable_complete, Biome_subset)
      }
      
    }
    
  }
  
  # exclude the grid-cell covered by ice
  Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded <- anti_join(Reconstructed_simulated_biome_timeslice_site_0_21ka, resulttable_complete, by = names(resulttable_complete))
  
  # Rename columns 5 to 11 in the data frame 
  names(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded)[5:11] <- c("Pollen-based reconstruction", "MPI-ESM_GLAC1D", "MPI-ESM_ICE6G", "CLIMBER-X_GLAC1D", "CLIMBER-X_ICE6G", "TRACE-21K-I_ICE5G", "TRACE-21K-II_ICE5G")
  
  # Save the final ensemble results
  write.csv(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded, file="result/Statistical analysis 3_result.2-Reconstructed and simulated biome at timeslice for site since 21 ka BP_ice excluded.csv", row.names=FALSE)
  
}


# 3.3 Assign reconstructed biome at timeslice of site to grid-cell (Asia as example)
# -------------------------------------------------------------------------------------------------

# Load dataset
Taxa_PFT_Asia            <- read.csv2("data/Statistical analysis 3_data.5-Taxa-PFTs_Asia.csv", sep=";", stringsAsFactors=FALSE, check.names=FALSE)  
PFT_Biome_Asia           <- read.csv2("data/Statistical analysis 3_data.6-PFTs-MegaBiomes_Asia.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep=";")    
LegacyPollen_overview    <- read.csv2("data/Statistical analysis 3_data.7-Overview of LegacyPollen2.0 dataset.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep=",")    

# Data Type Conversion
Taxa_PFT_Asia            <- type.convert(Taxa_PFT_Asia, as.is = TRUE)
PFT_Biome_Asia           <- type.convert(PFT_Biome_Asia, as.is = TRUE)
LegacyPollen_overview    <- type.convert(LegacyPollen_overview, as.is = TRUE)

# Rename the First Column
names(Taxa_PFT_Asia)[1] <- "Taxa"
names(PFT_Biome_Asia)[1] <- "Biome"

# Calculate the number of PFT and taxa
{
  # Create an Empty DataFrame for Storing Results
  Taxa_PFT_number_biome <- data.frame(matrix(NA, nrow=8, ncol=5, byrow=TRUE))
  colnames(Taxa_PFT_number_biome) <- c("Biome_code", "Biome", "PFT", "Taxa", "Priority")
  
  # Assign Biome Names
  Taxa_PFT_number_biome$Biome_code <- 1:8
  Taxa_PFT_number_biome$Biome <- PFT_Biome_Asia$Biome
  
  # Generate Unique Biome List
  Biome_list <- unique(PFT_Biome_Asia$Biome)
  
  # Calculate PFT and Taxa Counts for Each Biome
  for (i in Biome_list) {
    
    Taxa_PFT_number_biome[Taxa_PFT_number_biome$Biome == i, 3] <- rowSums(PFT_Biome_Asia[PFT_Biome_Asia$Biome == i, 2:ncol(PFT_Biome_Asia)] == 1)
    
    Taxa_PFT_number_biome[Taxa_PFT_number_biome$Biome == i, 4] <- sum(apply(Taxa_PFT_Asia[, c("Taxa", intersect(colnames(PFT_Biome_Asia)[PFT_Biome_Asia[PFT_Biome_Asia$Biome == i, -1] == 1], colnames(Taxa_PFT_Asia)))][, -1], 1, function(x) any(x == 1)))
    
  }
  
  # Rank Biomes by PFT and Taxa Counts
  Taxa_PFT_number_biome <- Taxa_PFT_number_biome %>%
    arrange(PFT, Taxa) %>%                     
    mutate(Priority = rank(PFT, ties.method = "first")) 
  
}


# Calculation
{
  # Define the timeslices for analysis
  timeslice_list <- seq(0,21,by=0.5)
  
  # Initialize an empty data frame to store reconstructed biome data across all timeslices
  Reconstructed_biome_timeslice_grid_0_21ka <- NULL
  
  # Loop through each timeslice
  for(t in 1:length(timeslice_list)){
    
    # Subset the data for the current timeslice (columns 1 to 5 of the data frame)
    timeslice_pollen_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,1:5], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
    
    # Initialize a result matrix for the current timeslice
    biome_pollen_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
    colnames(biome_pollen_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "Pollen-based reconstruction")
    
    # Populate longitude, latitude, and timeslice information
    biome_pollen_grid[ ,1:2] <- Coordinate_list_T31
    biome_pollen_grid[ ,3]   <- timeslice_list[t]
    
    # Loop through each grid cell
    for (i in 1:nrow(Coordinate_list_T31)) {
      
      print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
      
      # Define the boundaries of the current grid cell
      Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
      Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
      
      # Subset the pollen data for the current grid cell within the timeslice
      site_timeslice_pollen_subset <- na.omit(subset(timeslice_pollen_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
      
      if(nrow(site_timeslice_pollen_subset)>0){
        
        # Update the number of sites in the current grid cell
        biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(site_timeslice_pollen_subset)
        
        if(nrow(site_timeslice_pollen_subset) == 1) {
          
          # If there is only one matching site, assign its reconstruction value directly
          biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- site_timeslice_pollen_subset$`Pollen-based reconstruction`
          
        } else {
          
          # If there are multiple matching sites, handle the reconstruction values
          if(sum(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`) == max(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`))) == 1) {
            
            # No tie, assign the most frequent reconstruction value
            biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`))[which.max(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`))]
            
          } else {
            
            # Identify the most frequent biomes in the pollen-based reconstruction
            most_frequent_biomes <- names(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`))[table(site_timeslice_pollen_subset$`Pollen-based reconstruction`) == max(table(site_timeslice_pollen_subset$`Pollen-based reconstruction`))]
            
            # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
            biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 5] <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
            
          }
          
        }
        
      }else{
        
        # If no matching sites, set values to indicate no data
        biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
        biome_pollen_grid[biome_pollen_grid$Longitude == Coordinate_list_T31[i, 1] & biome_pollen_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
        
      }
      
    }
    
    # Append the result for the current timeslice to the main data frame
    Reconstructed_biome_timeslice_grid_0_21ka <- rbind(Reconstructed_biome_timeslice_grid_0_21ka, biome_pollen_grid)
    
  }
  
}


# 3.4 Creat simulated biome ensemble for grid-cell (Asia as example)
# -------------------------------------------------------------------------------------------------
{
  # Define the list of timeslices from 0 to 21 ka with 0.5 ka increments
  timeslice_list <- seq(0,21,by=0.5)
  
  # Initialize an empty data frame to store the simulation ensemble results for all timeslices
  Simulated_ensemble_biome_timeslice_grid_0_21ka <- NULL
  
  # Loop through each timeslice
  for(t in 1:length(timeslice_list)){
    
    # Subset the data for the current timeslice (columns 1 to 4 and 6 to 11)
    timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,6:11)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
    
    # Initialize a result matrix for the current timeslice
    Simulated_ensemble_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
    colnames(Simulated_ensemble_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "Simulation_Ensemble")
    
    # Populate longitude, latitude, and timeslice information
    Simulated_ensemble_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
    Simulated_ensemble_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
    
    # Loop through each grid cell
    for (i in 1:nrow(Coordinate_list_T31)) {
      
      print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
      
      # Define the longitude and latitude bounds for the current grid cell
      Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
      Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
      
      # Subset the data for the current grid cell
      subset_site <- subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
      
      if(nrow(subset_site)>0){
        
        # Update the number of sites in the current grid cell
        Simulated_ensemble_biome_timeslice_grid[Simulated_ensemble_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & Simulated_ensemble_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
        
        # Extract and count the simulation values from the subset
        {
          values <- list()
          
          for (j in 1:nrow(subset_site)) {
            
            row_values <- subset_site[j, 5:ncol(subset_site)] # Extract simulation values for the row
            
            row_values <- row_values[!is.na(row_values)] # Exclude NA values
            
            values <- c(values, row_values)
          }
          
          frequency <- table(unlist(values)) # Count the frequency of each simulation value
          
        }
        
        # Determine the most frequent simulation value
        if(sum(frequency == max(frequency)) == 1) {
          
          # If there is a clear most frequent value
          Simulated_ensemble_biome_timeslice_grid[Simulated_ensemble_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & Simulated_ensemble_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(frequency)[which.max(frequency)]
          
        } else {
        
          # If there are multiple most frequent values (ties), resolve using Priority ranking
          tied_biomes <- names(frequency)[frequency == max(frequency)]
          
          # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
          Simulated_ensemble_biome_timeslice_grid[Simulated_ensemble_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & Simulated_ensemble_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]  <- tied_biomes[which.min(Taxa_PFT_number_biome$Priority[match(tied_biomes, Taxa_PFT_number_biome$Biome_code)])]
        }
        
        
      }else{
        
        # If no data for the grid cell, set default values
        Simulated_ensemble_biome_timeslice_grid[Simulated_ensemble_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & Simulated_ensemble_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
        
        Simulated_ensemble_biome_timeslice_grid[Simulated_ensemble_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & Simulated_ensemble_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
        
      }
      
    }
    
    # Append the results for the current timeslice to the main data frame
    Simulated_ensemble_biome_timeslice_grid_0_21ka <- rbind(Simulated_ensemble_biome_timeslice_grid_0_21ka, Simulated_ensemble_biome_timeslice_grid)
    
  }
  
}


# 3.5 Assign simulated biome at timeslice of site to grid-cell (Asia as example)
# -------------------------------------------------------------------------------------------------
{
  # (1) MPI-ESM_GLAC1D
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    MPI_GLAC1D_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,6)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      MPI_GLAC1D_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(MPI_GLAC1D_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "MPI-ESM_GLAC1D")
      
      # Populate longitude, latitude, and timeslice information
      MPI_GLAC1D_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      MPI_GLAC1D_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`MPI-ESM_GLAC1D`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`MPI-ESM_GLAC1D`) == max(table(subset_site$`MPI-ESM_GLAC1D`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`MPI-ESM_GLAC1D`))[which.max(table(subset_site$`MPI-ESM_GLAC1D`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`MPI-ESM_GLAC1D`))[table(subset_site$`MPI-ESM_GLAC1D`) == max(table(subset_site$`MPI-ESM_GLAC1D`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]  <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
              
              
            }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]  <- 0
          
          MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      MPI_GLAC1D_biome_timeslice_grid_0_21ka <- rbind(MPI_GLAC1D_biome_timeslice_grid_0_21ka, MPI_GLAC1D_biome_timeslice_grid)
      
    }
    
  }
  
  # (2) MPI-ESM_ICE6G
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    MPI_ICE6G_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,7)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      MPI_ICE6G_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(MPI_ICE6G_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "MPI-ESM_ICE6G")
      
      # Populate longitude, latitude, and timeslice information
      MPI_ICE6G_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      MPI_ICE6G_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          MPI_ICE6G_biome_timeslice_grid[MPI_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            MPI_ICE6G_biome_timeslice_grid[MPI_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`MPI-ESM_ICE6G`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`MPI-ESM_ICE6G`) == max(table(subset_site$`MPI-ESM_ICE6G`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              MPI_ICE6G_biome_timeslice_grid[MPI_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`MPI-ESM_ICE6G`))[which.max(table(subset_site$`MPI-ESM_ICE6G`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`MPI-ESM_ICE6G`))[table(subset_site$`MPI-ESM_ICE6G`) == max(table(subset_site$`MPI-ESM_ICE6G`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              MPI_GLAC1D_biome_timeslice_grid[MPI_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]  <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
            
            }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          MPI_ICE6G_biome_timeslice_grid[MPI_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- 0
          
          MPI_ICE6G_biome_timeslice_grid[MPI_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & MPI_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      MPI_ICE6G_biome_timeslice_grid_0_21ka <- rbind(MPI_ICE6G_biome_timeslice_grid_0_21ka, MPI_ICE6G_biome_timeslice_grid)
      
    }
    
  }
  
  # (3) CLIMBER-X_GLAC1D
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    CLIMBER_X_GLAC1D_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,8)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      CLIMBER_X_GLAC1D_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(CLIMBER_X_GLAC1D_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "CLIMBER-X_GLAC1D")
      
      # Populate longitude, latitude, and timeslice information
      CLIMBER_X_GLAC1D_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      CLIMBER_X_GLAC1D_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`CLIMBER-X_GLAC1D`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`CLIMBER-X_GLAC1D`) == max(table(subset_site$`CLIMBER-X_GLAC1D`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`CLIMBER-X_GLAC1D`))[which.max(table(subset_site$`CLIMBER-X_GLAC1D`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`CLIMBER-X_GLAC1D`))[table(subset_site$`CLIMBER-X_GLAC1D`) == max(table(subset_site$`CLIMBER-X_GLAC1D`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
            }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
          
          CLIMBER_X_GLAC1D_biome_timeslice_grid[CLIMBER_X_GLAC1D_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_GLAC1D_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      CLIMBER_X_GLAC1D_biome_timeslice_grid_0_21ka <- rbind(CLIMBER_X_GLAC1D_biome_timeslice_grid_0_21ka, CLIMBER_X_GLAC1D_biome_timeslice_grid)
      
    }
    
  }
  
  
  # (4) CLIMBER-X_ICE6G
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    CLIMBER_X_ICE6G_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,9)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      CLIMBER_X_ICE6G_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(CLIMBER_X_ICE6G_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "CLIMBER-X_ICE6G")
      
      # Populate longitude, latitude, and timeslice information
      CLIMBER_X_ICE6G_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      CLIMBER_X_ICE6G_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`CLIMBER-X_ICE6G`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`CLIMBER-X_ICE6G`) == max(table(subset_site$`CLIMBER-X_ICE6G`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`CLIMBER-X_ICE6G`))[which.max(table(subset_site$`CLIMBER-X_ICE6G`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`CLIMBER-X_ICE6G`))[table(subset_site$`CLIMBER-X_ICE6G`) == max(table(subset_site$`CLIMBER-X_ICE6G`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]    <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
              
            }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- 0
          
          CLIMBER_X_ICE6G_biome_timeslice_grid[CLIMBER_X_ICE6G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & CLIMBER_X_ICE6G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      CLIMBER_X_ICE6G_biome_timeslice_grid_0_21ka <- rbind(CLIMBER_X_ICE6G_biome_timeslice_grid_0_21ka, CLIMBER_X_ICE6G_biome_timeslice_grid)
      
    }
    
  }
  
  
  # (5) TRACE-21K-I_ICE5G
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    TRACE_21K_I_ICE5G_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,10)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      TRACE_21K_I_ICE5G_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(TRACE_21K_I_ICE5G_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "TRACE-21K-I_ICE5G")
      
      # Populate longitude, latitude, and timeslice information
      TRACE_21K_I_ICE5G_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      TRACE_21K_I_ICE5G_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`TRACE-21K-I_ICE5G`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`TRACE-21K-I_ICE5G`) == max(table(subset_site$`TRACE-21K-I_ICE5G`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`TRACE-21K-I_ICE5G`))[which.max(table(subset_site$`TRACE-21K-I_ICE5G`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`TRACE-21K-I_ICE5G`))[table(subset_site$`TRACE-21K-I_ICE5G`) == max(table(subset_site$`TRACE-21K-I_ICE5G`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]    <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
              
            }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
          
          TRACE_21K_I_ICE5G_biome_timeslice_grid[TRACE_21K_I_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_I_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      TRACE_21K_I_ICE5G_biome_timeslice_grid_0_21ka <- rbind(TRACE_21K_I_ICE5G_biome_timeslice_grid_0_21ka, TRACE_21K_I_ICE5G_biome_timeslice_grid)
      
    }
    
  }
  
  
  # (6) TRACE-21K-II_ICE5G
  {
    # Define the list of timeslices from 0 to 21 ka with increments of 0.5 ka
    timeslice_list <- seq(0,21,by=0.5)
    
    # Initialize an empty data frame to store results for all timeslices
    TRACE_21K_II_ICE5G_biome_timeslice_grid_0_21ka <- NULL
    
    # Loop through each timeslice
    for(t in 1:length(timeslice_list)){
      
      # Subset the data for the current timeslice
      timeslice_subset <- subset(Reconstructed_simulated_biome_timeslice_site_0_21ka_excluded[ ,c(1:4,11)], Dataset_ID %in% LegacyPollen_overview$Dataset_ID[LegacyPollen_overview$Continent == "Asia"] & Timeslice == timeslice_list[t])
      
      # Initialize a result matrix for the current timeslice
      TRACE_21K_II_ICE5G_biome_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(Coordinate_list_T31), ncol=5, byrow=TRUE))
      colnames(TRACE_21K_II_ICE5G_biome_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "Site_number", "TRACE-21K-II_ICE5G")
      
      # Populate longitude, latitude, and timeslice information
      TRACE_21K_II_ICE5G_biome_timeslice_grid[ ,1:2] <- Coordinate_list_T31
      TRACE_21K_II_ICE5G_biome_timeslice_grid[ ,3]   <- timeslice_list[t]
      
      # Loop through each grid cell
      for (i in 1:nrow(Coordinate_list_T31)) {
        
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(Coordinate_list_T31), " +++++"))
        
        # Define the longitude and latitude boundaries for the current grid cell
        Longitude_subset <- c(Coordinate_list_T31[i, 1] - 360/Longitude_number/2, Coordinate_list_T31[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(Coordinate_list_T31[i, 2] - 180/Latitude_number/2,  Coordinate_list_T31[i, 2] + 180/Latitude_number/2)
        
        # Subset the data for the current grid cell
        subset_site <- na.omit(subset(timeslice_subset, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2]))
        
        if(nrow(subset_site)>0){
          
          # Update the number of sites in the current grid cell
          TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]     <- nrow(subset_site)
          
          if(nrow(subset_site) == 1) {
            
            # If there is only one site, assign its reconstruction value directly
            TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- subset_site$`TRACE-21K-II_ICE5G`
            
          } else {
            
            # Handle multiple sites in the grid cell
            if(sum(table(subset_site$`TRACE-21K-II_ICE5G`) == max(table(subset_site$`TRACE-21K-II_ICE5G`))) == 1) {
              
              # If no tie, assign the most frequent reconstruction value
              TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]     <- names(table(subset_site$`TRACE-21K-II_ICE5G`))[which.max(table(subset_site$`TRACE-21K-II_ICE5G`))]
              
            } else {
              
              # Identify the most frequent biomes 
              most_frequent_biomes <- names(table(subset_site$`TRACE-21K-II_ICE5G`))[table(subset_site$`TRACE-21K-II_ICE5G`) == max(table(subset_site$`TRACE-21K-II_ICE5G`))]
              
              # Assign the most frequent biome to the grid cell, resolving ties using Priority ranking
              TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]  <- most_frequent_biomes[which.min(Taxa_PFT_number_biome$Priority[match(most_frequent_biomes, Taxa_PFT_number_biome$Biome_code)])]
           }
            
          }
          
        }else{
          
          # If no data, set default values for the grid cell
          TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 4]   <- 0
          
          TRACE_21K_II_ICE5G_biome_timeslice_grid[TRACE_21K_II_ICE5G_biome_timeslice_grid$Longitude == Coordinate_list_T31[i, 1] & TRACE_21K_II_ICE5G_biome_timeslice_grid$Latitude == Coordinate_list_T31[i, 2], 5]   <- NA
          
        }
        
      }
      
      # Append results for the current timeslice to the main result container
      TRACE_21K_II_ICE5G_biome_timeslice_grid_0_21ka <- rbind(TRACE_21K_II_ICE5G_biome_timeslice_grid_0_21ka, TRACE_21K_II_ICE5G_biome_timeslice_grid)
      
    }
    
  }
  
  
}


# -------------------------------------------------------------------------------------------------
# Combine datasets

Reconstructed_simulated_biome_timeslice_grid_0_21ka_final  <- full_join(Reconstructed_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], MPI_GLAC1D_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(MPI_ICE6G_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(CLIMBER_X_GLAC1D_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(CLIMBER_X_ICE6G_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(TRACE_21K_I_ICE5G_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(TRACE_21K_II_ICE5G_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) %>%
  full_join(Simulated_ensemble_biome_timeslice_grid_0_21ka[ ,c(1:3,5)], by = c("Longitude", "Latitude", "Timeslice")) 

# Save the final ensemble results
write.csv(Reconstructed_simulated_biome_timeslice_grid_0_21ka_final, file="result/Statistical analysis 3_result.3-Reconstructed and simulated biome at timeslice for grid-cell since 21 ka BP.csv", row.names=FALSE)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 4. Statistical analysis 4: Validation of model climate bias ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# 4.1 Calculate bioclimatic variables of observation
# -------------------------------------------------------------------------------------------------

# (1) 1931-1970 C.E
{
  # Load datasets
  CRU_climate_1931_1970CE_grid_global   <- read.csv2("data/Statistical analysis 4_data.1-Observed climate of 1931-1970 from CRU_grid cell (0.5 degrees).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
  
  # Data Type Conversion
  CRU_climate_1931_1970CE_grid_global   <- type.convert(CRU_climate_1931_1970CE_grid_global, as.is = TRUE) 
  
  # Calculation
  {
    # Extract coordinate and year lists
    coordinate_list <- as.data.frame(unique(CRU_climate_1931_1970CE_grid_global[ ,1:2])) # The unique combination of longitude and latitude defines grid cells
    year_list       <- sort(unique(CRU_climate_1931_1970CE_grid_global$Year)) # Get a sorted list of unique years from the dataset
    
    # Calculate bioclimatic variables
    {
      # Set the size of each processing chunk
      chunk_size <- 100  # Process 100 grid cells at a time
      
      # Define the output file to save results
      output_file <- "result/Statistical analysis 4_result.1-Observed bioclimatic variables of 1931-1970 from CRU_grid cell (0.5 degrees).csv"
      
      # Loop through each processing chunk
      for (k in seq(1, nrow(coordinate_list), by = chunk_size)) {
        
        # Determine the end index of the current chunk
        chunk_end <- min(k + chunk_size - 1, nrow(coordinate_list))
        
        {
          # Extract a subset (batch) of the coordinate list for the current chunk
          coordinate_list_batch <- coordinate_list[k:chunk_end, ]
          
          # Initialize an empty data frame to store bioclimatic variables for the batch
          Bioclimatic_variables_1931_1970_batch <- NULL
          
          # Join the CRU climate data with the current batch of coordinates
          CRU_climate_1931_1970CE_grid_global_batch <- right_join(CRU_climate_1931_1970CE_grid_global, coordinate_list_batch, by = join_by(Long, Lat))
          
          # Loop through each grid cell in the batch
          for (j in 1:nrow(coordinate_list_batch)) {
            
            # Print progress for the current grid cell
            print(paste0("+++++ Batch-", chunk_end/chunk_size, "/", ceiling(nrow(coordinate_list)/chunk_size), " :grid-cell ", j, "/", nrow(coordinate_list_batch)," +++++"))
            
            # Extract data for the current grid cell
            subset_grid <- CRU_climate_1931_1970CE_grid_global_batch[CRU_climate_1931_1970CE_grid_global_batch$Long == coordinate_list_batch[j, 1] & CRU_climate_1931_1970CE_grid_global_batch$Lat == coordinate_list_batch[j, 2], ]
            
            # Initialize the result matrix for the grid cell
            Bioclimatic_variables_grid           <- data.frame(matrix(NA, nrow= length(year_list), ncol= 8, byrow=TRUE))
            colnames(Bioclimatic_variables_grid) <- c("Source", "Longitude", "Latitude", "Year", "BIO10", "BIO11", "BIO18", "BIO19")
            
            # Assign metadata for the grid cell
            Bioclimatic_variables_grid[ ,1] <- "CRU TS v. 4.08"
            Bioclimatic_variables_grid[ ,2] <- unique(subset_grid$Long)
            Bioclimatic_variables_grid[ ,3] <- unique(subset_grid$Lat)
            
            # Loop through each year to calculate bioclimatic variables
            for (i in 1:length(year_list)) {
              
              # Extract data for the specific year
              subset_grid_subset <- na.omit(subset(subset_grid, Year == year_list[i]))
              
              # Assign the year
              Bioclimatic_variables_grid[i, 4] <- year_list[i]
              
              
              if(nrow(subset_grid_subset) != 12){
                
                # If the data for the year is incomplete (not all 12 months present), assign NA
                Bioclimatic_variables_grid[i, 5:8] <- NA
                
              }else{
                
                # Calculate quarterly means for precipitation and temperature
                # -------------------------------------------------------------------------------------------------
                ## Precipitation
                subset_grid_subset$precipitation_1_2_3    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 1] + subset_grid_subset$PRE[subset_grid_subset$Month == 2] + subset_grid_subset$PRE[subset_grid_subset$Month == 3])
                subset_grid_subset$precipitation_2_3_4    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 2] + subset_grid_subset$PRE[subset_grid_subset$Month == 3] + subset_grid_subset$PRE[subset_grid_subset$Month == 4])
                subset_grid_subset$precipitation_3_4_5    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 3] + subset_grid_subset$PRE[subset_grid_subset$Month == 4] + subset_grid_subset$PRE[subset_grid_subset$Month == 5])
                subset_grid_subset$precipitation_4_5_6    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 4] + subset_grid_subset$PRE[subset_grid_subset$Month == 5] + subset_grid_subset$PRE[subset_grid_subset$Month == 6])
                subset_grid_subset$precipitation_5_6_7    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 5] + subset_grid_subset$PRE[subset_grid_subset$Month == 6] + subset_grid_subset$PRE[subset_grid_subset$Month == 7])
                subset_grid_subset$precipitation_6_7_8    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 6] + subset_grid_subset$PRE[subset_grid_subset$Month == 7] + subset_grid_subset$PRE[subset_grid_subset$Month == 8])
                subset_grid_subset$precipitation_7_8_9    <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 7] + subset_grid_subset$PRE[subset_grid_subset$Month == 8] + subset_grid_subset$PRE[subset_grid_subset$Month == 9])
                subset_grid_subset$precipitation_8_9_10   <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 8] + subset_grid_subset$PRE[subset_grid_subset$Month == 9] + subset_grid_subset$PRE[subset_grid_subset$Month == 10])
                subset_grid_subset$precipitation_9_10_11  <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 9] + subset_grid_subset$PRE[subset_grid_subset$Month == 10] + subset_grid_subset$PRE[subset_grid_subset$Month == 11])
                subset_grid_subset$precipitation_10_11_12 <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 10] + subset_grid_subset$PRE[subset_grid_subset$Month == 11] + subset_grid_subset$PRE[subset_grid_subset$Month == 12])
                subset_grid_subset$precipitation_11_12_1  <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 11] + subset_grid_subset$PRE[subset_grid_subset$Month == 12] + subset_grid_subset$PRE[subset_grid_subset$Month == 1])
                subset_grid_subset$precipitation_12_1_2   <- mean(subset_grid_subset$PRE[subset_grid_subset$Month == 12] + subset_grid_subset$PRE[subset_grid_subset$Month == 1] + subset_grid_subset$PRE[subset_grid_subset$Month == 2])
                
                ## Temperature
                subset_grid_subset$temperature_1_2_3    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 1] + subset_grid_subset$TMP[subset_grid_subset$Month == 2] + subset_grid_subset$TMP[subset_grid_subset$Month == 3])
                subset_grid_subset$temperature_2_3_4    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 2] + subset_grid_subset$TMP[subset_grid_subset$Month == 3] + subset_grid_subset$TMP[subset_grid_subset$Month == 4])
                subset_grid_subset$temperature_3_4_5    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 3] + subset_grid_subset$TMP[subset_grid_subset$Month == 4] + subset_grid_subset$TMP[subset_grid_subset$Month == 5])
                subset_grid_subset$temperature_4_5_6    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 4] + subset_grid_subset$TMP[subset_grid_subset$Month == 5] + subset_grid_subset$TMP[subset_grid_subset$Month == 6])
                subset_grid_subset$temperature_5_6_7    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 5] + subset_grid_subset$TMP[subset_grid_subset$Month == 6] + subset_grid_subset$TMP[subset_grid_subset$Month == 7])
                subset_grid_subset$temperature_6_7_8    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 6] + subset_grid_subset$TMP[subset_grid_subset$Month == 7] + subset_grid_subset$TMP[subset_grid_subset$Month == 8])
                subset_grid_subset$temperature_7_8_9    <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 7] + subset_grid_subset$TMP[subset_grid_subset$Month == 8] + subset_grid_subset$TMP[subset_grid_subset$Month == 9])
                subset_grid_subset$temperature_8_9_10   <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 8] + subset_grid_subset$TMP[subset_grid_subset$Month == 9] + subset_grid_subset$TMP[subset_grid_subset$Month == 10])
                subset_grid_subset$temperature_9_10_11  <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 9] + subset_grid_subset$TMP[subset_grid_subset$Month == 10] + subset_grid_subset$TMP[subset_grid_subset$Month == 11])
                subset_grid_subset$temperature_10_11_12 <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 10] + subset_grid_subset$TMP[subset_grid_subset$Month == 11] + subset_grid_subset$TMP[subset_grid_subset$Month == 12])
                subset_grid_subset$temperature_11_12_1  <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 11] + subset_grid_subset$TMP[subset_grid_subset$Month == 12] + subset_grid_subset$TMP[subset_grid_subset$Month == 1])
                subset_grid_subset$temperature_12_1_2   <- mean(subset_grid_subset$TMP[subset_grid_subset$Month == 12] + subset_grid_subset$TMP[subset_grid_subset$Month == 1] + subset_grid_subset$TMP[subset_grid_subset$Month == 2])
                
                # Identify the warmest and coldest quarters
                subset_grid_subset$Warmest_Quarter   <- apply(subset_grid_subset[ ,which(names(subset_grid_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
                subset_grid_subset$Coldest_Quarter   <- apply(subset_grid_subset[ ,which(names(subset_grid_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
                
                # Calculate bioclimatic variables
                # -------------------------------------------------------------------------------------------------
                
                ## BIO10 = Mean Temperature of Warmest Quarter
                Bioclimatic_variables_grid[i, 5] <- (subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid_subset$Warmest_Quarter)))] +
                                                       subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid_subset$Warmest_Quarter)))] +  
                                                       subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid_subset$Warmest_Quarter)))])/3
                
                ## BIO11 = Mean Temperature of Coldest Quarter
                Bioclimatic_variables_grid[i, 6] <- (subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid_subset$Coldest_Quarter)))] +
                                                       subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid_subset$Coldest_Quarter)))] +  
                                                       subset_grid_subset$TMP[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid_subset$Coldest_Quarter)))])/3
                
                ## BIO18 = Precipitation of Warmest Quarter
                Bioclimatic_variables_grid[i, 7] <- (subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid_subset$Warmest_Quarter)))] +
                                                       subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid_subset$Warmest_Quarter)))] +  
                                                       subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid_subset$Warmest_Quarter)))])
                
                ## BIO19 = Precipitation of Coldest Quarter
                Bioclimatic_variables_grid[i, 8] <- (subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid_subset$Coldest_Quarter)))] +
                                                       subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid_subset$Coldest_Quarter)))] +  
                                                       subset_grid_subset$PRE[subset_grid_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid_subset$Coldest_Quarter)))])
                
                
                
                
              }
              
            }
            
            # Append the bioclimatic variables for the current grid cell to the batch result
            Bioclimatic_variables_1931_1970_batch <- rbind(Bioclimatic_variables_1931_1970_batch, Bioclimatic_variables_grid)
            
          }
          
        } 
        
        # Append the current chunk to the CSV file
        fwrite(Bioclimatic_variables_1931_1970_batch, file = output_file, append = TRUE, col.names = !file.exists(output_file))
        
      }
      
    }
    
    # Remove the object from memory
    rm(CRU_climate_1931_1970CE_grid_global) 
  }
  
}


# (2) Calculate climatological mean of 1931-1970 C.E 
{
  # Load datasets
  CRU_Bioclimatic_variables_1931_1970CE_grid_global   <- read.csv2("result/Statistical analysis 4_result.1-Observed bioclimatic variables of 1931-1970 from CRU_grid cell (0.5 degrees).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
  
  # Data Type Conversion
  CRU_Bioclimatic_variables_1931_1970CE_grid_global   <- type.convert(CRU_Bioclimatic_variables_1931_1970CE_grid_global, as.is = TRUE) 
  
  # Extract unique longitude and latitude pairs
  coordinate_list <- unique(CRU_Bioclimatic_variables_1931_1970CE_grid_global[, 2:3])
  
  # Initialize a data frame to store the mean bioclimatic variables for each grid cell
  Bioclimatic_variables_mean_1931_1970 <- data.frame(
    Source = "CRU TS v. 4.08 Bioclimatic variables_mean of 1931-1970", # Add a source column with a descriptive label
    Longitude = coordinate_list$Longitude, # Longitude values from the coordinate list
    Latitude = coordinate_list$Latitude, # Latitude values from the coordinate list
    matrix(NA, nrow = nrow(coordinate_list), ncol = 4) # Initialize columns for bioclimatic variables
  )
  
  # Set the column names for the bioclimatic variables
  colnames(Bioclimatic_variables_mean_1931_1970)[4:7] <- names(CRU_Bioclimatic_variables_1931_1970CE_grid_global)[5:8]
  
  # Compute the mean bioclimatic variables for each grid cell
  Bioclimatic_variables_mean_1931_1970[, 4:7] <- CRU_Bioclimatic_variables_1931_1970CE_grid_global %>%
    group_by(Longitude, Latitude) %>%  # Group data by longitude and latitude
    summarise(across(BIO10:BIO19, mean, .names = "mean_{col}"), # Calculate the mean for each bioclimatic variable
              .groups = "drop") %>%  # Drop grouping structure after summarization
    select(-Longitude, -Latitude)  # Remove longitude and latitude from the result
  
  # Write the computed mean bioclimatic variables to a CSV file
  write.csv(Bioclimatic_variables_mean_1931_1970, file="result/Statistical analysis 4_result.2-Observed bioclimatic variables of 1931-1970 mean from CRU_grid cell (0.5 degrees).csv", row.names=FALSE)
  
}


# 4.2 Calculate bioclimatic variables of simulation at 0 cal. ka BP
# -------------------------------------------------------------------------------------------------

# Load datasets
Simulated_climate_0ka_grid_global   <- read.csv2("data/Statistical analysis 4_data.2-Simulated climate at 0 ka BP for grid-cell.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Simulated_climate_0ka_grid_global   <- type.convert(Simulated_climate_0ka_grid_global, as.is = TRUE) 

# Calculation
{
  # (1) MPI-ESM_GLAC1D
  {
    # Subset the simulated climate data
    MPI_GLAC1D_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "MPI-ESM_GLAC1D")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_MPI_GLAC1D <- unique(MPI_GLAC1D_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_MPI_GLAC1D_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_MPI_GLAC1D)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_MPI_GLAC1D)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(MPI_GLAC1D_climate_0ka_grid_global, Longitude == coordinate_list_MPI_GLAC1D[j, 1] & Latitude == coordinate_list_MPI_GLAC1D[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_MPI_GLAC1D_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_MPI_GLAC1D_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_MPI_GLAC1D_grid[ ,1] <- "MPI-ESM_GLAC1D"
      bioclimate_MPI_GLAC1D_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_MPI_GLAC1D_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_MPI_GLAC1D_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                             subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                             subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_MPI_GLAC1D_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                             subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                             subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_MPI_GLAC1D_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                             subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                             subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_MPI_GLAC1D_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                             subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                             subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_MPI_GLAC1D_grid_0ka_global <- rbind(Bioclimatic_variables_MPI_GLAC1D_grid_0ka_global, bioclimate_MPI_GLAC1D_grid)
      
    }
    
  }
  
  
  # (2) MPI-ESM_ICE6G
  {
    # Subset the simulated climate data
    MPI_ICE6G_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "MPI-ESM_ICE6G")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_MPI_ICE6G <- unique(MPI_ICE6G_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_MPI_ICE6G_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_MPI_ICE6G)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_MPI_ICE6G)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(MPI_ICE6G_climate_0ka_grid_global, Longitude == coordinate_list_MPI_ICE6G[j, 1] & Latitude == coordinate_list_MPI_ICE6G[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_MPI_ICE6G_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_MPI_ICE6G_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_MPI_ICE6G_grid[ ,1] <- "MPI-ESM_ICE6G"
      bioclimate_MPI_ICE6G_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_MPI_ICE6G_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_MPI_ICE6G_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                            subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                            subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_MPI_ICE6G_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                            subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                            subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_MPI_ICE6G_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                            subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                            subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_MPI_ICE6G_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                            subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                            subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_MPI_ICE6G_grid_0ka_global <- rbind(Bioclimatic_variables_MPI_ICE6G_grid_0ka_global, bioclimate_MPI_ICE6G_grid)
      
    }
    
    
  }
  
  
  # (3) CLIMBER-X_GLAC1D
  {
    # Subset the simulated climate data
    CLIMBER_GLAC1D_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "CLIMBER-X_GLAC1D")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_CLIMBER_GLAC1D <- unique(CLIMBER_GLAC1D_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_CLIMBER_GLAC1D_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_CLIMBER_GLAC1D)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_CLIMBER_GLAC1D)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(CLIMBER_GLAC1D_climate_0ka_grid_global, Longitude == coordinate_list_CLIMBER_GLAC1D[j, 1] & Latitude == coordinate_list_CLIMBER_GLAC1D[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_CLIMBER_GLAC1D_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_CLIMBER_GLAC1D_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_CLIMBER_GLAC1D_grid[ ,1] <- "CLIMBER-X_GLAC1D"
      bioclimate_CLIMBER_GLAC1D_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_CLIMBER_GLAC1D_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_CLIMBER_GLAC1D_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_CLIMBER_GLAC1D_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_CLIMBER_GLAC1D_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_CLIMBER_GLAC1D_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_CLIMBER_GLAC1D_grid_0ka_global <- rbind(Bioclimatic_variables_CLIMBER_GLAC1D_grid_0ka_global, bioclimate_CLIMBER_GLAC1D_grid)
      
    }
    
    
  }
  
  
  # (4) CLIMBER-X_ICE6G
  {
    # Subset the simulated climate data
    CLIMBER_ICE6G_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "CLIMBER-X_ICE6G")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_CLIMBER_ICE6G <- unique(CLIMBER_ICE6G_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_CLIMBER_ICE6G_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_CLIMBER_ICE6G)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_CLIMBER_ICE6G)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(CLIMBER_ICE6G_climate_0ka_grid_global, Longitude == coordinate_list_CLIMBER_ICE6G[j, 1] & Latitude == coordinate_list_CLIMBER_ICE6G[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_CLIMBER_ICE6G_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_CLIMBER_ICE6G_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_CLIMBER_ICE6G_grid[ ,1] <- "CLIMBER-X_ICE6G"
      bioclimate_CLIMBER_ICE6G_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_CLIMBER_ICE6G_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_CLIMBER_ICE6G_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_CLIMBER_ICE6G_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_CLIMBER_ICE6G_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_CLIMBER_ICE6G_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_CLIMBER_ICE6G_grid_0ka_global <- rbind(Bioclimatic_variables_CLIMBER_ICE6G_grid_0ka_global, bioclimate_CLIMBER_ICE6G_grid)
      
    }
    
    
  }
  
  
  # (5) TRACE-21K-I_ICE5G
  {
    # Subset the simulated climate data
    TRACE_I_ICE5G_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "TRACE-21K-I_ICE5G")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_TRACE_I_ICE5G <- unique(TRACE_I_ICE5G_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_TRACE_I_ICE5G_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_TRACE_I_ICE5G)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_TRACE_I_ICE5G)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(TRACE_I_ICE5G_climate_0ka_grid_global, Longitude == coordinate_list_TRACE_I_ICE5G[j, 1] & Latitude == coordinate_list_TRACE_I_ICE5G[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_TRACE_I_ICE5G_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_TRACE_I_ICE5G_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_TRACE_I_ICE5G_grid[ ,1] <- "TRACE-21K-I_ICE5G"
      bioclimate_TRACE_I_ICE5G_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_TRACE_I_ICE5G_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_TRACE_I_ICE5G_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_TRACE_I_ICE5G_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_TRACE_I_ICE5G_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_TRACE_I_ICE5G_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_TRACE_I_ICE5G_grid_0ka_global <- rbind(Bioclimatic_variables_TRACE_I_ICE5G_grid_0ka_global, bioclimate_TRACE_I_ICE5G_grid)
      
    }
    
    
  }
  
  
  # (6) TRACE-21K-II_ICE5G
  {
    # Subset the simulated climate data
    TRACE_II_ICE5G_climate_0ka_grid_global <- subset(Simulated_climate_0ka_grid_global, Source == "TRACE-21K-II_ICE5G")
    
    # Extract unique longitude and latitude pairs for grid cells
    coordinate_list_TRACE_II_ICE5G <- unique(TRACE_II_ICE5G_climate_0ka_grid_global[ ,1:2])
    
    # Initialize an empty data frame to store bioclimatic variables for all grid cells
    Bioclimatic_variables_TRACE_II_ICE5G_grid_0ka_global <- NULL
    
    # Loop through each grid cell
    for (j in 1:nrow(coordinate_list_TRACE_II_ICE5G)) {
      
      # Print progress for the current grid cell
      print(paste0("+++++ grid-cell ", j, "/", nrow(coordinate_list_TRACE_II_ICE5G)," +++++"))
      
      # Extract data for the current grid cell
      subset_grid <- subset(TRACE_II_ICE5G_climate_0ka_grid_global, Longitude == coordinate_list_TRACE_II_ICE5G[j, 1] & Latitude == coordinate_list_TRACE_II_ICE5G[j, 2]) 
      
      # Initialize a result matrix for the current grid cell
      bioclimate_TRACE_II_ICE5G_grid           <- data.frame(matrix(NA, nrow= 1, ncol= 7, byrow=TRUE))
      colnames(bioclimate_TRACE_II_ICE5G_grid) <- c("Source", "Longitude", "Latitude", "BIO10", "BIO11", "BIO18", "BIO19")
      
      # Assign metadata for the current grid cell
      bioclimate_TRACE_II_ICE5G_grid[ ,1] <- "TRACE-21K-II_ICE5G"
      bioclimate_TRACE_II_ICE5G_grid[ ,2] <- unique(subset_grid$Longitude)
      bioclimate_TRACE_II_ICE5G_grid[ ,3] <- unique(subset_grid$Latitude)
      
      # Calculate quarterly mean precipitation and temperature
      # -------------------------------------------------------------------------------------------------
      ## Precipitation
      subset_grid$precipitation_1_2_3    <- mean(subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3])
      subset_grid$precipitation_2_3_4    <- mean(subset_grid$Precipitation[subset_grid$Month == 2] + subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4])
      subset_grid$precipitation_3_4_5    <- mean(subset_grid$Precipitation[subset_grid$Month == 3] + subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5])
      subset_grid$precipitation_4_5_6    <- mean(subset_grid$Precipitation[subset_grid$Month == 4] + subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6])
      subset_grid$precipitation_5_6_7    <- mean(subset_grid$Precipitation[subset_grid$Month == 5] + subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7])
      subset_grid$precipitation_6_7_8    <- mean(subset_grid$Precipitation[subset_grid$Month == 6] + subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8])
      subset_grid$precipitation_7_8_9    <- mean(subset_grid$Precipitation[subset_grid$Month == 7] + subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9])
      subset_grid$precipitation_8_9_10   <- mean(subset_grid$Precipitation[subset_grid$Month == 8] + subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10])
      subset_grid$precipitation_9_10_11  <- mean(subset_grid$Precipitation[subset_grid$Month == 9] + subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11])
      subset_grid$precipitation_10_11_12 <- mean(subset_grid$Precipitation[subset_grid$Month == 10] + subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12])
      subset_grid$precipitation_11_12_1  <- mean(subset_grid$Precipitation[subset_grid$Month == 11] + subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1])
      subset_grid$precipitation_12_1_2   <- mean(subset_grid$Precipitation[subset_grid$Month == 12] + subset_grid$Precipitation[subset_grid$Month == 1] + subset_grid$Precipitation[subset_grid$Month == 2])
      
      ## Temperature
      subset_grid$temperature_1_2_3    <- mean(subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3])
      subset_grid$temperature_2_3_4    <- mean(subset_grid$Temperature[subset_grid$Month == 2] + subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4])
      subset_grid$temperature_3_4_5    <- mean(subset_grid$Temperature[subset_grid$Month == 3] + subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5])
      subset_grid$temperature_4_5_6    <- mean(subset_grid$Temperature[subset_grid$Month == 4] + subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6])
      subset_grid$temperature_5_6_7    <- mean(subset_grid$Temperature[subset_grid$Month == 5] + subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7])
      subset_grid$temperature_6_7_8    <- mean(subset_grid$Temperature[subset_grid$Month == 6] + subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8])
      subset_grid$temperature_7_8_9    <- mean(subset_grid$Temperature[subset_grid$Month == 7] + subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9])
      subset_grid$temperature_8_9_10   <- mean(subset_grid$Temperature[subset_grid$Month == 8] + subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10])
      subset_grid$temperature_9_10_11  <- mean(subset_grid$Temperature[subset_grid$Month == 9] + subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11])
      subset_grid$temperature_10_11_12 <- mean(subset_grid$Temperature[subset_grid$Month == 10] + subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12])
      subset_grid$temperature_11_12_1  <- mean(subset_grid$Temperature[subset_grid$Month == 11] + subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1])
      subset_grid$temperature_12_1_2   <- mean(subset_grid$Temperature[subset_grid$Month == 12] + subset_grid$Temperature[subset_grid$Month == 1] + subset_grid$Temperature[subset_grid$Month == 2])
      
      # Identify the warmest and coldest quarters
      subset_grid$Warmest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
      subset_grid$Coldest_Quarter   <- apply(subset_grid[ ,which(names(subset_grid) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
      
      # Calculate BIO variables
      # -------------------------------------------------------------------------------------------------
      
      ## BIO10 = Mean Temperature of Warmest Quarter
      bioclimate_TRACE_II_ICE5G_grid[1, 4] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])/3
      
      ## BIO11 = Mean Temperature of Coldest Quarter
      bioclimate_TRACE_II_ICE5G_grid[1, 5] <- (subset_grid$Temperature[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                 subset_grid$Temperature[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])/3
      
      
      
      ## BIO18 = Precipitation of Warmest Quarter
      bioclimate_TRACE_II_ICE5G_grid[1, 6] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Warmest_Quarter)))] +
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Warmest_Quarter)))] +  
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Warmest_Quarter)))])
      
      ## BIO19 = Precipitation of Coldest Quarter
      bioclimate_TRACE_II_ICE5G_grid[1, 7] <- (subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_grid$Coldest_Quarter)))] +
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_grid$Coldest_Quarter)))] +  
                                                 subset_grid$Precipitation[subset_grid$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_grid$Coldest_Quarter)))])
      
      # Append the bioclimatic variables of the current grid cell to the global results
      Bioclimatic_variables_TRACE_II_ICE5G_grid_0ka_global <- rbind(Bioclimatic_variables_TRACE_II_ICE5G_grid_0ka_global, bioclimate_TRACE_II_ICE5G_grid)
      
    }
    
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  # Combine datasets
  Bioclimatic_variables_all_model_grid_0ka_global <- rbind(Bioclimatic_variables_MPI_GLAC1D_grid_0ka_global, Bioclimatic_variables_MPI_ICE6G_grid_0ka_global, Bioclimatic_variables_CLIMBER_GLAC1D_grid_0ka_global, 
                                                           Bioclimatic_variables_CLIMBER_ICE6G_grid_0ka_global, Bioclimatic_variables_TRACE_I_ICE5G_grid_0ka_global, Bioclimatic_variables_TRACE_II_ICE5G_grid_0ka_global) 
  
  # Save the combined bioclimatic variables as a CSV file
  write.csv(Bioclimatic_variables_all_model_grid_0ka_global, file="result/Statistical analysis 4_result.3-Simulated bioclimatic variables at 0 ka BP for grid-cell.csv", row.names=FALSE)
  
}


# 4.3 Calculate bioclimatic variables of simulation at 0 cal. ka BP
# -------------------------------------------------------------------------------------------------

# Load datasets
Bioclimatic_variables_all_model_grid_0ka_global   <- read.csv2("result/Statistical analysis 4_result.3-Simulated bioclimatic variables at 0 ka BP for grid-cell.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimatic_variables_CRU_grid_0ka_global         <- read.csv2("result/Statistical analysis 4_result.2-Observed bioclimatic variables of 1931-1970 mean from CRU_grid cell (0.5 degrees).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_PNV_global_5min                            <- read.csv2("data/Statistical analysis 4_data.3-Modern potential natural biomes (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Land_sea_mask_0ka_all_model_global                <- read.csv2("data/Statistical analysis 4_data.4-Land-sea mask at 0ka.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Bioclimatic_variables_all_model_grid_0ka_global   <- type.convert(Bioclimatic_variables_all_model_grid_0ka_global, as.is = TRUE) 
Bioclimatic_variables_CRU_grid_0ka_global         <- type.convert(Bioclimatic_variables_CRU_grid_0ka_global, as.is = TRUE) 
Modern_PNV_global_5min                            <- type.convert(Modern_PNV_global_5min, as.is = TRUE) 
Land_sea_mask_0ka_all_model_global                <- type.convert(Land_sea_mask_0ka_all_model_global, as.is = TRUE) 

# Subset datasets for specific models
Cliamte_0ka_MPI_GLAC1D_global      <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "MPI-ESM_GLAC1D") %>% dplyr::select(-1) 
Cliamte_0ka_MPI_ICE6G_global       <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "MPI-ESM_ICE6G") %>% dplyr::select(-1)
Cliamte_0ka_CLIMBER_GLAC1D_global  <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "CLIMBER-X_GLAC1D") %>% dplyr::select(-1) 
Cliamte_0ka_CLIMBER_ICE6G_global   <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "CLIMBER-X_ICE6G") %>% dplyr::select(-1) 
Cliamte_0ka_TRACE_I_global         <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "TRACE-21K-I_ICE5G") %>% dplyr::select(-1) 
Cliamte_0ka_TRACE_II_global        <- subset(Bioclimatic_variables_all_model_grid_0ka_global, Source == "TRACE-21K-II_ICE5G") %>% dplyr::select(-1) 

# Subset the land-sea mask datasets for specific models
Land_sea_mask_0ka_MPI_GLAC1D_global     <- subset(Land_sea_mask_0ka_all_model_global, Source == "MPI-ESM_GLAC1D") %>% select(-4) 
Land_sea_mask_0ka_MPI_ICE6G_global      <- subset(Land_sea_mask_0ka_all_model_global, Source == "MPI-ESM_ICE6G") %>% select(-4) 
Land_sea_mask_0ka_CLIMBER_GLAC1D_global <- subset(Land_sea_mask_0ka_all_model_global, Source == "CLIMBER-X_GLAC1D") %>% select(-4) 
Land_sea_mask_0ka_CLIMBER_ICE6G_global  <- subset(Land_sea_mask_0ka_all_model_global, Source == "CLIMBER-X_ICE6G") %>% select(-4) 
Land_sea_mask_0ka_TRACE_I_global        <- subset(Land_sea_mask_0ka_all_model_global, Source == "TRACE-21K-I_ICE5G") %>% select(-4) 
Land_sea_mask_0ka_TRACE_II_global       <- subset(Land_sea_mask_0ka_all_model_global, Source == "TRACE-21K-II_ICE5G") %>% select(-4)


# Extract grid-cells that contain potential natural vegetation (PNV)
{
  # Extract unique coordinates from the observed bioclimatic variables
  coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global[ ,2:3])
  
  # Initialize a result matrix to store coordinates with PNV
  CRU_PNV_coordinate           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= 4, byrow=TRUE))
  colnames(CRU_PNV_coordinate) <- c(names(Bioclimatic_variables_CRU_grid_0ka_global)[2:3], "PNV", "PNV-grid_number")
  
  CRU_PNV_coordinate[ ,1:2] <- coordinate_list
  
  # Loop through each coordinate to identify PNV presence
  for (i in 1:nrow(coordinate_list)) {
    
    print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
    
    # Define a spatial subset for the grid cell
    Longitude_subset <- c(coordinate_list[i, 1] - 0.25,  coordinate_list[i, 1] + 0.25)
    Latitude_subset  <- c(coordinate_list[i, 2] - 0.25,  coordinate_list[i, 2] + 0.25)
    
    # Extract rows matching the grid cell from the PNV dataset
    subset_PNV <- subset(Modern_PNV_global_5min, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
    
    # Count the number of matching rows and determine PNV presence
    CRU_PNV_coordinate[i, 4] <- nrow(subset_PNV)
    
    if(nrow(subset_PNV) == 0){
      
      CRU_PNV_coordinate[i, 3] <- "No"
      
    }else{
      
      CRU_PNV_coordinate[i, 3] <- "Yes"
      
    }
    
  }
  
  # Filter coordinates with PNV present
  CRU_PNV_coordinate_final <- na.omit(CRU_PNV_coordinate[CRU_PNV_coordinate$PNV == "Yes", ])
  
  # Join the observed bioclimatic variables with the coordinates containing PNV
  Bioclimatic_variables_CRU_grid_0ka_global_final <- right_join(Bioclimatic_variables_CRU_grid_0ka_global, CRU_PNV_coordinate_final[ ,1:2], by = c("Longitude", "Latitude"))
  
  # Save the final dataset as a CSV file
  write.csv(Bioclimatic_variables_CRU_grid_0ka_global_final, file="result/Statistical analysis 4_result.4-Observed bioclimatic variables of 1931-1970 mean_CRU (0.5 degrees) with PNV.csv", row.names=FALSE)
  
}

# Bias calculation
{
  # (1) MPI-ESM_GLAC1D
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_MPI_GLAC1D_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_MPI_GLAC1D_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_MPI_GLAC1D_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_MPI_GLAC1D_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D) <- c(names(Cliamte_0ka_MPI_GLAC1D_global)[1:2], "Grid_number", names(Cliamte_0ka_MPI_GLAC1D_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias[ ,1]   <- "MPI-ESM_GLAC1D - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        MPI_GLAC1D_subset <- subset(Cliamte_0ka_MPI_GLAC1D_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(MPI_GLAC1D_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias[i, j+1] <- MPI_GLAC1D_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias), Land_sea_mask_0ka_MPI_GLAC1D_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (2) MPI-ESM_ICE6G
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_MPI_ICE6G_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_MPI_ICE6G_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_MPI_ICE6G_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_MPI_ICE6G_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G) <- c(names(Cliamte_0ka_MPI_ICE6G_global)[1:2], "Grid_number", names(Cliamte_0ka_MPI_ICE6G_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias[ ,1]   <- "MPI-ESM_ICE6G - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        MPI_ICE6G_subset <- subset(Cliamte_0ka_MPI_ICE6G_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(MPI_ICE6G_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias[i, j+1] <- MPI_ICE6G_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias), Land_sea_mask_0ka_MPI_ICE6G_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (3) CLIMBER-X_GLAC1D
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_CLIMBER_GLAC1D_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_CLIMBER_GLAC1D_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_CLIMBER_GLAC1D_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_CLIMBER_GLAC1D_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D) <- c(names(Cliamte_0ka_CLIMBER_GLAC1D_global)[1:2], "Grid_number", names(Cliamte_0ka_CLIMBER_GLAC1D_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias[ ,1]   <- "CLIMBER-X_GLAC1D - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        CLIMBER_GLAC1D_subset <- subset(Cliamte_0ka_CLIMBER_GLAC1D_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(CLIMBER_GLAC1D_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias[i, j+1] <- CLIMBER_GLAC1D_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias), Land_sea_mask_0ka_CLIMBER_GLAC1D_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (4) CLIMBER-X_ICE6G
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_CLIMBER_ICE6G_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_CLIMBER_ICE6G_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_CLIMBER_ICE6G_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_CLIMBER_ICE6G_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G) <- c(names(Cliamte_0ka_CLIMBER_ICE6G_global)[1:2], "Grid_number", names(Cliamte_0ka_CLIMBER_ICE6G_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias[ ,1]   <- "CLIMBER-X_ICE6G - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        CLIMBER_ICE6G_subset <- subset(Cliamte_0ka_CLIMBER_ICE6G_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(CLIMBER_ICE6G_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias[i, j+1] <- CLIMBER_ICE6G_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias), Land_sea_mask_0ka_CLIMBER_ICE6G_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (5) TRACE-21K-I_ICE5G
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_TRACE_I_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_TRACE_I_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_TRACE_I_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_TRACE_I_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I) <- c(names(Cliamte_0ka_TRACE_I_global)[1:2], "Grid_number", names(Cliamte_0ka_TRACE_I_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias[ ,1]   <- "TRACE-21K-I_ICE5G - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        TRACE_I_subset <- subset(Cliamte_0ka_TRACE_I_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(TRACE_I_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias[i, j+1] <- TRACE_I_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias), Land_sea_mask_0ka_TRACE_I_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (6) TRACE-21K-II_ICE5G
  {
    # Downscale for CRU Data
    {
      # Extract unique coordinates
      coordinate_list <- unique(Cliamte_0ka_TRACE_II_global[ ,1:2])
      
      # Get the number of unique longitude and latitude grid points
      Longitude_number <- length(unique(Cliamte_0ka_TRACE_II_global$Longitude))
      Latitude_number  <- length(unique(Cliamte_0ka_TRACE_II_global$Latitude))
      
      # Initialize a result matrix for downscaled CRU data
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Cliamte_0ka_TRACE_II_global)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II) <- c(names(Cliamte_0ka_TRACE_II_global)[1:2], "Grid_number", names(Cliamte_0ka_TRACE_II_global)[3:6])
      
      # Assign coordinates to the result matrix
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II[ ,1:2] <- coordinate_list
      
      # Loop through each grid cell
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Define the bounding box for the grid cell
        Longitude_subset <- c(coordinate_list[i, 1] - 360/Longitude_number/2, coordinate_list[i, 1] + 360/Longitude_number/2)
        Latitude_subset  <- c(coordinate_list[i, 2] - 180/Latitude_number/2,  coordinate_list[i, 2] + 180/Latitude_number/2)
        
        # Subset the CRU data for the current grid cell
        subset_grid <- subset(Bioclimatic_variables_CRU_grid_0ka_global_final, Longitude >= Longitude_subset[1] & Longitude <= Longitude_subset[2] & Latitude >= Latitude_subset[1] & Latitude <= Latitude_subset[2])
        
        # Record the number of overlapping CRU grid cells
        Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II[i, 3]    <- nrow(subset_grid)
        
        # If no CRU data exists for the grid cell, set values to NA
        if(nrow(subset_grid) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II[i, 4:7] <- NA
          
        }else{
          
          # Compute the median value for each bioclimatic variable in the subset
          for (j in 4:7) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II[i, j] <- median(subset_grid[ ,j])
            
          }
          
        }
        
      }
      
    }  
    
    
    # Bias calculation
    {
      # Remove rows with missing values and extract coordinates
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_final <- na.omit(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II[ ,c(1:2,4:7)])
      
      coordinate_list <- unique(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_final[ ,1:2])
      
      # Initialize a result matrix for bias calculations
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias           <- data.frame(matrix(NA, nrow= nrow(coordinate_list), ncol= ncol(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_final)+1, byrow=TRUE))
      colnames(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias) <- c("Source" ,names(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_final))
      
      # Set source label for bias matrix
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias[ ,1]   <- "TRACE-21K-II_ICE5G - CRU"
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias[ ,2:3] <- coordinate_list
      
      # Loop through each coordinate
      for (i in 1:nrow(coordinate_list)) {
        
        print(paste0("+++++ grid-cell ", i, "/", nrow(coordinate_list)," +++++"))
        
        # Subset model and CRU data for the current grid cell
        TRACE_II_subset <- subset(Cliamte_0ka_TRACE_II_global, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        CRU_subset        <- subset(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_final, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # If no data exists for either model or CRU, set bias values to NA
        if(nrow(TRACE_II_subset) == 0 | nrow(CRU_subset) == 0){
          
          Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias[i, 4:7] <- NA
          
        }else{
          
          # Calculate the bias as the difference between model and CRU values
          for (j in 3:6) {
            
            Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias[i, j+1] <- TRACE_II_subset[ ,j] - CRU_subset[ ,j]
            
          }
          
        }
        
      }
      
      # Join bias results with land-sea mask and remove rows with missing values
      Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias_final <- full_join(na.omit(Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias), Land_sea_mask_0ka_TRACE_II_global, by =c("Longitude", "Latitude")) %>%
        na.omit() %>%
        select(-8)
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # Combine datasets
  Bioclimatic_variables_CRU_grid_0ka_global_all_models_bias_final <- rbind(Bioclimatic_variables_CRU_grid_0ka_global_MPI_GLAC1D_bias_final, Bioclimatic_variables_CRU_grid_0ka_global_MPI_ICE6G_bias_final, 
                                                                           Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_GLAC1D_bias_final, Bioclimatic_variables_CRU_grid_0ka_global_CLIMBER_ICE6G_bias_final, 
                                                                           Bioclimatic_variables_CRU_grid_0ka_global_TRACE_I_bias_final, Bioclimatic_variables_CRU_grid_0ka_global_TRACE_II_bias_final)
  
  # Save the combined bioclimatic variables as a CSV file
  write.csv(Bioclimatic_variables_CRU_grid_0ka_global_all_models_bias_final, file="result/Statistical analysis 4_result.5-Bias of bioclimatic variables between Observation and simulation at 0 ka BP for grid-cell.csv", row.names=FALSE)
  
}



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 5. Statistical analysis 5: Space Constrained Clustering ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
EMD_spatial_pollen_simulations_grid_global_0_21ka_continent   <- read.csv2("data/Statistical analysis 5_data.1-Summary of EMD for grid by continent.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
EMD_spatial_pollen_simulations_grid_global_0_21ka_continent   <- type.convert(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, as.is = TRUE) 

# Calculation
{
  # (1) Europe
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "Europe")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_Europe <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_Europe <- rbind(EMD_spatial_pollen_simulations_grid_clustering_Europe, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_Europe_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_Europe, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Europe[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_Europe_final <- EMD_spatial_pollen_simulations_grid_clustering_Europe_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (2) Asia
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "Asia")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_Asia <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_Asia <- rbind(EMD_spatial_pollen_simulations_grid_clustering_Asia, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_Asia_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_Asia, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Asia[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_Asia_final <- EMD_spatial_pollen_simulations_grid_clustering_Asia_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  
  # (3) North America
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "North America")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_North_America <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_North_America <- rbind(EMD_spatial_pollen_simulations_grid_clustering_North_America, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_North_America_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_North_America, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_North_America[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_North_America_final <- EMD_spatial_pollen_simulations_grid_clustering_North_America_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  # (4) Africa
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "Africa")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_Africa <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_Africa <- rbind(EMD_spatial_pollen_simulations_grid_clustering_Africa, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_Africa_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_Africa, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Africa[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_Africa_final <- EMD_spatial_pollen_simulations_grid_clustering_Africa_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  # (5) Indopacific
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "Indopacific")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_Indopacific <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_Indopacific <- rbind(EMD_spatial_pollen_simulations_grid_clustering_Indopacific, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_Indopacific_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_Indopacific, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_Indopacific[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_Indopacific_final <- EMD_spatial_pollen_simulations_grid_clustering_Indopacific_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  # (6) South America
  {
    # Subset the dataset
    EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America <- subset(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent, Continent == "South America")
    
    # Compute the distance matrix using dynamic time warping (DTWARP) method
    dist_EMD_space <- TSclust::diss(SERIES = as.matrix(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[ , 3:8]), METHOD = "DTWARP") # note the dataframe needs to be converted to a matrix
    
    # Convert the distance matrix to an unnamed format
    dist_EMD_space <- unname(dist_EMD_space)
    
    # Generate a list of neighbors based on Delaunay triangulation
    listW <- nb2listw(tri2nb(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2], row.names = NULL), style="B")
    
    # Convert the neighbor list to a matrix
    links.mat.dat <- listw2mat(listW)
    
    # Extract neighbor pairs
    neighbors <- listw2sn(listW)[,1:2]
    
    # Plot the points representing the grid cells
    plot(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2], type='n',asp=1)
    title("Delaunay triangulation")
    text(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2], pos=3)  # Annotate the grid cells
    
    # Draw lines to show the connections (triangulation) between neighboring points
    for(i in 1:nrow(neighbors))
      lines(rbind(EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2][neighbors[i,1],],
                  EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2][neighbors[i,2],]))
    
    # Perform hierarchical clustering with contiguity constraints using Ward's method
    grpWD2cst_constr_hclust <- constr.hclust(dist_EMD_space, method="ward.D2",
                                             links = neighbors, coords = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[, 1:2], chron = F)
    
    ## Generic functions from hclust can be used, for instance to obtain a list of members of each cluster:
    {
      # Initialize a dataframe to store clustering results
      EMD_spatial_pollen_simulations_grid_clustering_South_America <- NULL
      
      # Loop through different numbers of clusters (k)
      for (i in 2:10){
        
        print(paste0("+++++ Grid-cell_based k = ",i, " +++++"))
        
        # Cut the hierarchical clustering tree to create `k` clusters
        hclus <- data.frame(Grid_ID = EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America$Grid_ID, Cluster = stats::cutree(grpWD2cst_constr_hclust, k=i)) 
        
        # Add the number of clusters (k) as a column
        hclus$K <- i
        
        # Append the results for the current `k` to the final clustering dataframe
        EMD_spatial_pollen_simulations_grid_clustering_South_America <- rbind(EMD_spatial_pollen_simulations_grid_clustering_South_America, hclus)
        
      }
      
      # Merge the clustering results with metadata about the grid cells
      EMD_spatial_pollen_simulations_grid_clustering_South_America_final <- left_join(EMD_spatial_pollen_simulations_grid_clustering_South_America, EMD_spatial_pollen_simulations_grid_global_0_21ka_continent_South_America[ ,c(1:2,9:10)], by = "Grid_ID")
      
      # Rearrange columns for better readability
      EMD_spatial_pollen_simulations_grid_clustering_South_America_final <- EMD_spatial_pollen_simulations_grid_clustering_South_America_final[ ,c(1,4:6,2:3)]
      
    }
    
  }
  # -------------------------------------------------------------------------------------------------
  
  # Combine datasets
  EMD_spatial_pollen_simulations_grid_clustering_global_final <- rbind(subset(EMD_spatial_pollen_simulations_grid_clustering_Europe_final, K == 2),
                                                                       subset(EMD_spatial_pollen_simulations_grid_clustering_Asia_final, K == 3),
                                                                       subset(EMD_spatial_pollen_simulations_grid_clustering_North_America_final, K == 3), 
                                                                       subset(EMD_spatial_pollen_simulations_grid_clustering_Africa_final, K == 2), 
                                                                       subset(EMD_spatial_pollen_simulations_grid_clustering_Indopacific_final, K == 2), 
                                                                       subset(EMD_spatial_pollen_simulations_grid_clustering_South_America_final, K == 3))           
  
  # Save the combined clustering result as a CSV file
  write.csv(EMD_spatial_pollen_simulations_grid_clustering_global_final, file="result/Statistical analysis 5_result.1-Clustering of data-model at each timeslice since 21 ka for grid by continent.csv", row.names=FALSE)
  
}



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 6. Statistical analysis 6: EMD between pollen-based reconstructions and ESM-based simulations ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka   <- read.csv2("data/Statistical analysis 6_data.1-Affinity score of reconstruction and simulation ensemble for site_Standardized.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka   <- type.convert(Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka, as.is = TRUE) 

# Define the Earth Mover's Distance (EMD) weight matrix for calculating similarity between categories
EMD_weights <- matrix(c(0,	1,	2,	3,	1,	2,	3,	4,
                        1,	0,	1,	2,	2,	2,	3,	3,
                        2,	1,	0,	1,	3,	2,	3,	2,
                        3,	2,	1,	0,	4,	3,	2,	1,
                        1,	2,	3,	4,	0,	1,	2,	4,
                        2,	2,	2,	3,	1,	0,	1,	1,
                        3,	3,	3,	2,	2,	1,	0,	1,
                        4,	3,	2,	1,	4,	1,	2,	0), 
                      ncol = 8, nrow = 8, byrow = TRUE) 


# Calculation
{
  # Add a new column for storing EMD values
  Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka$EMD <- apply(
    Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka, 
    1, 
    function(row) {
      # Print progress
      cat("Processing Dataset_ID:", row["Dataset_ID"], "\n")
      
      # Compute the Earth Mover's Distance (EMD)
      paleotools::EMD(
        as.numeric(row[5:12]),   # Reconstruction values
        as.numeric(row[13:20]), # Simulation values
        weight.m = EMD_weights
      )
    }
  )
  
  # Save the results to a CSV file
  # write.csv(Affinity_score_reconstructions_simulations_timeslice_standardized_global_0_21ka, file="result/Statistical analysis 6_result.1-EMD between pollen-based reconstructions and ESM-based simulations at timeslice for site.csv", row.names=FALSE)
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# END IN HERE #
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------





