library('tidyverse')
library('ggplot2')
library("readxl")
library("gganimate")
library("gifski")
library("lmtest")
library("hexbin")
library("qpcR")
library(RTransferEntropy)
library(lsa)
library(igraph)
library(dplyr)
# Set to where you've stored the data
setwd("~/Thesis/chennai data from Kanagaraj")

# Set plot theme
windowsFonts(name=windowsFont("Times New Roman"))
par(family = "Times New Roman")
theme_set(theme_minimal(base_size = 14) +
            theme(text = element_text(family = "Times New Roman"),
                  plot.title = element_text(size = 18, face = "bold"),
                  axis.title = element_text(size = 20, face = "bold"),
                  axis.text = element_text(size = 17),
                  axis.line = element_line(size = 1.5),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 16)))

# Read file to get data
dataset <- read_excel('ChennaiTrajectoryData2.45-3.00PM.xlsx')
testset <- read_excel('ChennaiTrajectoryData3.00-3.15PM.xlsx')


mean(testset$`Long Speed (m/sec)`)

###########################
#Revised Lateral Distance Plots

firstcar <- dataset %>% slice(which(dataset$'Vehicle Number' == 16 | dataset$'Vehicle Number' == 17)) 
firstcar
# Group by time steps
grouped_data <- firstcar %>%
  group_by(`Time (sec)`) 

# Calculate difference in Lat Distance per second
result <- grouped_data %>%
  mutate(diff_lat_distance = ifelse(n() == 2, (last(`Lat Distance (m)`) - first(`Lat Distance (m)`)), NA))

# Filter rows where there are two vehicles in the same time step
result <- result %>%
  filter(!is.na(diff_lat_distance))

# Set the first synchronous timestep to zero

min_time <- min(result$`Time (sec)`)
result$`Time (sec)` <- result$`Time (sec)` - min_time

# Calculate the maximum absolute value of the difference in Lat Distance per second
max_abs_diff <- max(abs(result$diff_lat_distance))

# Calculate the y-axis limits to center the plot on y=0
y_axis_limits <- c(-max_abs_diff, max_abs_diff)

#Plot difference in lateral distance over time

latdiffplot <- ggplot(result, aes(x = `Time (sec)`, y = diff_lat_distance)) +
  geom_line(color = "#0072B2", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Difference in Lateral Distance Between Vehicles over Time",
       x = "Time (sec)",
       y = "Difference in Lateral Distance (m)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  ylim(y_axis_limits)  # Set y-axis limits

latdiffplot

# Doing for all remaining pairs in 16 17, 16 18, 16 19, 17 18, 17 19, 18 19

selected_vehicles <- 16:19
filtered_data <- dataset %>%
  filter(`Vehicle Number` %in% selected_vehicles)

# Create combinations of vehicles
vehicle_combinations <- combn(selected_vehicles, 2)

# Initialize an empty list to store dataframes for each pair
results_list <- list()

for (i in 1:ncol(vehicle_combinations)) {
  vehicle_pair <- vehicle_combinations[, i]
  
  # Filter data for the current vehicle pair
  temp_data <- filtered_data %>%
    filter(`Vehicle Number` %in% vehicle_pair) %>%
    group_by(`Time (sec)`) %>%
    mutate(diff_lat_distance = ifelse(n() == 2, (last(`Lat Distance (m)`) - first(`Lat Distance (m)`)), NA))
  
  # Filter rows where there are two vehicles in the same time step
  temp_data <- temp_data %>%
    filter(!is.na(diff_lat_distance))
  
  # Create a column name for diff_lat_distance based on the vehicle pair
  pair_name <- paste(vehicle_pair, collapse = '-')
  
  # Rename the diff_lat_distance column using select and rename
  temp_data <- temp_data %>%
    dplyr::select(`Time (sec)`, !!pair_name := diff_lat_distance)
  
  # Store the results in a named list
  results_list[[pair_name]] <- temp_data
}

results_list

combined_data <- reduce(results_list, full_join, by = "Time (sec)")

custom_colors <- c("16-17" = "#1f77b4", "16-18" = "#ff7f0e", "16-19" = "#2ca02c",
                   "17-18" = "#d62728", "17-19" = "#9467bd", "18-19" = "#8c564b")

# Create the plot
gg <- ggplot(combined_data, aes(x = `Time (sec)`)) +
  geom_line(aes(y = `16-17`, color = "16-17"), size = 1.2) +
  geom_line(aes(y = `16-18`, color = "16-18"), size = 1.2) +
  geom_line(aes(y = `16-19`, color = "16-19"), size = 1.2) +
  geom_line(aes(y = `17-18`, color = "17-18"), size = 1.2) +
  geom_line(aes(y = `17-19`, color = "17-19"), size = 1.2) +
  geom_line(aes(y = `18-19`, color = "18-19"), size = 1.2) +
  labs(x = "Time (s)", y = "Lateral Distance Between Leader and Follower") +
  scale_color_manual(values = custom_colors, name = "Leader-Follower\nPairs") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5) +
  ylim(c(-5, 7)) +
  theme(legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12, face = "bold"),   # Legend text size
        panel.grid.major = element_blank(),      # Remove major gridlines
        panel.grid.minor = element_blank())      # Remove minor gridlines
gg

ggsave("latdist_plot.png", gg, width = 4, height = 6, dpi = 800)


############################

times <- function(df, vn) {
  index <- na.omit(df[c(which(df["Vehicle Number"]==vn)),]['Time (sec)'])
  return(index)
}
#1b. Get all other cars that share x timeslices with that car

firsts <- dataset %>% group_by(`Vehicle Number`)
firsts <- firsts %>% summarize(`Vehicle Number`, first = min(`Time (sec)`)) %>% unique()

vehicles <- function(df, vn) {
  last <- as.numeric(tail(times(df,vn),n=1))
  result <- which(firsts$`first`-last < -10 & firsts$`Vehicle Number` > vn)
  return(result)
}
#2. Run Granger tests on all reasonable orders for them, aligned to time
#2a. Make a vector for each and put dummy entries at the end

maketables <- function(df, vn){
  if(length(vehicles(df,vn))==0){
    print("Vehicle has no other eligible vehicles in sliding window",)
  }
  else{
  data <- data.frame(df$`Lat Distance (m)`[df$`Vehicle Number` == vn])
  names <- c(paste0("Vehicle", vn))
  for(i in vehicles(df, vn)) {                                   # Head of for-loop
    new <- df$`Lat Distance (m)`[df$`Vehicle Number` == i]                      # Create new column
    data <- qpcR:::cbind.na(data, new)               # Append new column
    names <- append(names, paste0("Vehicle", i))
  }
  colnames(data) <- names
  return(data)
  }
}
dataset

#2b. Run Granger test on each vector and output results

runGtest <- function(df, vn){
  table <- maketables(df, vn)
  data <- data.frame() %>% rbind.data.frame(c('Lead','Follower','leaddata','followdata','lag','pvalue'))
  for(i in 2:ncol(table)){
    for(j in 2:4){
      new <- grangertest(table[,1], table[,i], order = j)
      if(as.numeric(new$'Pr(>F)'[2])<0.001){
        data <- rbind.data.frame(data, c(colnames(table)[1], colnames(table)[i], new$'Res.Df',new$'Df'[2],new$'Pr(>F)'[2]))
      }
    }
  }
  return(data)
}

unique_vehicle_numbers <- unique(dataset$`Vehicle Number`)
vehicle_types <- vector("character", length(unique_vehicle_numbers))

vehicle_types <- sapply(unique_vehicle_numbers, function(num) {
  subset_data <- subset(dataset, `Vehicle Number` == num)
  if (nrow(subset_data) > 0) {
    return(subset_data$`Vehicle Type`[1])
  } else {
    return(NA)  # Return NA if the vehicle number is not found in the dataset
  }
})

vehicle_types[1+239]

results = data.frame()

# Loop through each vehicle number
for (vehicle_number in unique_vehicle_numbers) {
  # Filter the dataset for the current vehicle number
  
  test <- maketables(dataset, vehicle_number)
  if (is.character(test)){
  } else{
    for (n in 2:ncol(test)){
      acf_result <- acf(test[,c(1,n)], plot = FALSE, na.action = na.pass)
      acf_result_series <- acf_result$acf[, 1, 1]
      
      N <- acf_result$n.used
      
      threshold <- (exp(2*1.96/sqrt(N-3)-1))/(exp(2*1.96/sqrt(N-3)+1))
      
      # Initialize the significant_lag variable
      significant_lag <- 0
      
      # Find the indices of significant lags
      significant_indices <- which(abs(acf_result_series) > threshold)
      
      # Check if there are any significant lags
      if (length(significant_indices) > 0) {
        # Find the first non-significant lag
        first_non_significant <- min(which(abs(acf_result_series) <= threshold))
        
        # Check if a non-significant lag exists
        if (!is.na(first_non_significant)) {
          # Select the continuous significant lags up to the first non-significant lag
          continuous_significant <- significant_indices[significant_indices <= first_non_significant]
          
          # Find the highest index of continuous significant lag
          significant_lag <- max(continuous_significant)
        } else {
          # If there are no non-significant lags, set to the highest significant index
          significant_lag <- max(significant_indices)
        }
      }
      
      # Store the highest significant lag for the current vehicle
      results = rbind(results, c(significant_lag, vehicle_types[vehicle_number+n]))
    }
  }
}
results

laghist <- ggplot(data = results, aes(x = X13L/2 + 1)) +
  geom_histogram(binwidth = 1, fill = "light blue", color = "black", center = 5) +
  labs(
       x = "Time (s)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(0, max(results$X13L/2+1), by = 1), labels = function(x) as.integer(x))

ggsave("LagHist.png", laghist, width = 8, height = 6, dpi = 300)

results$X13

for(type in 1:6){
  laghist_granular <- ggplot(data = subset(results, `X1` == type), aes(x = X13/2 + 1)) +
    geom_histogram(binwidth = 1, fill = "light blue", color = "black", center = 5) +
    labs(
      x = "Time (s)",
      y = "Frequency"
    ) +
    scale_x_continuous(breaks = seq(0, max(results$X13), by = 1), labels = function(x) as.integer(x))
  ggsave(paste("LagHist", type,".png"), laghist_granular, width = 8, height = 6, dpi = 300)
}

laghist_granular
mean(results$X13L/2)

######################
# Transfer Entropy #

test <- maketables(testset, 1002)
test
nrow(test)

entropy_results <- data.frame(timestep = 1:50)
for (n in 2:ncol(test)){
  entropy <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("timestep", "TE"))))
  if(length(test[,n]>2)){
    for(i in 2:min(length(test[,n]),length(test[,1]))){
      entropy <- rbind(entropy, c(i, coef(RTransferEntropy::transfer_entropy(test[1:i,1], test[1:i,n],nboot=100),quiet = T)[1,2]))
    }
    col_name <- paste("TE -> ", colnames(test)[n],  sep = "")
    colnames(entropy) <- c("timestep", col_name)
    entropy_results <- merge(entropy_results, entropy, by = "timestep", all = T)
  }
}

df_filtered <- entropy_results[!is.na(entropy_results$timestep), ]

# Create a line plot
line_colors <- c("red", "blue", "green", "purple", "orange", "pink", "brown", "gray")
df_filtered

# Create a line plot with a for loop
p <- ggplot(df_filtered, aes(x = timestep / 2)) +
  labs(title = "Line Plot of TE vs. timestep",
       x = "Time (s)",
       y = "TE")

for (i in 2:ncol(df_filtered)) {
  col_name <- colnames(df_filtered)[i]
  p <- p + geom_line(aes_string(y = as.name(col_name)), color = line_colors[i - 1], size = 1.5)
}

p

df_filtered

df_long <- df_filtered %>%
  pivot_longer(cols = -timestep, names_to = "variable", values_to = "value")

# Create a faceted line plot
ggplot(df_long, aes(x = timestep / 2, y = value)) +
  geom_line(color = "lightblue", size = 1.5) +
  labs(title = "Line Plot of TE vs. timestep",
       x = "Time (s)",
       y = "TE") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(1.5, "lines"))

base_size <- 12  # Adjust this as needed

# Increase the size of the text elements
plot <- ggplot(df_long, aes(x = timestep / 2, y = value)) +
  geom_line(color = "lightblue", size = 1.5) +
  labs(
       x = "Time (s)",
       y = "TE") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = base_size) +  # Use a minimal theme
  theme(
    text = element_text(size = base_size),  # Increase text size
    plot.title = element_text(size = base_size + 2, face = "bold"),  # Title text size and bold
    axis.title = element_text(size = base_size + 1),  # Axis title size
    strip.text = element_text(size = base_size + 1),  # Facet label size
    panel.border = element_rect(color = "black", fill = NA),
    panel.spacing = unit(1.5, "lines")
  )

plot

# Print the plot
ggsave("TEfacet.png", plot, width = 8, height = 6, dpi = 300)


############### LOCAL TRANSFER ENTROPY USING DEPRECATED RINFORM LIBRARY
library("rinform")
rinform::transfer_entropy(
  na.omit(test$Vehicle1006[1:28]),
  na.omit(test$Vehicle1002[1:28]),
  k = 1,
  local = TRUE
)



rinform::transfer_entropy(
  na.omit(test$Vehicle1006[1:28]),
  na.omit(test$Vehicle1002[1:28]),
  k = 1,
  local = TRUE
)


test[,1]


runLocalEtest <- function(df, vn){
  entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())
  test <- maketables(df, vn)
  main_name <- paste('Vehicle',vn,sep ="")
  for (n in 2:ncol(test)) {
  col_name <- names(test)[n]  # Get the column name
  len <- min(length(na.omit(test[,1])), length(na.omit(test[, n])))
    # Calculate entropy
    entropy <- rinform::transfer_entropy(
      na.omit(test[,1][1:len]),
      na.omit(test[,n][1:len]),
      k = 1,
      local = TRUE
    )
    
    # Create a data frame for the current variable
    var_df <- data.frame(
      Vehicle = main_name,
      Variable = col_name,
      Entropy = entropy
    )
    
    # Append to the main data frame
    entropy_df <- bind_rows(entropy_df, var_df)
  }
  return(entropy_df)
  }


runLocalEtestAll <- function(df) {
  entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())
  
  for (vn in unique(df$`Vehicle Number`)) {
    test <- maketables(df, vn)
    main_name <- paste('Vehicle', vn, sep = "")
    if(is.data.frame(test)){
      for (n in 2:ncol(test)) {
        col_name <- names(test)[n]  # Get the column name
        len <- min(length(na.omit(test[, 1])), length(na.omit(test[, n])))
        
        # Calculate entropy
        entropy <- rinform::transfer_entropy(
          na.omit(test[, 1][1:len]),
          na.omit(test[, n][1:len]),
          k = 1,
          local = TRUE
        )
        
        # Create a data frame for the current variable
        var_df <- data.frame(
          Vehicle = main_name,
          Variable = col_name,
          Entropy = entropy
        )
        
        # Append to the main data frame
        entropy_df <- rbind(entropy_df, var_df)
      }
    }
  }
  
  return(entropy_df)
}

result <- runLocalEtestAll(testset)

filtered_result <- result[result$Entropy > 0.1, ]

# Get all unique pairs of Vehicle values
pairs <- unique(cbind(filtered_result$Vehicle, filtered_result$Variable))

# Calculate the duration of each pair by counting the rows
pair_durations <- sapply(1:nrow(pairs), function(i) {
  pair <- pairs[i, ]
  duration <- nrow(filtered_result[filtered_result$Vehicle == pair[1] & filtered_result$Variable == pair[2], ])
  return(duration)
})

# Create a data frame with pairs and their durations
result_df <- data.frame(Pair = sapply(1:nrow(pairs), function(i) {
  pair <- pairs[i, ]
  return(paste(pair[1], "-", pair[2]))
}), Duration = pair_durations)

mean(result_df$Duration)

num_bins <-  max(result_df$Duration)/2
x_range <- c(0, max(result_df$Duration)/2)

breaks <- seq(min(result_df$Duration), max(result_df$Duration), length.out = num_bins + 1)

# Create a histogram
hist_data <- hist(result_df$Duration, breaks = breaks, main = "Duration Histogram", xlab = "Duration", ylab = "Frequency", col = "blue", border = "black")

# Plot the histogram
durationhistplot<- plot(hist_data, col = "light blue", border = "black", main = "Histogram of Duration of Leader-Follower Relationships", xlab = "Duration of Leader-Follower Relationship (s)", ylab = "Frequency")
  
  x_ticks <- seq(0, max(result_df$Duration), 1)
  axis(1, at = x_ticks, labels = x_ticks)

  ggsave("duration_histogram.png", plot = durationhistplot, width = 8, height = 6, dpi = 800)
  

# Create a histogram
hist_data <- hist(result_df$Duration/2, breaks = num_bins, plot = FALSE)

# Plot the histogram with more frequent x-axis ticks and labels
hist_plot <- barplot(hist_data$counts, xlim = x_range, col = "light blue", border = "black", 
                     main = "Duration Histogram", xlab = "Duration of Leader-Follower Relationship (s)", ylab = "Frequency")
# Define the positions for x-axis ticks and labels
x_ticks <- seq(0, max(result_df$Duration), 1)  # You can adjust the "by" value for the desired frequency


# Add x-axis ticks and labels
axis(1, at = x_ticks, labels = x_ticks)

# Add a grid
grid()

test <- maketables(dataset, 1002)

testset

# Create an empty data frame to store the entropy values
entropy_df <- data.frame(Vehicle = character(), Variable = character(), Entropy = numeric())

# Iterate over columns of the test dataframe (starting from column 2)
for (n in 2:ncol(test)) {
  col_name <- names(test)[n]  # Get the column name
  len <- min(length(na.omit(test$Vehicle1002)), length(na.omit(test[, n])))
  
  # Calculate entropy
  entropy <- rinform::transfer_entropy(
    na.omit(test$Vehicle1002)[1:len],
    na.omit(test[, n])[1:len],
    k = 1,
    local = TRUE
  )
  
  # Create a data frame for the current variable
  var_df <- data.frame(
    Vehicle = "Vehicle1002",
    Variable = col_name,
    Entropy = entropy
  )
  
  # Append to the main data frame
  entropy_df <- bind_rows(entropy_df, var_df)
}

entropy_df <- entropy_df %>%
  group_by(Vehicle) %>%
  mutate(Index = row_number()%%11)

# Create a ggplot object
p <- ggplot(entropy_df, aes(x = Index, y = Entropy)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_smooth(color = "light blue", size = 1.25, se = FALSE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(x = "Time (s)", y = "Local Transfer Entropy from Vehicle 1002 to Candidate Followers") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Adjust size here
    panel.spacing = unit(1, "lines")
  ) + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10))
p
entropy_df

p

shaded_regions_df <- entropy_df %>%
  group_by(Variable) %>%
  filter(any(Entropy > 0.1))

# Create the plot with shaded region

p <- ggplot(entropy_df, aes(x = Index, y = Entropy)) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_smooth(color = "light blue", size = 1.25, se = FALSE) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(x = "Time (s)", y = "Local Transfer Entropy from Vehicle 1002 to Candidate Followers") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines")
  ) + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  geom_rect(
    aes(xmin = 2, xmax = 6.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1003")
  ) +
  geom_rect(
    aes(xmin = 2, xmax = 6.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1004")
  ) +
  geom_rect(
    aes(xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1005")
  ) +
  geom_rect(
    aes(xmin = 0.5, xmax = 6, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1006")
  ) +
  geom_rect(
    aes(xmin = 2, xmax = 9, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1007")
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1008")
  ) +
  geom_rect(
    aes(xmin = 9.5, xmax = 10, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1008")
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 4.5, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1009")
  ) +
  geom_rect(
    aes(xmin = 7, xmax = 10, ymin = -Inf, ymax = Inf),
    fill = "gray70", alpha = 0.01,
    data = subset(entropy_df, `Variable` == "Vehicle1009")
  ) +
  coord_cartesian(ylim = c(-0.05, 0.4))  # Set the y-axis limits
  
# Print the plot
print(p)
ggsave("TEfacet.png", p, width = 8, height = 8, dpi = 800)


##################

library('tidyverse')
library('ggplot2')
library("readxl")
library("gganimate")
library("gifski")
library("lmtest")
library("hexbin")
library("qpcR")
library(RTransferEntropy)
library(lsa)
library(igraph)
library(dplyr)
library(deSolve)
library(viridis)
setwd("~/Thesis/chennai data from Kanagaraj")

# Set plot theme
par(family = "Times New Roman")
theme_set(theme_minimal(base_size = 14) +
            theme(text = element_text(family = "Times New Roman"),
                  plot.title = element_text(size = 18, face = "bold"),
                  axis.title = element_text(size = 20, face = "bold"),
                  axis.text = element_text(size = 18),
                  axis.line = element_line(size = 1.5),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 16)))

# Read file to get data
dataset <- read_excel('ChennaiTrajectoryData2.45-3.00PM.xlsx')

## Trajectory Plots
subset_dataset <- dataset %>%
  filter(`Vehicle Number` %in% 16:19)

plot <- ggplot(subset_dataset, aes(y = `Long Distance (m)`, x = `Time (sec)`)) +
  geom_line(aes(group = `Vehicle Number`, color = paste("Vehicle", `Vehicle Number`)), size = 1.5) +
  labs(y = "Space (m)", x = "Time (s)") +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "Vehicle ID") +
  theme(legend.title = element_text(size = 16, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12, face = "bold"),   # Legend text size
        panel.grid.major = element_blank(),      # Remove major gridlines
        panel.grid.minor = element_blank())      # Remove minor gridlines

ggsave("Trajectory_Plot.png", plot, width = 4, height = 6, dpi = 800)
plot
