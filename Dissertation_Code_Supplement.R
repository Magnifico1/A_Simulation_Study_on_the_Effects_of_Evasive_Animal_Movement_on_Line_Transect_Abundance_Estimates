#-------------------------------------------------------------------------------
# A SIMULATION STUDY ON THE EFFECT OF EVASIVE ANIMAL MOVEMENT ON LINE TRANSECT -
# --------------------------   ABUNDANCE ESTIMATES   ---------------------------
#-------------------------------------------------------------------------------

# Author: Alexander Ross
# Last updated: 29th August 2023

# Load packages
library(mrds)
library(latex2exp)

# Define general survey parameters
w <- 200    # width of survey area
L <- 2000   # length of survey area
n <- 1000   # number of detectable animals in survey area

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#----------------------------                          -------------------------
#---------------------------   SINGLE OBSERVER (Boat)   ------------------------
#----------------------------                          -------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# § Simulating Uniform Distribution for Single Observer ------------------------

# Single Survey

set.seed(2000)

# Randomly generate animals in survey area from uniform distribution
x <- runif(n=n, min=0, max=L)  # Distance along transect 
y <- runif(n=n, min=-w, max=w) # Perpendicular distance from transect

# Coördinates of randomly-generated animals
observations_unif <- matrix(c(x,y), nrow=2, ncol=n, byrow=TRUE, 
                            dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate detection probabilities for animals based on half-normal function
p_line <- 0.8    # probability of boat detecting animals on the transect line
σ <- sqrt(-w^2/(2*log(0.1/p_line)))    # Scale parameter giving p(w)=0.1
probs_unif <- p_line*exp(-y^2/(2*σ^2))

# Generate detections from Bernoulli trials
detections_unif <- rbinom(n=n, size=1, prob=probs_unif)
# Detected animals' distances
data_unif <- y[which(detections_unif==1)]

# Plot histogram of detection distances
hist(abs(data_unif), xlab="Detection Distance (m)", 
     main="Histogram of Detection Distances\nunder Uniform Distribution",
     col="black", density=25, angle=60)

# Plot detected animals in red
plot(x, y, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="Uniform Distribution Detections by
     Single Observer")
abline(0, 0, lty=3, col="red")
points(x[detections_unif==1], y[detections_unif==1], col="blue")

# Data from uniform distribution 
detection_data_unif <- data.frame("object"=seq(1:n), 
                                  "observer"=rep(1,n), 
                                  "distance"=abs(y),
                                  "p_detection"=probs_unif,
                                  "detected"=detections_unif)

# Fit detection function with distance sampling method and half-normal key
ds_unif <- ddf(method='ds', dsmodel=~cds(key='hn'), 
               data=detection_data_unif, meta.data=list(width=w))

# Produce goodness of fit statistics and a qq plot
gof_unif <- ddf.gof(ds_unif, main="Goodness of Fit under Distance Sampling
                                   Method for Uniform Distribution")
chisq.p_unif <- gof_unif$chisquare$chi1$p

# Obtain estimates of abundance and its standard error
region <- data.frame("Region.Label"=1, 
                     "Area"=2*w)
samples <- data.frame("Sample.Label"=1, 
                      "Region.Label"=1, 
                      "Effort"=1)
obs_unif <- data.frame("object"=unique(detection_data_unif$object),
                       "Region.Label"=1, 
                       "Sample.Label"=1)
abundance_unif <- dht(model=ds_unif, region.table=region, sample.table=samples, 
                      obs.table=obs_unif)
Nhat_unif <- abundance_unif$individuals$N$Estimate
Nhat.se_unif <- abundance_unif$individuals$N$se



# § Simulating Evasive Movement for Single Observer ----------------------------

# Single survey

set.seed(2000)

# Define responsive movement parameters
prop <- 0.95          # proportion of all animals responding
z <- 0.5*w            # distance from line at which animals do not respond
ρ <- 10               # trackline evasive movement distance

# Randomly generate animals showing evasive movement
x_rm <- runif(n=n, min=0, max=L)
y <- runif(n=n, min=-w, max=w)
indices_prop <- sample(seq_along(y), size=prop*length(y))
y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                  y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
y_rm <- y
y_rm[indices_prop] <- y_shift[indices_prop]

# Plot distances moved in response
plot(abs(y), abs(y_shift-y)*20/ρ, ylim=c(0,20),
     xlab="Original Distance from Transect Line (m)", 
     ylab="Distance Moved in Response (m)", main="Evasive Movement Profile",
     type="p", col="darkblue")
points(abs(y), abs(y_shift-y)*1/ρ, col="cyan", type="p")
points(abs(y), abs(y_shift-y)*5/ρ, col="deepskyblue", type="p")
points(abs(y), abs(y_shift-y)*10/ρ, col="dodgerblue1", type="p")
points(abs(y), abs(y_shift-y)*15/ρ, col="mediumblue", type="p")
legend("topright", legend=c("ρ=20", "ρ=15", "ρ=10", "ρ=5", "ρ=1"),
       col=c("darkblue", "mediumblue", "dodgerblue1", "deepskyblue", "cyan"), 
       lty=1, lwd=1)

# Coördinates of randomly-generated animals
observations_rm <- matrix(c(x_rm,y_rm), nrow=2, ncol=n, byrow=TRUE, 
                          dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate detection probabilities based on half-normal function
probs_rm <- p_line*exp(-y_rm^2/(2*σ^2))
# Generate detections from Bernoulli trials
detections_rm <- rbinom(n=n, size=1, prob=probs_rm)
# Detected animals' distances
data_rm <- y_rm[which(detections_rm==1)]

# Difference in number of animals detected after evasive movement
diff_data <- length(data_unif) - length(data_rm)

# Plot histogram of detection distances
hist(abs(data_rm), xlab="Detection Distance (m)", 
     main=c("Histogram of Detection Distances under\nEvasive Movement (ρ=",
            ρ, ")"),
     col="black", density=25, angle=60)

# Plot detected animals in red
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main=c("Evasive Movement Detections 
           by Single Observer (ρ=", ρ, ")"))
abline(0, 0, lty=3, col="red")
points(x_rm[detections_rm==1], y_rm[detections_rm==1], col="blue")

# Data from evasive movement
detection_data_rm <- data.frame("object"=seq(1:n), 
                                "distance"=abs(y_rm), 
                                "p_detection"=probs_rm, 
                                "detected"=detections_rm)

# Fit detection function with distance sampling method and half-normal key
ds_rm <- ddf(method='ds', dsmodel=~cds(key='hn'), 
             data=detection_data_rm, meta.data=list(width=w))

# Produce goodness of fit statistics and a qq plot
gof_rm <- ddf.gof(ds_rm, main="Goodness of Fit under Distance Sampling
                  Method for Evasive Movement")
chisq.p_rm <- gof_rm$chisquare$chi1$p

# Obtain estimates of abundance and its standard error
obs_rm <- data.frame("object"=unique(detection_data_rm$object),
                     "Region.Label"=1, 
                     "Sample.Label"=1)
abundance_rm <- dht(model=ds_rm, region.table=region, sample.table=samples, 
                    obs.table=obs_rm)
Nhat_rm <- abundance_rm$individuals$N$Estimate
Nhat.se_rm <- abundance_rm$individuals$N$se


# § Simulating Multiple Surveys for varying Evasive Movement -------------------

set.seed(2000)

# Define range of responsive movement parameter
ρ_values <- seq(from=0, to=20, by=2)

# Number of simulated surveys per evasive movement value
iterations <- 200

# Create empty vectors to store results
mean_Nhats_rm <- numeric(length(ρ_values))
mean_Nhats.se_rm <- numeric(length(ρ_values))
mean_chisq.ps_rm <- numeric(length(ρ_values))

# Loop across evasive movement values
for (i in 1:length(ρ_values)) {
  
  # Create an empty vector to store the results for each ρ
  iteration_Nhats <- numeric(iterations)
  iteration_Nhats.se <- numeric(iterations)
  iteration_chisq.ps <- numeric(iterations)
  
  # Set ρ value
  ρ <- ρ_values[i]
  
  # Loop across surveys
  for (j in 1:iterations) {
    
    # Simulate detections in each survey as above
    x_rm <- runif(n=n, min=0, max=L)
    y <- runif(n=n, min=-w, max=w)
    indices_prop <- sample(seq_along(y), size=round(prop*length(y)))
    y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                      y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
    y_rm <- y
    y_rm[indices_prop] <- y_shift[indices_prop]
    probs_rm <- p_line*exp(-y_rm^2/(2*σ^2))
    detections_rm <- rbinom(n=n, size=1, prob=probs_rm)
    detection_data_rm <- data.frame("object"=seq(1:n), 
                                    "distance"=abs(y_rm), 
                                    "p_detection"=probs_rm,
                                    "detected"=detections_rm)
    
    # Fit the detection function
    ds_rm <- ddf(method='ds', dsmodel=~cds(key='hn'), 
                 data=detection_data_rm, meta.data=list(width=w))
    
    # Produce goodness of fit statistics
    gof_rm <- ddf.gof(ds_rm, qq=FALSE)
    chisq.p <- gof_rm$chisquare$chi1$p
    
    # Obtain estimates of abundance and its standard error
    obs_rm <- data.frame("object"=unique(detection_data_rm$object),
                         "Region.Label"=1, 
                         "Sample.Label"=1)
    abundance_rm <- dht(model=ds_rm, region.table=region, sample.table=samples, 
                        obs.table=obs_rm)
    Nhat_rm <- abundance_rm$individuals$N$Estimate
    Nhat.se_rm <- abundance_rm$individuals$N$se
    
    # Store estimated statistics
    iteration_Nhats[j] <- Nhat_rm
    iteration_Nhats.se[j] <- Nhat.se_rm
    iteration_chisq.ps[j] <- chisq.p
  }
  
  # Store the mean and standard error
  mean_Nhats_rm[i] <- mean(iteration_Nhats)
  mean_Nhats.se_rm[i] <- mean(iteration_Nhats.se)
  mean_chisq.ps_rm[i] <- mean(iteration_chisq.ps)
}

# Mean, standard error, and chi-squared p-value for abundance estimates 
for (i in 1:length(ρ_values)) {
  ρ <- ρ_values[i]
  cat("Results for ρ =", ρ, "\n")
  print(paste("Mean Abundance:", mean_Nhats_rm[i]))
  print(paste("Standard Error:", mean_Nhats.se_rm[i]))
  print(paste("χ-squared p-value:", mean_chisq.ps_rm[i]))
}

# Plot abundance results on numerical scale
plot(ρ_values, mean_Nhats_rm, pch=16, xlab="ρ", ylab=TeX(r'($\hat{N}$)'),
     ylim=c(min(mean_Nhats_rm - mean_Nhats.se_rm),
            ifelse(max(mean_Nhats_rm + mean_Nhats.se_rm)>n, 
                   max(mean_Nhats_rm + mean_Nhats.se_rm), n)),
     main="Evasive Movement Effect on 
     Estimated Abundance by Single Observer")
abline(n, 0, lty=2, col="red")
abline(lm(mean_Nhats_rm ~ ρ_values), col="blue")
abline(mean_Nhats_rm[1], 0, lty=2)
arrows(ρ_values, mean_Nhats_rm - mean_Nhats.se_rm, 
       ρ_values, mean_Nhats_rm + mean_Nhats.se_rm, angle=90, code=3,length=0.05)

# Plot abundance results on bias scale
bias_single <- 100*(mean_Nhats_rm - n)/n
plot(ρ_values, bias_single, pch=16, xlab="ρ", 
     ylab=TeX(r'($\%$ Bias in $\hat{N}$)'), 
     ylim=c((100*(min(mean_Nhats_rm - mean_Nhats.se_rm) - n)/n), 
            ifelse((100*(max(mean_Nhats_rm + mean_Nhats.se_rm) - n)/n)>0, 
                   (100*(max(mean_Nhats_rm + mean_Nhats.se_rm) - n)/n), 0)),
     main="Evasive Movement Effect on 
     Estimated Abundance by Single Observer")
abline(0,0,lty=2, col="red")
abline(lm(bias_single ~ ρ_values), col="blue")
abline(bias_single[1], 0, lty=2)
arrows(ρ_values, 100*(mean_Nhats_rm - mean_Nhats.se_rm - n)/n, 
       ρ_values, 100*(mean_Nhats_rm + mean_Nhats.se_rm - n)/n, 
       angle=90, code=3, length=0.05)

# Quantification of relationship between abundance and evasive movement
fit_single <- lm(mean_Nhats_rm ~ ρ_values)
summary(fit_single)
coefs_single <- fit_single$coefficients
a_s <- coefs_single[[1]]
b_s <- coefs_single[[2]]



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#----------                                                        -------------
#---------   DOUBLE OBSERVER (Boat Observer 1 & Drone Observer 2)   ------------
#----------                                                        -------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# INDEPENDENCE CONFIGURATION

# § Simulating Uniform Distribution for Double Observers -----------------------

# Single survey

set.seed(2000)

# Randomly generate animals
x <- runif(n=n, min=0, max=L)
y <- runif(n=n, min=-w, max=w)

# Coördinates of randomly-generated animals
observations_unif <- matrix(c(x,y), nrow=2, ncol=n, byrow=TRUE, 
                            dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate boat detection probabilities based on half-normal function
σ1 <- sqrt(-w^2/(2*log(0.1/p_line)))
probs_unif1 <- p_line*exp(-y^2/(2*σ1^2))
# Generate boat detections from Bernoulli trials
detections_unif1 <- rbinom(n=n, size=1, prob=probs_unif1)

# Find half-normal scale parameter giving effective strip half-width of 28m
hn_eshw <- function(hn_scale, w) {
  # Calculate area under half normal to truncation distance w
  hn_area <- 1 - pnorm(w, 0, hn_scale, lower.tail = FALSE)*2
  # Calculate effective strip half-width
  eshw <- hn_area / dnorm(0, 0, hn_scale)*0.5
  return(eshw)
}
# Find SSE
obj.func <- function(hn_scale, target_eshw, w) {
  eshw <- hn_eshw(hn_scale, w)
  return((eshw - target_eshw)^2)
}
target_eshw <- 28
res <- optimize(obj.func, c(0, w), target_eshw, w)
hn_scale <- res$minimum

# Calculate drone detection probabilities based on half-normal function
σ2 <- round(hn_scale, digits=2) 
probs_unif2 <- exp(-y^2/(2*σ2^2))
# Generate drone detections from Bernoulli trials
detections_unif2 <- rbinom(n=n, size=1, prob=probs_unif2)

# Plot uniformly-distributed animals with boat-detected animals coloured
par(mfrow=c(1,2))
plot(x, y, xlab="Distance along Transect (m)",ylab="Perpendicular distance (m)", 
     main="UD Boat Detections")
abline(0, 0, lty=3, col="red")
points(x[detections_unif1==1], y[detections_unif1==1], col="blue")

# Plot uniformly-distributed  animals with drone-detected animals coloured
plot(x, y, xlab="Distance along Transect (m)",ylab="Perpendicular distance (m)", 
     main="UD Drone Detections")
abline(0, 0, lty=3, col="red")
points(x[detections_unif2==1], y[detections_unif2==1], col="red")

# Plot animals with those detected by both observers coloured 
par(mfrow=c(1,1))
plot(x, y, xlab="Distance along Transect (m)",ylab="Perpendicular distance (m)", 
     main="Uniform Distribution Duplicates")
abline(0, 0, lty=3, col="red")
points(x[detections_unif1==1 & detections_unif2==1], 
       y[detections_unif1==1 & detections_unif2==1], col="purple")

# Data from uniform distribution
detection_data_unif <- data.frame(
                        "object"=rep(1:n, each=2), 
                        "observer"=as.factor(rep(1:2, times=n)), 
                        "distance"=rep(abs(y), each=2),
                        "p_detection"=c(rbind(probs_unif1, probs_unif2)),
                        "detected"=c(rbind(detections_unif1, detections_unif2)))

# Remove animals not detected by either observer from dataframe
zeros_unif <- detection_data_unif$detected==0
frequency_unif <- table(detection_data_unif$object[zeros_unif])
non_detections_unif <- as.numeric(names(frequency_unif[frequency_unif>1]))
detections_unif <- detection_data_unif[!(detection_data_unif$object %in% 
                                           non_detections_unif),]

# Fit detection function under independence configuration and point independence
detection_func_unif <- ddf(method='io', 
                           mrmodel=~glm(link='logit', 
                                        formula=~distance + observer),
                           dsmodel=~cds(key='hn'), data=detections_unif, 
                           meta.data=list(width=w))

# Plot histograms of double-observer detections
detection_tables_unif <- det.tables(detection_func_unif)
print(detection_tables_unif)
plot(detection_tables_unif, col2="blue")

# Plot unconditional and conditional detection functions
plot(detection_func_unif, xlab="Perpendicular Distance (m)",
     col="black", angle=60, density=25)

# Produce goodness of fit statistics and qq plot
gof_unif <- ddf.gof(detection_func_unif, 
                    main="Goodness of Fit for Uniform Distribution 
                    (Independence)")
chisq.p_unif <- gof_unif$chisquare$chi1$p

# Extract estimate of primary detection probability on transect line
summary_unif <- summary(detection_func_unif)
p.0_unif <- summary_unif$mr.summary$average.p0.1
p.0.se_unif <- summary_unif$mr.summary$average.p0.1.se

# Obtain estimates of abundance and its standard error
obs_unif <- data.frame("object"=unique(detections_unif$object),
                       "Region.Label"=1, 
                       "Sample.Label"=1)
abundance_unif <- dht(model=detection_func_unif, region.table=region, 
                      sample.table=samples, obs.table=obs_unif)
Nhat_unif <- abundance_unif$individuals$N$Estimate
Nhat.se_unif <- abundance_unif$individuals$N$se



# § Simulating Evasive Movement for Double Observers ---------------------------

# Single survey

set.seed(2000)

# Define responsive movement parameter
ρ <- 20   # Example value

# Randomly generate animals showing evasive movement
x_rm <- runif(n=n, min=0, max=L)
y <- runif(n=n, min=-w, max=w)
indices_prop <- sample(seq_along(y), size=round(prop*length(y)))
y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                  y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
y_rm <- y
y_rm[indices_prop] <- y_shift[indices_prop]

# Coördinates of randomly-generated animals
observations_rm <- matrix(c(x_rm, y_rm),nrow=2,ncol=n, byrow=TRUE, 
                 dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate boat detection probabilities based on half-normal function
probs_rm1 <- p_line*exp(-y_rm^2/(2*σ1^2))
# Generate boat detections from Bernoulli trials
detections_rm1 <- rbinom(n=n, size=1, prob=probs_rm1)
# Detected animals by boat
data_rm1 <- y_rm[which(detections_rm1==1)]

# Calculate drone detection probabilities based on half-normal function
probs_rm2 <- exp(-y_rm^2/(2*σ2^2))
# Generate drone detections from Bernoulli trials
detections_rm2 <- rbinom(n=n, size=1, prob=probs_rm2)
# Detected animals by drone
data_rm2 <- y_rm[which(detections_rm2==1)]

# Plot evasive movement animals with boat-detected animals coloured
par(mfrow=c(1,2))
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="EM Boat Detections")
abline(0,0, lty=3, col="red")
points(x_rm[detections_rm1==1], y_rm[detections_rm1==1], col="blue")

# Plot evasive movement animals with drone-detected animals coloured
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="EM Drone Detections")
abline(0,0, lty=3, col="red")
points(x_rm[detections_rm2==1], y_rm[detections_rm2==1], col="red")

# Plot animals with those detected by both observers coloured 
par(mfrow=c(1,1))
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="Evasive Movement Duplicates")
abline(0,0, lty=3, col="red")
points(x_rm[detections_rm1==1 & detections_rm2==1], 
       y_rm[detections_rm1==1 & detections_rm2==1], col="purple")

# Data from uniform distribution
detection_data_rm <- data.frame(
                            "object"=rep(1:n, each=2), 
                            "observer"=as.factor(rep(1:2, times=n)), 
                            "distance"=rep(abs(y_rm), each=2),
                            "p_detection"=c(rbind(probs_rm1, probs_rm2)),
                            "detected"=c(rbind(detections_rm1, detections_rm2)))

# Remove animals not detected by either observer from dataframe
zeros_rm <- detection_data_rm$detected==0
frequency_rm <- table(detection_data_rm$object[zeros_rm])
non_detections_rm <- as.numeric(names(frequency_rm[frequency_rm>1]))
detections_rm <- detection_data_rm[!(detection_data_rm$object %in% 
                                       non_detections_rm),]

# Fit detection function under independence configuration and point independence
detection_func_rm <- ddf(method='io', 
                         mrmodel=~glm(link='logit', 
                                      formula=~distance + observer),
                         dsmodel=~cds(key='hn'), data=detections_rm, 
                         meta.data=list(width=w))

# Plot histograms of double-observer detections
detection_tables_rm <- det.tables(detection_func_rm)
print(detection_tables_rm)
par(mfrow=c(1,1))
plot(detection_tables_rm, col2="blue")

# Plot unconditional and conditional detection functions
plot(detection_func_rm, xlab="Perpendicular Distance (m)",
     col="black", angle=60, density=25)

# Produce goodness of fit statistics and qq plot
gof_rm <- ddf.gof(detection_func_rm, 
                  main="Goodness of Fit for Evasive Movement
                    (Independence)")
chisq.p_rm <- gof_rm$chisquare$chi1$p

# Extract estimate of primary detection probability on transect line
summary_rm <- summary(detection_func_rm)
p.0_rm <- summary_rm$mr.summary$average.p0.1
p.0.se_rm <- summary_rm$mr.summary$average.p0.1.se

# Obtain estimates of abundance and its standard error
obs_rm <- data.frame("object"=unique(detections_rm$object), 
                     "Region.Label"=1, 
                     "Sample.Label"=1)
abundance_rm <- dht(model=detection_func_rm, region.table=region, 
                    sample.table=samples, obs.table=obs_rm)
Nhat_rm <- abundance_rm$individuals$N$Estimate
Nhat.se_rm <- abundance_rm$individuals$N$se



# § Simulating Multiple Surveys for varying Evasive Movement -------------------

set.seed(2000)

# Define range of responsive movement parameter
ρ_values <- seq(from=0, to=20, by=2)

# Define the number of surveys per evasive movement value
iterations <- 200

# Create empty vectors to store the means
mean_chisq.p_io <- numeric(length(ρ_values))
mean_Nhat_rm_io <- numeric(length(ρ_values))
mean_Nhat.se_rm_io <- numeric(length(ρ_values))
mean_p.0_rm_io <- numeric(length(ρ_values))
mean_p.0.se_rm_io <- numeric(length(ρ_values))

# Loop across evasive movement values
for (i in seq_along(ρ_values)) {
  # Set ρ value
  ρ <- ρ_values[i]
  
  # Create empty vectors to store the values for each iteration
  iteration_chisq.p <- numeric(iterations)
  iteration_Nhat <- numeric(iterations)
  iteration_Nhat.se <- numeric(iterations)
  iteration_p.0 <- numeric(iterations)
  iteration_p.0.se <- numeric(iterations)
  
  # Loop across surveys
  for (j in 1:iterations) {
    
    # Simulate detections in each survey as above
    x_rm <- runif(n=n, min=0, max=L)
    y <- runif(n=n, min=-w, max=w)
    indices_prop <- sample(seq_along(y), size=round(prop*length(y)))
    y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                      y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
    y_rm <- y
    y_rm[indices_prop] <- y_shift[indices_prop]
    probs_rm1 <- p_line*exp(-y_rm^2/(2*σ1^2))
    detections_rm1 <- rbinom(n=n, size=1, prob=probs_rm1)
    probs_rm2 <- exp(-y_rm^2/(2*σ2^2))
    detections_rm2 <- rbinom(n=n, size=1, prob=probs_rm2)
    detection_data_rm <- data.frame(
                            "object"=rep(1:n, each=2), 
                            "observer"=as.factor(rep(2:1, times=n)), 
                            "distance"=rep(abs(y_rm), each=2), 
                            "p_detection"=c(rbind(probs_rm1, probs_rm2)), 
                            "detected"=c(rbind(detections_rm1, detections_rm2)))
    
    # Remove animals not detected by either observer from dataframe
    zeros_rm <- detection_data_rm$detected==0
    frequency_rm <- table(detection_data_rm$object[zeros_rm])
    non_detections_rm <- as.numeric(names(frequency_rm[frequency_rm>1]))
    detections_rm <- detection_data_rm[!(detection_data_rm$object %in% 
                                           non_detections_rm),]
    
    # Fit detection function
    detection_func_rm <- ddf(method='io', 
                             mrmodel=~glm(link='logit', 
                                          formula=~distance + observer),
                             dsmodel=~cds(key='hn'), 
                             data=detections_rm, meta.data=list(width=w))
    
    # Produce goodness of fit statistics
    gof <- ddf.gof(detection_func_rm, qq=FALSE)
    chisq.p <- gof$chisquare$chi1$p
    
    # Extract estimate of primary detection probability on transect line
    summary_rm <- summary(detection_func_rm)
    p.0 <- summary_rm$mr.summary$average.p0.1
    p.0.se <- summary_rm$mr.summary$average.p0.1.se
    
    # Obtain abundance estimate
    obs_rm <- data.frame("object"=unique(detections_rm$object),
                         "Region.Label"=1, 
                         "Sample.Label"=1)
    abundance <- dht(model=detection_func_rm, region.table=region, 
                     sample.table=samples, obs.table=obs_rm)
    Nhat <- abundance$individuals$N$Estimate
    Nhat.se <- abundance$individuals$N$se
    
    # Store the values for each iteration
    iteration_chisq.p[j] <- chisq.p
    iteration_Nhat[j] <- Nhat
    iteration_Nhat.se[j] <- Nhat.se
    iteration_p.0[j] <- p.0
    iteration_p.0.se[j] <- p.0.se
  }
  # Calculate mean statistics for each evasive movement value
  mean_chisq.p_io[i] <- mean(iteration_chisq.p)
  mean_Nhat_rm_io[i] <- mean(iteration_Nhat)
  mean_Nhat.se_rm_io[i] <- mean(iteration_Nhat.se)
  mean_p.0_rm_io[i] <- mean(iteration_p.0)
  mean_p.0.se_rm_io[i] <- mean(iteration_p.0.se)
}

# Print mean values for each evasive movement value
for (i in seq_along(ρ_values)) {
  cat("ρ =", ρ_values[i], "\n")
  cat("Mean χ^2 p-value:", mean_chisq.p_io[i], "\n")
  cat("Mean Nhat:", mean_Nhat_rm_io[i], "\n")
  cat("Mean Nhat.se:", mean_Nhat.se_rm_io[i], "\n")
  cat("Mean p.0:", mean_p.0_rm_io[i], "\n")
  cat("Mean p.0.se:", mean_p.0.se_rm_io[i], "\n")
  cat("------------------------\n")
}

# Plot abundance results on numerical scale
par(mfrow=c(1,1))
plot(ρ_values, mean_Nhat_rm_io, pch=16, xlab="ρ", ylab=TeX(r'($\hat{N}$)'), 
     ylim=c(min(mean_Nhat_rm_io - mean_Nhat.se_rm_io),
            ifelse(max(mean_Nhat_rm_io + mean_Nhat.se_rm_io)>n, 
                   max(mean_Nhat_rm_io + mean_Nhat.se_rm_io), n)),
     main="Evasive Movement Effect on 
     Estimated Abundance (Independence Configuration)")
abline(n, 0, lty=2, col="red")
abline(lm(mean_Nhat_rm_io ~ ρ_values), col="blue")
abline(mean_Nhat_rm_io[1], 0, lty=2)
arrows(ρ_values, mean_Nhat_rm_io - mean_Nhat.se_rm_io, ρ_values, 
       mean_Nhat_rm_io + mean_Nhat.se_rm_io, angle=90, code=3, length=0.05)

# Plot abundance results on bias scale
Nhat_bias_io <- 100*(mean_Nhat_rm_io - n)/n
plot(ρ_values, Nhat_bias_io, pch=16, xlab="ρ", 
     ylab=TeX(r'($\%$ Bias in $\hat{N}$)'), 
     ylim=c((100*(min(mean_Nhat_rm_io - mean_Nhat.se_rm_io) - n)/n), 
            ifelse((100*(max(mean_Nhat_rm_io + mean_Nhat.se_rm_io) - n)/n)<0, 0,
                   (100*(max(mean_Nhat_rm_io + mean_Nhat.se_rm_io) - n)/n))),
     main="Evasive Movement Effect on 
     Estimated Abundance (Independence Configuration)")
abline(0, 0, lty=2, col="red")
abline(lm(Nhat_bias_io~ρ_values), col="blue")
abline(Nhat_bias_io[1], 0, lty=2)
arrows(ρ_values, 100*(mean_Nhat_rm_io - mean_Nhat.se_rm_io - n)/n, 
       ρ_values, 100*(mean_Nhat_rm_io + mean_Nhat.se_rm_io - n)/n, 
       angle = 90, code = 3, length = 0.05)

# Plot p(0) results
plot(ρ_values, mean_p.0_rm_io, pch=16, xlab="ρ", ylab=TeX(r'($\hat{p}_1(0)$)'), 
     ylim=c(p_line, max(mean_p.0_rm_io + mean_p.0.se_rm_io)),
     main="Evasive Movement Effect\non Estimated Primary p(0) 
     (Independence Configuration)")
arrows(ρ_values, mean_p.0_rm_io - mean_p.0.se_rm_io, ρ_values, 
       mean_p.0_rm_io + mean_p.0.se_rm_io, angle = 90, code = 3, length = 0.05)
abline(p_line, 0, lty=2, col="red")

# Plot p(0) results on bias scale
plot(ρ_values, 100*((mean_p.0_rm_io - p_line)/p_line), pch=16, xlab="ρ", 
     ylab=TeX(r'($\%$ Bias in $\hat{p}_1(0)$)'), 
     ylim=c(0, max(100*(mean_p.0_rm_io + mean_p.0.se_rm_io - p_line)/p_line)),
     main="Evasive Movement Effect\non Estimated Primary p(0)
     (Independence Configuration)",
     arrows(ρ_values, 100*(mean_p.0_rm_io - mean_p.0.se_rm_io - p_line)/p_line, 
            ρ_values, 100*(mean_p.0_rm_io + mean_p.0.se_rm_io - p_line)/p_line, 
            angle=90, code=3, length=0.05))
abline(0, 0, lty=2, col="red")

# Quantification of relationship between abundance and evasive movement
fit_io <- lm(mean_Nhat_rm_io ~ ρ_values)
summary(fit_io)
coefs_io <- fit_io$coefficients
a_i <- coefs_io[[1]]
b_i <- coefs_io[[2]]



# TRIAL CONFIGURATION

# § Simulating Uniform Distribution for Double Observers -----------------------

# Single survey

set.seed(2000)

# Randomly generate uniformly-distributed animals
x <- runif(n=n, min=0, max=L)
y <- runif(n=n, min=-w, max=w)

# Coördinates of randomly-generated animals
observations_unif <- matrix(c(x, y), nrow=2, ncol=n, byrow=TRUE, 
                            dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate boat detection probabilities based on half-normal function
probs_unif1 <- p_line*exp(-y^2/(2*σ1^2))
# Generate boat detections from Bernoulli trials
detections_unif1 <- rbinom(n=n, size=1, prob=probs_unif1)
# Detected animals by boat
data_unif1 <- y[which(detections_unif1==1)]

# Calculate drone detection probabilities based on step function
d <- 28   # Half-width of drone view
probs_unif2 <- ifelse(abs(y)<=d, 0.99, 0)
# Generate drone detections
detections_unif2 <- rbinom(n=n, size=1, prob=probs_unif2)
# Detected animals by drone
data_unif2 <- y[which(detections_unif2==1)]

# Plot uniformly-distributed animals with boat-detected animals coloured
par(mfrow=c(1,2))
plot(x, y, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="UD Boat Detections")
abline(0, 0, lty=3, col="red")
points(x[detections_unif1==1], y[detections_unif1==1], col="blue")

# Plot uniformly-distributed animals with drone-detected objects coloured
plot(x, y, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="UD Drone Detections")
abline(0, 0, lty=3, col="red")
points(x[detections_unif2==1], y[detections_unif2==1], col="red")

# Plot animals with those detected by both observers coloured 
par(mfrow=c(1,1))
plot(x, y, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="Uniform Distribution Duplicates")
abline(0, 0, lty=3, col="red")
points(x[detections_unif1==1 & detections_unif2==1], 
       y[detections_unif1==1 & detections_unif2==1], col="purple")

# Data from uniform distribution
detection_data_unif <- data.frame(
                        "object"=rep(1:n, each=2), 
                        "observer"=as.factor(rep(1:2, times=n)), 
                        "distance"=rep(abs(y), each=2),
                        "p_detection"=c(rbind(probs_unif1, probs_unif2)),
                        "detected"=c(rbind(detections_unif1, detections_unif2)))

# Remove animals not detected by either observer from dataframe
zeros_unif <- detection_data_unif$detected==0
frequency_unif <- table(detection_data_unif$object[zeros_unif])
non_detections_unif <- as.numeric(names(frequency_unif[frequency_unif>1]))
detections_unif <- detection_data_unif[!(detection_data_unif$object %in% 
                                           non_detections_unif),]

# Fit detection function under trial configuration and point independence
detection_func_unif <- ddf(method='trial', 
                           mrmodel=~glm(link='logit', formula=~distance),
                           dsmodel=~cds(key='hn'), data=detections_unif, 
                           meta.data=list(width=w))

# Plot histograms of double-observer detections
detection_tables_unif <- det.tables(detection_func_unif)
print(detection_tables_unif)
plot(detection_tables_unif, col2="blue")

# Plot unconditional and conditional detection functions
plot(detection_func_unif, xlab="Perpendicular Distance (m)",
     col="black", angle=60, density=25)

# Produce goodness of fit statistics and qq plot
gof_unif <- ddf.gof(detection_func_unif, 
                    main="Goodness of Fit for Uniform Animal Distribution
                    (Trial)")
chisq.p_unif <- gof_unif$chisquare$chi1$p

# Extract estimate of primary detection probability on transect line
summary_unif <- summary(detection_func_unif)
p.0_unif <- summary_unif$mr.summary$average.p0.1
p.0.se_unif <- summary_unif$mr.summary$average.p0.1.se

# Obtain estimates of abundance and its standard error
obs_unif <- data.frame("object"=unique(detections_unif$object),
                       "Region.Label"=1, 
                       "Sample.Label"=1)
abundance_unif <- dht(model=detection_func_unif, region.table=region, 
                      sample.table=samples, obs.table=obs_unif)
Nhat <- abundance_unif$individuals$N$Estimate
Nhat.se <- abundance_unif$individuals$N$se



# § Simulating Evasive Movement for Double Observers ---------------------------

# Single survey

set.seed(2000)

# Define responsive movement parameter
ρ <- 20    # Example value

# Randomly generate animals showing evasive movement
x_rm <- runif(n=n, min=0, max=L)
y <- runif(n=n, min=-w, max=w)
indices_prop <- sample(seq_along(y), size=round(prop*length(y)))
y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                        y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
y_rm <- y
y_rm[indices_prop] <- y_shift[indices_prop]

# Coördinates of randomly-generated animals
obs_rm <- matrix(c(x_rm, y_rm), nrow=2, ncol=n, byrow=TRUE, 
                 dimnames=list(c('x', 'y'), seq(1:n)))

# Calculate boat detection probabilities based on half-normal function
probs_rm1 <- p_line*exp(-y_rm^2/(2*σ1^2))
# Generate boat detections from Bernoulli trials
detections_rm1 <- rbinom(n=n, size=1, prob=probs_rm1)

# Calculate drone detection probabilities based on step function
probs_rm2 <- ifelse(abs(y_rm)<=d, 0.99, 0)
# Generate drone detections 
detections_rm2 <- rbinom(n=n, size=1, prob=probs_rm2)

# Plot uniformly-distributed animals with boat-detected animals coloured
par(mfrow=c(1,2))
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="EM Boat Detections")
abline(0, 0, lty=3, col="red")
points(x_rm[detections_rm1==1], y_rm[detections_rm1==1], col="blue")

# Plot evasive movement animals with drone-detected animals coloured
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="EM Drone Detections")
abline(0, 0, lty=3, col="red")
points(x_rm[detections_rm2==1], y_rm[detections_rm2==1], col="red")

# Plot animals with those detected by both observers coloured 
par(mfrow=c(1,1))
plot(x_rm, y_rm, xlab="Distance along Transect (m)", 
     ylab="Perpendicular distance (m)", main="EM Duplicates")
abline(0, 0, lty=3, col="red")
points(x_rm[detections_rm1==1 & detections_rm2==1], 
       y_rm[detections_rm1==1 & detections_rm2==1], col="purple")

# Data from uniform distribution
detection_data_rm <- data.frame(
                            "object"=rep(1:n, each=2), 
                            "observer"=as.factor(rep(1:2, times=n)), 
                            "distance"=rep(abs(y_rm), each=2),
                            "p_detection"=c(rbind(probs_rm1, probs_rm2)),
                            "detected"=c(rbind(detections_rm1, detections_rm2)))

# Remove animals not detected by either observer from dataframe
zeros_rm <- detection_data_rm$detected==0
frequency_rm <- table(detection_data_rm$object[zeros_rm])
non_detections_rm <- as.numeric(names(frequency_rm[frequency_rm>1]))
detections_rm <- detection_data_rm[!(detection_data_rm$object %in% 
                                       non_detections_rm),]

# Fit detection function under trial configuration and point independence
detection_func_rm <- ddf(method='trial', 
                         mrmodel=~glm(link='logit', formula=~distance),
                         dsmodel=~cds(key='hn'), data=detections_rm, 
                         meta.data=list(width=w))

# Plot histograms of double-observer detections
detection_tables_rm <- det.tables(detection_func_rm)
print(detection_tables_rm)
par(mfrow=c(1,1))
plot(detection_tables_rm, col2="blue")

# Plot unconditional and conditional detection functions
plot(detection_func_rm, xlab="Perpendicular Distance (m)",
     col="black", angle=60, density=25)

# Produce goodness of fit statistics and qq plot
gof_rm <- ddf.gof(detection_func_rm, 
                  main="Goodness of Fit for Evasive Movement (Trial)")
chisq.p_rm <- gof_rm$chisquare$chi1$p

# Extract estimate of primary detection probability on transect line
summary_rm <- summary(detection_func_rm)
p.0_rm <- summary_rm$mr.summary$average.p0.1
p.0.se_rm <- summary_rm$mr.summary$average.p0.1.se

# Obtain estimates of abundance and its standard error
obs_rm <- data.frame("object"=unique(detections_rm$object), 
                     "Region.Label"=1, 
                     "Sample.Label"=1)
abundance_rm <- dht(model=detection_func_rm, region.table=region, 
                    sample.table=samples, obs.table=obs_rm)
Nhat_rm <- abundance_rm$individuals$N$Estimate
Nhat.se_rm <- abundance_rm$individuals$N$se


# § Simulating Multiple Surveys for varying Evasive Movement -------------------

set.seed(2000)

# Define range of responsive movement parameter
ρ_values <- seq(from=0, to=20, by=2)

# Define the number of surveys per ρ value
iterations <- 200

# Create empty vectors to store the means
mean_chisq.p_trial <- numeric(length(ρ_values))
mean_Nhat_rm_trial <- numeric(length(ρ_values))
mean_Nhat.se_rm_trial <- numeric(length(ρ_values))
mean_p.0_rm_trial <- numeric(length(ρ_values))
mean_p.0.se_rm_trial <- numeric(length(ρ_values))

# Loop across evasive movement values
for (i in seq_along(ρ_values)) {
  # Set ρ value
  ρ <- ρ_values[i]
  
  # Create empty vectors to store the values for each iteration
  iteration_chisq.p <- numeric(iterations)
  iteration_Nhat <- numeric(iterations)
  iteration_Nhat.se <- numeric(iterations)
  iteration_p.0 <- numeric(iterations)
  iteration_p.0.se <- numeric(iterations)
  
  # Loop across surveys
  for (j in 1:iterations) {
    
    # Simulate detections in each survey as above
    x_rm <- runif(n=n, min=0, max=L)
    y <- runif(n=n, min=-w, max=w)
    indices_prop <- sample(seq_along(y), size=round(prop*length(y)))
    y_shift <- ifelse(y>=0, y + ifelse(y<z, ρ/z*(z-y), 0), 
                      y - ifelse(abs(y)<z, ρ/z*(z-abs(y)), 0))
    y_rm <- y
    y_rm[indices_prop] <- y_shift[indices_prop]
    probs_rm1 <- p_line*exp(-y_rm^2/(2*σ1^2))
    detections_rm1 <- rbinom(n=n, size=1, prob=probs_rm1)
    probs_rm2 <- ifelse(abs(y_rm)<=d, 0.99, 0)
    detections_rm2 <- rbinom(n=n, size=1, prob=probs_rm2)
    detection_data_rm <- data.frame(
      "object"=rep(1:n, each=2), 
      "observer"=as.factor(rep(1:2, times=n)), 
      "distance"=rep(abs(y_rm), each=2), 
      "p_detection"=c(rbind(probs_rm1, probs_rm2)), 
      "detected"=c(rbind(detections_rm1, detections_rm2)))
    
    # Remove animals not detected by either observer from dataframe
    zeros_rm <- detection_data_rm$detected==0
    frequency_rm <- table(detection_data_rm$object[zeros_rm])
    non_detections_rm <- as.numeric(names(frequency_rm[frequency_rm>1]))
    detections_rm <- detection_data_rm[!(detection_data_rm$object %in% 
                                           non_detections_rm),]
    
    # Fit detection function
    detection_func_rm <- ddf(method='trial', 
                             mrmodel=~glm(link='logit', formula=~distance),
                             dsmodel=~cds(key='hn'), 
                             data=detections_rm, meta.data=list(width=w),
                             control=list(refit=TRUE))
    
    # Perform the second fit without distance covariate if coefficient >0
    if(isTRUE(detection_func_rm$mr$mr$coefficients[2]>0)){
      detection_func_rm <- ddf(method='trial', 
                               mrmodel=~glm(link='logit', formula=~1),
                               dsmodel=~cds(key='hn'), 
                               data=detections_rm, meta.data=list(width=w),
                               control=list(refit=TRUE))
    }
    
    # Produce goodness of fit statistics
    gof <- ddf.gof(detection_func_rm, qq=FALSE)
    chisq.p <- gof$chisquare$chi1$p
    
    # Extract estimate of primary detection probability at distance 0
    summary_rm <- summary(detection_func_rm)
    p.0 <- summary_rm$mr.summary$average.p0.1
    p.0.se <- summary_rm$mr.summary$average.p0.1.se
    
    # Obtain abundance estimate
    region <- data.frame("Region.Label"=1, "Area"=2*w)
    samples <- data.frame("Sample.Label"=1, "Region.Label"=1, "Effort"=1)
    obs_rm <- data.frame("object"=unique(detections_rm$object),"Region.Label"=1, 
                         "Sample.Label"=1)
    abundance <- dht(model=detection_func_rm, region.table=region, 
                     sample.table=samples, obs.table=obs_rm)
    Nhat <- abundance$individuals$N$Estimate
    Nhat.se <- abundance$individuals$N$se
    
    # Store the values for each iteration
    iteration_chisq.p[j] <- chisq.p
    iteration_Nhat[j] <- Nhat
    iteration_Nhat.se[j] <- Nhat.se
    iteration_p.0[j] <- p.0
    iteration_p.0.se[j] <- p.0.se
  }
  
  # Calculate mean statistics for each evasive movement value
  mean_chisq.p_trial[i] <- mean(iteration_chisq.p)
  mean_Nhat_rm_trial[i] <- mean(iteration_Nhat)
  mean_Nhat.se_rm_trial[i] <-mean(iteration_Nhat.se)
  mean_p.0_rm_trial[i] <- mean(iteration_p.0[iteration_Nhat!=0])
  mean_p.0.se_rm_trial[i] <- mean(iteration_p.0.se[iteration_Nhat!=0])
}

# Print the mean values for each evasive movement value
for (i in seq_along(ρ_values)) {
  cat("ρ =", ρ_values[i], "\n")
  cat("Mean χ^2 p-value:", mean_chisq.p_trial[i], "\n")
  cat("Mean Nhat:", mean_Nhat_rm_trial[i], "\n")
  cat("Mean Nhat.se:", mean_Nhat.se_rm_trial[i], "\n")
  cat("Mean p.0:", mean_p.0_rm_trial[i], "\n")
  cat("Mean p.0.se:", mean_p.0.se_rm_trial[i], "\n")
  cat("------------------------\n")
}

# Plot abundance results on numerical scale
plot(ρ_values, mean_Nhat_rm_trial, pch=16, xlab="ρ", ylab=TeX(r'($\hat{N}$)'), 
     ylim=c(min(mean_Nhat_rm_trial - mean_Nhat.se_rm_trial),
            ifelse(max(mean_Nhat_rm_trial + mean_Nhat.se_rm_trial)<n, n, 
                   max(mean_Nhat_rm_trial + mean_Nhat.se_rm_trial))),
     main="Evasive Movement Effect on 
     Estimated Abundance (Trial Configuration)")
arrows(ρ_values, mean_Nhat_rm_trial - mean_Nhat.se_rm_trial, ρ_values, 
       mean_Nhat_rm_trial + mean_Nhat.se_rm_trial, angle=90, code=3,length=0.05)
abline(n, 0, lty=2, col="red")
abline(lm(mean_Nhat_rm_trial ~ ρ_values), col="blue")
abline(mean_Nhat_rm_trial[1], 0, lty=2)

# Plot abundance results on bias scale
Nhat_bias_trial <- 100*((mean_Nhat_rm_trial - n)/n)
plot(ρ_values, Nhat_bias_trial, pch=16, xlab="ρ", 
     ylab=TeX(r'($\%$ Bias in $\hat{N}$)'), 
     ylim=c((100*(min(mean_Nhat_rm_trial - mean_Nhat.se_rm_trial) - n)/n), 
            ifelse((100*(max(mean_Nhat_rm_trial+mean_Nhat.se_rm_trial)-n)/n)<0,0,
                   (100*(max(mean_Nhat_rm_trial+mean_Nhat.se_rm_trial)-n)/n))),
     main="Evasive Movement Effect on 
     Estimated Abundance (Trial Configuration)")
arrows(ρ_values, 100*(mean_Nhat_rm_trial - mean_Nhat.se_rm_trial - n)/n, 
       ρ_values, 100*(mean_Nhat_rm_trial + mean_Nhat.se_rm_trial - n)/n, 
       angle=90, code=3, length=0.05)
abline(0, 0, lty=2, col="red")
abline(lm(Nhat_bias_trial ~ ρ_values), col="blue")
abline(Nhat_bias_trial[1], 0, lty=2)

# Plot p(0) results on numerical scale
plot(ρ_values, mean_p.0_rm_trial, pch=16, xlab="ρ", 
     ylab=TeX(r'($\hat{p}_1(0)$)'), 
     ylim=c(min(mean_p.0_rm_trial - mean_p.0.se_rm_trial), 
            max(mean_p.0_rm_trial + mean_p.0.se_rm_trial)),
     main="Evasive Movement Effect\non Estimated Primary p(0) 
     (Trial Configuration)")
arrows(ρ_values, mean_p.0_rm_trial - mean_p.0.se_rm_trial, 
       ρ_values, mean_p.0_rm_trial + mean_p.0.se_rm_trial, 
       angle=90, code=3, length=0.05)
abline(p_line, 0, lty=2, col="red")

# Plot p(0) results on bias scale
plot(ρ_values, 100*((mean_p.0_rm_trial - p_line)/p_line), pch=16, xlab="ρ", 
     ylab=TeX(r'($\%$ Bias in $\hat{p}_1(0)$)'), 
     ylim=c(min(100*(mean_p.0_rm_trial - mean_p.0.se_rm_trial - p_line)/p_line), 
           max(100*(mean_p.0_rm_trial + mean_p.0.se_rm_trial - p_line)/p_line)),
     main="Evasive Movement Effect\non Estimated Primary p(0) 
     (Trial Configuration)",
arrows(ρ_values, 100*(mean_p.0_rm_trial - mean_p.0.se_rm_trial - p_line)/p_line, 
       ρ_values, 100*(mean_p.0_rm_trial + mean_p.0.se_rm_trial - p_line)/p_line, 
       angle=90, code=3, length=0.05))
abline(0, 0, lty=2, col="red")

# Quantification of relationship between abundance and evasive movement
fit_trial <- lm(mean_Nhat_rm_trial ~ ρ_values)
summary(fit_trial)
coefs_trial <- fit_trial$coefficients
a_t <- coefs_trial[[1]]
b_t <- coefs_trial[[2]]



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------                                                   ----------------
#-----------   DETECTION FUNCTION SCALE PARAMETER INVESTIGATION  ---------------
#------------                                                   ----------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

set.seed(2000)

# Number of simulated surveys
iterations <- 500

# Range of detection function scale parameter values
σ_values <- seq(from=0.1*w, to=0.5*w, by=0.05*w)

# Create empty vectors to store results
mean_Nhat_sigma <- numeric(length(σ_values))
mean_Nhat.se_sigma <- numeric(length(σ_values))

# Loop across the scale parameter values
for (i in 1:length(σ_values)) {
  
  # Create an empty vector to store the results for the current sigma
  iteration_Nhats <- numeric(iterations)
  iteration_Nhats.se <- numeric(iterations)
  
  # Set σ value
  σ <- σ_values[i]
  
  # Loop across surveys
  for (j in 1:iterations) {
    
    # Simulate each survey
    x <- runif(n=n, min=0, max=L)
    y <- runif(n=n, min=-w, max=w)
    probs_unif <- exp(-y^2/(2*σ^2))
    detections_unif <- rbinom(n=n, size=1, prob=probs_unif)
    detection_data_unif <- data.frame(
      "object" = seq(1:n),
      "distance" = abs(y),
      "p_detection" = probs_unif,
      "detected" = detections_unif
    )
    
    # Fit the detection functions with distance sampling method
    ds_unif <- ddf(method='ds', dsmodel=~cds(key='hn'), 
                   data=detection_data_unif, meta.data=list(width=w))
    
    # Obtain estimates of abundance and its standard error
    obs_unif <- data.frame("object"=unique(detection_data_unif$object),
                           "Region.Label"=1, 
                           "Sample.Label"=1)
    abundance_unif <- dht(model=ds_unif, region.table=region,
                          sample.table=samples, obs.table=obs_unif)
    Nhat_single <- abundance_unif$individuals$N$Estimate
    Nhat.se_single <- abundance_unif$individuals$N$se
    
    # Store estimated abundance
    iteration_Nhats[j] <- Nhat_single
    iteration_Nhats.se[j] <- Nhat.se_single
  }
  
  # Store the mean and standard error for the current σ value
  mean_Nhat_sigma[i] <- mean(iteration_Nhats)
  mean_Nhat.se_sigma[i] <- mean(iteration_Nhats.se)
}

# Mean and standard deviation of abundance estimates from each σ value
for (i in 1:length(σ_values)) {
  σ <- σ_values[i]
  cat("Results for σ =", σ, "\n")
  print(paste("Mean:", mean_Nhat_sigma[i]))
  print(paste("Standard Error:", mean_Nhat.se_sigma[i]))
}

# Plot results on numerical scale
plot(σ_values, mean_Nhat_sigma, pch=16, xlab="σ", 
     ylab=TeX(r'($\hat{N}$)'), 
     ylim=c(min(mean_Nhat_sigma - mean_Nhat.se_sigma),
            max(mean_Nhat_sigma + mean_Nhat.se_sigma)),
     main="Detection Function Scale Effect on 
     Estimated Abundance by Single Observer")
abline(n, 0, lty=2, col="red")
arrows(σ_values, mean_Nhat_sigma - mean_Nhat.se_sigma, σ_values, 
       mean_Nhat_sigma + mean_Nhat.se_sigma, angle=90, code=3, length=0.05)

# Plot results on bias scale
plot(σ_values, 100*((mean_Nhat_sigma-n)/n), pch=16, xlab="σ", 
     ylab=TeX(r'($\%$ Bias in $\hat{N}$)'), 
     ylim=c((100*(min(mean_Nhat_singma - mean_Nhat.se_sigma) - n)/n), 
            (100*(max(mean_Nhat_sigma + mean_Nhat.se_sigma) - n)/n)),
     main="Detection Function Scale Effect on 
     Estimated Abundance by Single Observer")
abline(0, 0, lty=2, col="red")
arrows(σ_values, 100*(mean_Nhat_sigma - mean_Nhat.se_sigma - n)/n, 
       σ_values, 100*(mean_Nhat_sigma + mean_Nhat.se_sigma - n)/n, 
       angle=90, code=3, length=0.05)

# Signal that code has finished running
beepr::beep(sound=8)
