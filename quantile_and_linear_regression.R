################################################################################
#                                  Fabian SL Yii                               #
#                               fabian.yii@ed.ac.uk                            #
################################################################################

library(dplyr)
library(quantreg)
library(vctrs)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(sjPlot)
library(car)

## Clear workspace
rm(list=ls())

## Set working directory to parent directory
setwd("/Users/fabianyii/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Projects/UKB_cornea/")

## Read the cleaned tabular data in long format
# 137016 eyes (68508 RE & 68508 LE)
d_long  <- read.csv('data/cleaned_data_long_all.csv')

## CSV with IDs to be removed due to congenital malformations such as
## Marfan and Ehlers-Danlos syndromes
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=132568
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=132584
removeDF    <- read.csv("data/ids_to_remove.csv")
removeIDs   <- removeDF$id[removeDF$id %in% d_long$id]
d_long$congenitalMalformations <- FALSE
d_long[d_long$id %in% removeDF$id, ]$congenitalMalformations <- TRUE 


#############################################################################################
########################## Quality control & participant selection ##########################
#############################################################################################
## Include eyes with spherical refractive error (SER) data
# 135019 eyes (67607 RE & 67412 LE)
d_long <- d_long[!is.na(d_long$SER), ] 

## Include eyes with corneal hysteresis (CH) and corneal resistance factor (CRF) data
# 130722 eyes (65546 RE & 65176 LE)
d_long <- d_long[!(is.na(d_long$CH) | is.na(d_long$CRF)), ] 

## Include eyes with corneal radius (CR) data
# 124864 eyes (62475 RE & 62389 LE)
d_long <- d_long[!is.na(d_long$meanCornealRadius), ]

## Include eyes with Goldmann correlated intraocular pressure (IOP) data
# 124864 eyes (62475 RE & 62389 LE)
d_long <- d_long[!is.na(d_long$IOP), ]

## Remove bottom and top 0.5% of CH, CRF, CR, IOP
CH_out  <- quantile(d_long$CH, probs=c(0.005,0.995), na.rm=TRUE)
CRF_out <- quantile(d_long$CRF, probs=c(0.005,0.995), na.rm=TRUE)
CR_out  <- quantile(d_long$meanCornealRadius, probs=c(0.005,0.995), na.rm=TRUE)
IOP_out <- quantile(d_long$IOP, probs=c(0.005,0.995), na.rm=TRUE)
# 3538 eyes to be removed
remove_bool <- d_long$CH<CH_out[1] | d_long$CH>CH_out[2] | d_long$CRF<CRF_out[1] | d_long$CRF>CRF_out[2] | d_long$IOP<IOP_out[1] | d_long$IOP>IOP_out[2] | d_long$meanCornealRadius<CR_out[1]| d_long$meanCornealRadius>CR_out[2]
# 121326 eyes (60775 RE & 60551 LE)
d_long <- d_long[!remove_bool,]

## Include eyes with VA data
# 121001 eyes (60612 eyes & 60389 eyes)
d_long <- d_long[!is.na(d_long$VA), ] 

## Include eyes with VA > 0.00LogMAR 
# 69800 eyes (34587 RE & 35213 LE)
d_long <- d_long[d_long$VA <= 0, ] 

## Include participants with good systemic health 
# 48465 eyes (23976 RE & 24489 LE)
# Of which, 16390 eyes have hypertension;
# 3018 eyes have diabetes;
# 1208 eyes have myocardial infarction.
# Together, they represent 17964 (84.2%) out of 21335 removals 
d_long <- d_long[rowSums(d_long[,c(23:44,79)]) == 0, ] 

## Include participants with good ocular health 
# 47222 eyes (23372 RE & 23850 LE)
# Of which, 491 eyes have glaucoma;
# 431 eyes have chorioretinal diseases;
# 154 eyes have globe/scleral disorder including degenerative myopia and staphyloma.
# Together, they represent 1050 (84.5%) out of 1243 of removals. 
# Note 1: 96 eyes have strabismus
# Note 2: 51 eyes have some form of optic nerve disorder
# Note 3: 50 eyes have some form of corneal disorder
# Note 4: 8 eyes have nystagmus
d_long <- d_long[rowSums(d_long[,45:52]) == 0, ] 

## Include only myopes
# 16099 eyes (8005 RE & 8094 LE)
d_long$SER <- round(d_long$SER, 2)
d_long <- subset(d_long, SER <= -0.5)

## 4877 individuals have one eye only; create a new dataframe "d" containing only these IDs
RE <- subset(d_long, eye=="RE")
LE <- subset(d_long, eye=="LE")
oneEyeIDs <- c(setdiff(RE$id, LE$id), setdiff(LE$id, RE$id))
oneEyeIDs <- unique(oneEyeIDs)
d <- d_long[d_long$id %in% oneEyeIDs,]

## 5611 individuals have both eyes present (select one eye at random)
twoEyesIDs <- unique(d_long[!d_long$id %in% oneEyeIDs,]$id) # unique IDs with both eyes present
set.seed(101010)                                            # initialise a pseudorandom number generator
selectedEyes <- rbinom(length(twoEyesIDs), 1, 0.5)          # sample 0 or 1 at random for each unique ID
selectedEyes <- ifelse(selectedEyes==0, "RE", "LE")         # convert 0 to "RE" and 1 to "LE"
twoEyesData  <- d_long[d_long$id %in% twoEyesIDs,]          # extract rows corresponding to IDs with both eyes present
# For each unique ID extract the row corresponding to the eye randomly 
# selected above and append it to the new dataframe "d" created above
index <- 1                                 
for (id in twoEyesData$id){
  sub_d <- twoEyesData[twoEyesData$id==id & twoEyesData$eye==selectedEyes[index], ]
  d     <- rbind(d, sub_d)
  index <- index + 1 
}
# Remove empty rows; 10488 eyes of 10488 unique IDs finally included (5208 RE & 5280 LE)
d <- d[!is.na(d$id), ]


#############################################################################################
##################### Ordinary least squares multiple linear regression #####################
#############################################################################################
CH_lm  <- lm(SER ~ scale(CH) + scale(meanCornealRadius) + scale(IOP) + scale(age) + factor(sex), d)  # corneal hysteresis model
CRF_lm <- lm(SER ~ scale(CRF) + scale(meanCornealRadius) + scale(IOP) + scale(age) + factor(sex), d) # corneal resistance factor model

# Check variance inflation factor for each variable (all < 2)
vif(CH_lm) 
vif(CRF_lm)

# Very strong correlation between (pearson's r: 0.85), so the diagnostic plots are similar for both models
cor.test(d$CH, d$CRF)

############################# Residual analysis: diagnostic plots ###########################
my_theme <- theme_blank() + theme(plot.subtitle=element_text(size=10, hjust=1),
                                  panel.grid.major=element_line(color="gray94", linewidth=0.3),
                                  axis.text=element_text(size=7.5),
                                  axis.title=element_text(size=8),
                                  plot.subtitle=element_text(size=5)) 
## Normal Q-Q plot
# CH
p1 <- gg_qqplot(CH_lm, scale.factor=0.6) + 
      labs(title="", y="Observed residuals", x="Theoretical residuals", subtitle="Normal Q-Q plot") + 
      ylim(-20,5) +
      my_theme
# CRF
p2 <- gg_qqplot(CRF_lm, scale.factor=0.6) + 
      labs(title="", y="Observed residuals", x="Theoretical residuals", subtitle="Normal Q-Q plot") + 
      ylim(-20,5) +
      my_theme
      
## Residuals vs Fitted plot
# CH
p3 <- gg_resfitted(CH_lm, scale.factor=0.5) + 
      labs(title="", y="", x="Predicted SER (D)", subtitle="Residuals vs fitted plot") + 
      scale_x_continuous(breaks=seq(-1,-4.5,-0.5), labels=seq(-1,-4.5,-0.5) ) +
      ylim(-20, 5) + 
      my_theme
# CRF
p4 <- gg_resfitted(CRF_lm, scale.factor=0.5) + 
      labs(title="", y="", x="Predicted SER (D)", subtitle="Residuals vs fitted plot") + 
      scale_x_continuous(breaks=seq(-1,-4.5,-0.5), labels=seq(-1,-4.5,-0.5) ) +
      ylim(-20, 5) + 
      my_theme

# CH diagnostic plots
ggarrange(p1, p3, 
          hjust=-0.08,
          ncol=2, 
          nrow=1,
          font.label=list(size=11),
          labels=c("Ordinary least squares model:"))

# Save diagnostic plots
ggsave("figures/cornea_biomechanics_lm_residual.png", width=9, height=5.5, units="in", bg="white")

# ####   Combined (CH and CRF) diagnostic plots       ####
# # Very strong correlation between (pearson's r: 0.85), #
# # so the diagnostic plots are similar for both models. #
# ggarrange(p1, p3, p2, p4, hjust=-0.05,
#           ncol = 2, nrow = 2,
#           font.label = list(size=12),
#           labels = c("Corneal hysteresis (CH) OLS model:", "", "Corneal resistance factor (CRF) OLS model:", ""))


#############################################################################################
################################## Quantile regression ######################################
#############################################################################################
# 49 equally spaced SER quantiles ranging from 0.02 (towards myopia) to 0.98 (towards emmetropia)
taus   <- seq(0.02, 0.98, 0.02)
n_taus <- length(taus)

# CH model
CH_model <- rq(SER ~ scale(CH) + scale(meanCornealRadius) + scale(IOP) + scale(age) + factor(sex), 
               tau = taus, 
               data = d) 
# CRF model
CRF_model <- rq(SER ~ scale(CRF) + scale(meanCornealRadius) + scale(IOP) + scale(age) + factor(sex), 
                tau = taus, 
                data = d) 

## Internal function: extract refractive quantiles (taus), beta coefficients, 
## standard error & p values associated with each variable from the raw quantile
## regression output 
extract_results <- function(summaryObject, feature){
  
  # Variable names to be extracted
  variables <- c("Intercept", feature, "CR", "IOP", "Age", "Male")
  n_variables <- length(variables)
  
  # Create empty dataframe
  df <- data.frame("variables"=vec_rep_each(variables, n_taus),
                   "taus"=vec_rep(taus, n_variables),
                   "coef"=rep(NA, n_variables*n_taus),
                   "SE"=rep(NA, n_variables*n_taus),
                   "p_val"=rep(NA, n_variables*n_taus))
  
  # Start extracting results
  for(variable in variables){
    for(tau in taus){
      df_row <- (df$variables == variable) & (df$taus == tau)
      tau_ind <- match(tau, taus)
      variable_ind <- match(variable, variables)
      df[df_row,]$coef <- summaryObject[tau_ind][[1]]$coefficient[variable_ind, 1]
      df[df_row,]$SE <- summaryObject[tau_ind][[1]]$coefficient[variable_ind, 2]
      df[df_row,]$p_val <- summaryObject[tau_ind][[1]]$coefficient[variable_ind, 4]
    }
  }
  # Return dataframe
  return(df) }

# Extract pertinent data from the raw quantile regression output
CH_df  <- extract_results(summary(CH_model), "CH" )
CRF_df <- extract_results(summary(CRF_model), "CRF" )

################################### Coefficient plot ########################################
# Create a new dataframe "plot_coef" to save extracted parameters for plotting
plot_df <- rbind(CH_df, CRF_df)
plot_coef <- filter(plot_df, variables %in% c("CH", "CRF"))
plot_coef$intercept <- NA
plot_coef[which(plot_coef$variables=="CH"),]$intercept <- plot_df[1:49,]$coef
plot_coef[which(plot_coef$variables=="CRF"),]$intercept <- plot_df[295:343,]$coef
plot_coef$intercept <- round(plot_coef$intercept,2)

# Specify colour parameters
fg_col <- rgb(0.8,0.4,0, alpha=0.03)
point_cols <- rep("black", nrow(plot_coef))
point_cols[which(plot_coef$p_val<0.05 & plot_coef$variables == "CH")] <- "red"
point_cols[which(plot_coef$p_val<0.05 & plot_coef$variables == "CRF")] <- "lightblue"

# Add a new column to the dataframe to save dummy coefficients for plotting
plot_coef$dummy_coef <- ifelse(plot_coef$coef>0, plot_coef$coef+0.1, plot_coef$coef+0.1)

# Start plotting
ggplot(data=plot_coef, aes(x=taus, y=dummy_coef)) +
  geom_point(size=1.5, alpha=0.8, colour=point_cols) +
  geom_smooth(aes(group=variables, colour=variables), size=0.8, alpha=0.5, se=FALSE ) +
  geom_ribbon(aes(ymin=dummy_coef-SE*1.96, ymax=dummy_coef+SE*1.96, group=variables, fill=variables), alpha=0.08) +
  geom_hline(yintercept=0.1, size=0.5, colour="black") +
  theme_wsj(color="white") +
  labs(x="Refractive quantile", caption="Quantile regression") +
  scale_y_continuous(name="Standardised beta",
                     breaks=seq(0.1, 1.0, 0.05),
                     labels=sprintf("%.2f", seq(0.0, 0.9, 0.05)),
                     sec.axis=sec_axis(~.*17, 
                                       breaks=seq(0,-10,-1),
                                       labels=seq(0,-10,-1),
                                       name="SER (D)")) +
  geom_segment(aes(x=taus, xend=taus, y=0, yend=intercept/17), 
               colour=hcl.colors(n_taus*2, palette = "Heat 2"), 
               lwd=0.5, alpha=0.5) +
  geom_point(aes(y=intercept/17), colour=hcl.colors(n_taus*2, palette = "Heat 2") ) +
  annotate("text", label="p > 0.05 (black points)", 0.85, 0.4, size=4, colour="gray78") +
  annotate("text", label="Line of null effect", 0.08, 0.13, size=3, colour="gray48") +
  theme(legend.title=element_blank(), 
        plot.caption=element_text(size=9, hjust=-0.2),
        plot.background = element_rect(fill=fg_col[1]),
        axis.title=element_text(size=11),
        axis.title.x=element_text(margin=margin(t=10)),
        axis.title.y.left=element_text(hjust=0.75, margin=margin(r=10)),
        axis.title.y.right=element_text(hjust=0.83, margin=margin(l=15), colour=hcl.colors(1, palette="Heat 2") ),
        axis.text=element_text(size=9),
        axis.text.y.right=element_text(colour=hcl.colors(17, palette="Heat 2", rev=TRUE)[8:17]),
        panel.grid.major=element_blank(), 
        axis.ticks.y=element_line(color="black")) + 
  scale_x_continuous(breaks=seq(0,1, 0.1), limits=c(0, 1), labels= seq(0, 1, 0.1) ) +
  scale_colour_discrete(labels= c("Corneal hysteresis (CH)", "Corneal resistance factor (CRF)")) +
  scale_fill_discrete(labels= c("Corneal hysteresis (CH)", "Corneal resistance factor (CRF)")) 
ggsave("figures/cornea_biomechanics_qr.png", width=6, height=6.5, units="in", bg="white")

