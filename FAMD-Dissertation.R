# Load required libraries
library(ggplot2)   # For data visualization
library("survival")
library("dplyr")
library("htmltools")
library("FactoMineR")
library("factoextra")


#Reading in the dataset
df <- read.csv('heart_failure_clinical_records_dataset.csv')
str(df)


# Previewing the dataset
head(df)     # View the first few rows of the dataset
str(df)      # Display the structure of the dataset
summary(df)  # Show summary statistics of the dataset

# Exploratory data analysis

# Data visualization
# Creating a histogram of a numeric variable
ggplot(df, aes(x = age)) +
  geom_histogram() +
  labs(title = "Histogram of Numeric Variable")

# Performing data summarisation or aggregation
# Calculating the mean of a numeric variable by a categorical variable
df_summary <- df %>%
  group_by(DEATH_EVENT) %>%
  summarise(mean_age= mean(age))

# Viewing the summarized data
head(df_summary)


str(df)
#Factor Analysis
df2 <- df
df4 <- df

df3 <- sapply(df,as.factor)
df3 <- as.data.frame(df3)
summary(df3)
df3$age <- as.numeric(df3$age)
df3$creatinine_phosphokinase <- as.numeric(df3$creatinine_phosphokinase)
df3$ejection_fraction <- as.numeric(df3$ejection_fraction)
df3$platelets <- as.numeric(df3$platelets)
df3$serum_creatinine <- as.numeric(df3$serum_creatinine)
df3$serum_sodium <- as.numeric(df3$serum_sodium)
df3$time <- as.numeric(df3$time)
#Convert other reqd columns to numerics 

# Dataframe has columns containing 0s and 1s
# Specify the column names to convert
columns_to_convert <- c("smoking", "anaemia", "diabetes", "high_blood_pressure", "sex", "DEATH_EVENT")

# Iterate over the columns and convert them to categorical variables
for (col in columns_to_convert) {
  if (any(df2[[col]] %in% c(0, 1))) {
    df2[[col]] <- factor(df2[[col]], levels = c(0, 1), labels = c("No", "Yes"))
  }
}

for (i in 1:ncol(df4)) {
  col_values <- unique(df4[, i])
  
  # Check if the column contains only 0 and 1 values
  if (all(col_values %in% c(0, 1))) {
    col_name <- colnames(df4)[i]
    df4[, i] <- ifelse(df4[, i] == 0, paste0("No", i), paste0("Yes", i))
    colnames(df4)[i] <- col_name
  }
}
# View the updated dataframe
head(df2)
head(df3)
head(df4)

str(df2)
#removing target variable before famd
df2 <- df2[, -which(names(df2) == "DEATH_EVENT")]
df4 <- df4[, -which(names(df4) == "DEATH_EVENT")]

#FAMD function
res.famd <- FAMD(df4, ncp=12, graph = FALSE)
print(res.famd)

res.famd$eig

library(corrplot)
corrplot(res.famd$var$coord)

#Eigen values
eig.val <- get_eigenvalue(res.famd)
eig.val

#Scree plot
fviz_screeplot(res.famd)

# The coordinates, the cos2 and the contribution of all variables
var <- get_famd_var(res.famd)
var

# Coordinates of variables
var$coord
# Cos2: quality of representation on the factore map
var$cos2
# Contributions to the  dimensions
var$contrib 

#To show the correlation between variables (both quantitative and qualitative variables) and the principal dimensions
# Plot of variables
fviz_famd_var(res.famd, axes=c(1,2), repel = TRUE)

fviz_famd_var(res.famd, axes=c(2,3), repel = TRUE)

fviz_famd_var(res.famd, axes=c(3,4), repel = TRUE)

fviz_famd_var(res.famd, axes=c(4,5), repel = TRUE)

fviz_famd_var(res.famd, axes=c(5,6), repel = TRUE)

#look betwwwn other dimensions like 2 and 3 etc
# Contribution to the first dimension
fviz_contrib(res.famd, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2)

fviz_contrib(res.famd, "var", axes = 3)
#Quantitative variables

quanti.var <- get_famd_var(res.famd, "quanti.var")
quanti.var 

#Visualising / Factor map
fviz_famd_var(res.famd, "quanti.var", repel = TRUE, col.var = "red")

#To find the most contributing quantitative variables
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

# Color by cos2 values: quality of the factor map
fviz_famd_var(res.famd, "quanti.var", col.var = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
              repel = TRUE)

#Qualitative variable

row.names(df3) = make.names(row.names(df3),unique = TRUE)

quali.var <- get_famd_var(res.famd, "quali.var")
quali.var 

#Visualisation  - Error for this visualisation
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

ind <- get_famd_ind(res.famd)
ind

fviz_famd_ind(res.famd, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)


fviz_mfa_ind(res.famd, 
             habillage = df$DEATH_EVENT, # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 

# Above code Showing some error




















df1 <- df
# Convert time from days to months
df1$time_in_months <- df1$time / 30

# Round the converted values to the nearest whole number
df1$time_in_months <- round(df1$time_in_months)

# View the updated dataset
head(df1)


time_sur <- df$time
event <- df$DEATH_EVENT
time_surv <- df1$time_in_months

par(mfrow = c(1, 2))
# Creating a survival object
Ex_Surv <- Surv(time_sur, event)
Ex_survfit <- survfit(Ex_Surv ~ 1)

# Creating a Kaplan-Meier plot
plot(Ex_survfit, conf.int="none", mark.time=TRUE, xlab="Time (days)", ylab="Proportion Survival", main= "Kaplan-Meier plot")


# Creating a survival object
Ex_Surv1 <- Surv(time_surv, event)
Ex_survfit1 <- survfit(Ex_Surv1 ~ 1)

# Creating a Kaplan-Meier plot
plot(Ex_survfit1, conf.int="none", mark.time=TRUE, xlab="Time (months)", ylab="Proportion Survival", main= "Kaplan-Meier plot")




































