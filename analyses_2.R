run_everything = FALSE
library(haven)
library(psych)
library(motley)
library(MplusAutomation)
library(tidySEM)
source("c:/git_repositories/miniMAC/factor_syntax.R")
source("c:/git_repositories/miniMAC/rfunctions.R")
source("parceling.r")
df <- read_dta("sem_update.dta")

names(df) <- tolower(names(df))
names(df) <- gsub("social_cohesion", "coh", names(df))
names(df) <- gsub("^dist_", "", names(df))
names(df) <- gsub("^noise_", "", names(df))
names(df) <- gsub("^phq9_", "phq", names(df))

scales_names <- unique(gsub("\\d+$", "", grep("\\d+$", names(df), value = TRUE)))

scales_list <- lapply(scales_names, function(x){
  grep(paste0(x, "\\d+$"), names(df), value = TRUE)
})
names(scales_list) <- scales_names
scales_list <- scales_list[sapply(scales_list, length) > 1]
scales_list$alpha <- NULL
scales_list$blue <- NULL
scales_list$green <- NULL
scales_list$plea <- paste0("alpha", 1:4)
scales_list$safe <- paste0("alpha", 6:11)
scales_list$green <- c("park", "play", "sport", "forest")
scales_list$blue <- c("river", "lake", "beach")
scales_list$noi <- c("traffic", "inside", "outside")

# Reverse code variables --------------------------------------------------

rev_codes <- c("coh4", "coh5", "alpha1", "alpha3", "traffic", "pss4", "pss5", "pss7", "pss8")

#describe(df[, rev_codes])
#head(df[, rev_codes])
df[rev_codes] <- lapply(df[rev_codes], function(x){(max(x, na.rm = TRUE)+1)-x})
#describe(df[, rev_codes])
#head(df[, rev_codes])

#grep("rev$", names(df), value = T)
factor_vars <- sapply(df, function(x){!is.null(attr(x, "labels"))}) & !names(df) %in% unlist(scales_list)

df[factor_vars] <- lapply(df[factor_vars], haven::as_factor)

labelled_vars <- sapply(df, class) == "haven_labelled"
df[labelled_vars] <- lapply(df[labelled_vars], as.numeric)

df <- as.data.frame(df)


# Make parcels ------------------------------------------------------------

parcels <- parcel_items(df, scales_list)

df <- cbind(df, parcels$df_parcels)
df$phys <- as.numeric(df$physical_activity_7days)


# Add extra variables -----------------------------------------------------

# Cluster

df$clus <- df$idpc4i
df$work <- df$working_hours
#df$Dwork <- factor(df$working_hours != "Non-working", labels = c("unemployed", "working"))
df$Dwork <- factor(df$emp_status == "Working with paid work")
levels(df$Dwork) <- c("unemployed","working")
table(df$Dwork, df$emp_status)
table(df$Dwork, df$work_location)
df$Dworkhome <- df$Dwork
df$Dworkhome[which(df$work_location == "At own home")] <- "unemployed"
levels(df$Dworkhome) <- c("home", "away")
table(df$Dworkhome, df$work_location)

#table(df$Dwork, useNA = "always")  

# Exposure

df$exp <- df$traveltime
df$exp <- factor(df$exp == "Less than 15 minutes", labels = c("LO", "HI"))

# Urbanity

df$urb <- log(df$address_density)

#age sex, education, income, household type, marital status
df$edu <- df$education_level
levels(df$edu) <- c("L", "M", "H", "Unknown")
df$edu[df$edu == "Unknown"] <- NA
#table(df$edu)
df$dep <- df$depri16^(1/3)
df$frag <- df$frag16

df$inc <- df$income_quintile
df$inc[df$inc == "Unknown"] <- NA
levels(df$inc) <- c("1", "2", "3", "4", "5", "Unknown")
#df$income <- as.numeric(df$income)

#df$household <- ordered(df$household_type_recode, levels = levels(df$household_type_recode)[c(2,1,3,4)])
df$hous <- relevel(df$household_type_recode, ref="Couple without children")
levels(df$hous) <- c("2", "2c", "sp", "o")

df$mar <- df$marital_status
levels(df$mar) <- c("M", "S", "W", "N")

control_dummies <- c("edu", "inc", "hous", "mar")
control_vars <- c("age", control_dummies)
#table(df$household)
missing <- rowSums(is.na(df[, control_dummies])) > 0
sum(missing)

modmat <- model.matrix(as.formula(paste0("~", paste0(control_dummies, collapse = "+"))), data = df)[, -1]
dummies <- matrix(NA, nrow = nrow(df), ncol = ncol(modmat))
dummies[!missing, ] <- modmat
dummies <- data.frame(dummies)
names(dummies) <- colnames(modmat)
dummies <- dummies[, -which(apply(dummies, 2, sd, na.rm = TRUE) == 0)]

df <- data.frame(df, dummies)
                              
# Randomly select training sample of clusters -----------------------------

set.seed(79)
train <- sample(unique(df$idpc4i), size = .5*length(unique(df$idpc4i)), replace = FALSE)

df_train <- df[df$idpc4i %in% train, ]

#df_train <- df_train[df_train$idpc4i %in% names(table(df_train$idpc4i)[table(df_train$idpc4i) > #11]), ]




# Measurement model -------------------------------------------------------

prepareMplusData(df_train, "try_update.dat", inpfile = TRUE)

if(run_everything){
  cfa_all <- mplusModeler(mplusObject(
    MODEL = syntax_cfa_mplus(parcels$scales_list),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = unlist(parcels$scales_list)
  ), modelout = "cfa_all_update.inp", run = 1L)$results
} else {
  cfa_all <- readModels("cfa_all_update.out")
}

SummaryTable(cfa_all, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


# Measurement with clustered SE -------------------------------------------

if(run_everything){
  cfa_clus <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = syntax_cfa_mplus(parcels$scales_list),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), "clus")
  ), modelout = "cfa_clus_update.inp", run = 1L)$results
} else {
  cfa_clus <- readModels("cfa_clus_update.out")
}

SummaryTable(cfa_clus, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


# 1. Get measurement model working
# 2. Add regression model
# 3. Add area level control variables

# Regression model --------------------------------------------------------

vars_x <- c("noi", "plea", "safe", "coh", "green", "blue")

vars_y <- c("phq", "pss", "phys")

regs <- c(paste0(rep(vars_y, each = length(vars_x)), " ON ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", 2), " ON ", c("pss", "phys"), sep = ";\n", collapse = ""),
          "MODEL INDIRECT:\n\n",
          paste0(rep("phq", each = length(vars_x)), " IND pss ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", each = length(vars_x)), " IND phys ", vars_x, ";\n", collapse = "")
          )

if(run_everything){
  cfa_reg <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = c(syntax_cfa_mplus(parcels$scales_list), regs),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus"))
  ), modelout = "cfa_reg_update.inp", run = 1L)$results
} else {
  cfa_reg <- readModels("cfa_reg_update.out")
}

SummaryTable(cfa_reg, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


# Control variables -------------------------------------------------------
control_vars <- c("dep", "frag", "age", names(dummies))
control <- paste0(rep(c("phq", "pss", "phys"), each = length(control_vars)), " ON ", rep(control_vars, 3), sep = ";\n", collapse = "")

  
if(run_everything){
  cfa_contr <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = c(syntax_cfa_mplus(parcels$scales_list), control, regs),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars))
  ), modelout = "cfa_contr_update.inp", run = 1L)$results
} else {
  cfa_contr <- readModels("cfa_contr_update.out")
}

SummaryTable(cfa_contr, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")



# Model moderated ---------------------------------------------------------


model <- paste0(c(syntax_cfa_mplus(parcels$scales_list), control, regs), collapse = "\n")

model <- strsplit(model, ";\n")[[1]]
model <- gsub("\n", "", model)
model <- model[1:(grep("^MODEL INDIRECT", model)-1)]
model <- paste0(model, ";")
model_label <- !grepl("@", model)
model[model_label] <- syntax_label(model[model_label])

tmp <- c(model, 
         regs[3:length(regs)],
         "MODEL LO:\n", 
         gsub("\\)", "L\\)", model[82:101]),
         "MODEL HI:\n", 
         gsub("\\)", "H\\)", model[82:101]))

if(run_everything){
  sem_mod <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = exp (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = tmp,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "exp")
  ), modelout = "sem_mod_update.inp", run = 1L)$results
} else {
  sem_mod <- readModels("sem_mod_update.out")
}

SummaryTable(sem_mod, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


# Model moderated continuous ---------------------------------------------------

vars_x <- c("noi", "plea", "safe", "coh", "green", "blue")
ints <- paste0("I", vars_x, " | work XWITH ", vars_x, ";")
vars_x <- c(vars_x, "work", paste0("I", vars_x))



vars_y <- c("phq", "pss", "phys")

regs_int <- c(paste0(rep(vars_y, each = length(vars_x)), " ON ", vars_x, ";"),
              paste0(rep("phq", 2), " ON ", c("pss", "phys"), sep = ";")
)

model <- c(syntax_cfa_mplus(parcels$scales_list), control, ints, regs_int)

if(run_everything){
  sem_mod_work <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;",
    ANALYSIS = "TYPE = complex random;",
    MODEL = model,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "work")
  ), modelout = "sem_mod_work_update.inp", run = 1L)$results
} else {
  sem_mod <- readModels("sem_mod_update.out")
}

SummaryTable(sem_mod, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")



# Binary moderated by work hours ------------------------------------------

# OK – let’s stick to just controlling for income. Deprivation and fragmentation should also be included as indicators of area-level SES.

# So, could you please run the model as before but using working hours as a binary variable, comparing non-working and working (0 hours v all other categories)?
vars_x <- c("noi", "plea", "safe", "coh", "green", "blue")

vars_y <- c("phq", "pss", "phys")

regs <- c(paste0(rep(vars_y, each = length(vars_x)), " ON ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", 2), " ON ", c("pss", "phys"), sep = ";\n", collapse = ""),
          "MODEL INDIRECT:\n\n",
          paste0(rep("phq", each = length(vars_x)), " IND pss ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", each = length(vars_x)), " IND phys ", vars_x, ";\n", collapse = "")
)

control_vars <- c("dep", "frag", "age", names(dummies))
control <- paste0(rep(c("phq", "pss", "phys"), each = length(control_vars)), " ON ", rep(control_vars, 3), sep = ";\n", collapse = "")

model <- paste0(c(syntax_cfa_mplus(parcels$scales_list), control, regs), collapse = "\n")

model <- strsplit(model, ";\n")[[1]]
model <- gsub("\n", "", model)
model <- model[1:(grep("^MODEL INDIRECT", model)-1)]
model <- paste0(model, ";")
model_label <- !grepl("@", model)
model[model_label] <- syntax_label(model[model_label])

constraints <- toupper(gsub(" ", "\\.", gsub(" \\(.*$", "", model[grepl("^(phq|pss|phys) ON", model)])))
names(constraints) <- paste0("p", 1:length(constraints))
attr(constraints, "labs") <- gsub("^.+\\((.+)\\);$", "\\1", model[grepl("^(phq|pss|phys) ON", model)])
model_constraints <- c("MODEL CONSTRAINT:", "NEW (", names(constraints), ");", paste0(names(constraints), " = ", paste0(attr(constraints, "labs"), "H"), " - ", paste0(attr(constraints, "labs"), "L"), ";"))

tmp <- c(model, 
         regs[3:length(regs)],
         "MODEL LO:\n", 
         gsub("\\)", "L\\)", model[grepl("^(phq|pss|phys) ON", model)]),
         "MODEL HI:\n", 
         gsub("\\)", "H\\)", model[grepl("^(phq|pss|phys) ON", model)]),
         model_constraints)

if(run_everything){
  sem_mod_dworkconstraints <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dwork (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = tmp,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dwork")
  ), modelout = "sem_mod_dworkconstraints.inp", run = 1L)$results
} else {
  sem_mod_dwork <- readModels("sem_mod_dwork.out")
}

SummaryTable(sem_mod_dworkconstraints, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


results <- printResultsTable(sem_mod_dworkconstraints, keepCols = c("label", "est_sig", "confint", "Group"))
eval_constr <- results[grepl("^New.Add", results$label), ]
results <- results[!grepl("^New.Add", results$label), ]
results <- cbind(results[results$Group == "LO", 1:3], results[results$Group == "HI", 2:3])

eval_constr$label <- constraints
names(results) <- c("Parameter", "est_un", "CI_un", "est_emp", "CI_emp")
results$D_CI <- results$D_est <- NA
results[match(eval_constr$label, results$Parameter), c("D_est", "D_CI")] <- eval_constr[, c("est_sig", "confint")]

write.csv(results, "moderated_Dwork.csv", na = "")

write.csv(SummaryTable(list(cfa_contr, sem_mod_dworkconstraints), keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename"), "fits.csv")



# 30-11-2019 additions ----------------------------------------------------

# -	A covariance between the two mediators (stress + physical activity) needs to be added
# -	Urbanity (as a control variable) is missing – this is the variable ‘address_density’
# -	Please use ‘emp_status’ with 1 vs all other categories as the two groups with which to run the model. Since we decided to work with two categories, I find it cleaner to use employment status rather than working hours
# 
# I also want to check a couple of points. Did we reverse code the ‘traffic’ variable? I know we checked the stress variables, but I can’t recall about traffic. And also, regarding the physical activity mediator, which variable are you using? There are two in the file I sent – ‘physical_activity_normal’ and ‘physical_activity_7days’. Can you let me know?
#   
#   Lastly, can you send me both the standardized and unstandardized results when it’s done. It’s easier for my interpretation, and I can finalise the tables.


vars_x <- c("noi", "plea", "safe", "coh", "green", "blue")

vars_y <- c("phq", "pss", "phys")

regs <- c(paste0(rep(vars_y, each = length(vars_x)), " ON ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", 2), " ON ", c("pss", "phys"), sep = ";\n", collapse = ""),
          "MODEL INDIRECT:\n\n",
          paste0(rep("phq", each = length(vars_x)), " IND pss ", vars_x, ";\n", collapse = ""),
          paste0(rep("phq", each = length(vars_x)), " IND phys ", vars_x, ";\n", collapse = "")
)

control_vars <- c("dep", "frag", "age", "urb", names(dummies))
control <- paste0(rep(c("phq", "pss", "phys"), each = length(control_vars)), " ON ", rep(control_vars, 3), sep = ";\n", collapse = "")

model <- paste0(c(syntax_cfa_mplus(parcels$scales_list), control, regs), collapse = "\n")

model <- strsplit(model, ";\n")[[1]]
model <- gsub("\n", "", model)
model <- model[1:(grep("^MODEL INDIRECT", model)-1)]
model <- paste0(model, ";")
model_label <- !grepl("@", model)
model[model_label] <- syntax_label(model[model_label])

constraints <- toupper(gsub(" ", "\\.", gsub(" \\(.*$", "", model[grepl("^(phq|pss|phys) ON", model)])))
names(constraints) <- paste0("p", 1:length(constraints))
attr(constraints, "labs") <- gsub("^.+\\((.+)\\);$", "\\1", model[grepl("^(phq|pss|phys) ON", model)])
model_constraints <- c("MODEL CONSTRAINT:", "NEW (", names(constraints), ");", paste0(names(constraints), " = ", paste0(attr(constraints, "labs"), "H"), " - ", paste0(attr(constraints, "labs"), "L"), ";"))
saveRDS(constraints, "constraints.RData")

model3011 <- c(model, 
         "pss WITH phys;",
         regs[3:length(regs)],
         "MODEL LO:\n", 
         gsub("\\)", "L\\)", model[grepl("^(phq|pss|phys) ON", model)]),
         "MODEL HI:\n", 
         gsub("\\)", "H\\)", model[grepl("^(phq|pss|phys) ON", model)]),
         model_constraints)

if(run_everything){
  sem_mod_dworkconstraints3011 <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dwork (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = model3011,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dwork")
  ), modelout = "sem_mod_dworkconstraints_3011.inp", run = 1L)$results
} else {
  sem_mod_dworkconstraints3011 <- readModels("sem_mod_dworkconstraints_3011.out")
}

SummaryTable(sem_mod_dworkconstraints3011, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")


# Check modification indices ----------------------------------------------

sem_mod_dworkconstraints3011mod <- readModels("sem_mod_dworkconstraints_3011_modindices.out")

mod <- sem_mod_dworkconstraints3011mod$mod_indices
head(mod[order(mod$MI, decreasing = T), ], 20)

model3011_mod <- append(model3011, "safe BY traffic (safeBtraffic);", after = 26)
model3011_mod <- append(model3011_mod, "PPHQ4 WITH PPHQ3 (phq4Wphq3);", after = 19)

if(run_everything){
  sem_mod_dworkconstraints3011_mod <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dwork (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = model3011_mod,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dwork")
  ), modelout = "sem_mod_dworkconstraints_3011_mod.inp", run = 1L)$results
} else {
  sem_mod_dworkconstraints3011_mod <- readModels("sem_mod_dworkconstraints_3011_mod.out")
}

SummaryTable(sem_mod_dworkconstraints3011_mod, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")


results <- printResultsTable(sem_mod_dworkconstraints3011_mod, keepCols = c("label", "est_sig", "confint", "Group"))
eval_constr <- results[grepl("^New.Add", results$label), ]
results <- results[!grepl("^New.Add", results$label), ]
results <- cbind(results[results$Group == "LO", 1:3], results[results$Group == "HI", 2:3])

eval_constr$label <- constraints
names(results) <- c("Parameter", "est_un", "CI_un", "est_emp", "CI_emp")
results$D_CI <- results$D_est <- NA
results[match(eval_constr$label, results$Parameter), c("D_est", "D_CI")] <- eval_constr[, c("est_sig", "confint")]

write.csv(results, "moderated_Dwork_unstandardized.csv", na = "")

results_std <- printResultsTable(sem_mod_dworkconstraints3011_mod, parameters = "stdyx.standardized", keepCols = c("label", "est_sig", "confint", "Group"))
results_std <- results_std[!grepl("^New.Add", results_std$label), ]
results_std <- cbind(results_std[results_std$Group == "LO", 1:3], results_std[results_std$Group == "HI", 2:3])
names(results_std) <- c("Parameter", "est_un_std", "CI_un_std", "est_emp_std", "CI_emp_std")

results_std_out <- results
results_std_out[, c(2:5)] <- results_std[, c(2:5)]

write.csv(results_std_out, "moderated_Dwork_STDYXstandardized.csv", na = "")

write.csv(SummaryTable(list(cfa_contr, sem_mod_dworkconstraints3011_mod), keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename"), "fits.csv")

# Add work-from-home to unemployed group ----------------------------------

if(run_everything){
  sem_mod_1112 <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dworkhome (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = model3011_mod,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dworkhome")
  ), modelout = "sem_mod_1112.inp", run = 1L)$results
} else {
  sem_mod_1112 <- readModels("sem_mod_1112.out")
}

SummaryTable(sem_mod_1112, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")

results <- printResultsTable(sem_mod_1112, keepCols = c("label", "est_sig", "confint", "Group"))
eval_constr <- results[grepl("^New.Add", results$label), ]
results <- results[!grepl("^New.Add", results$label), ]
results <- cbind(results[results$Group == "LO", 1:3], results[results$Group == "HI", 2:3])

eval_constr$label <- constraints
names(results) <- c("Parameter", "est_un", "CI_un", "est_emp", "CI_emp")
results$D_CI <- results$D_est <- NA
results[match(eval_constr$label, results$Parameter), c("D_est", "D_CI")] <- eval_constr[, c("est_sig", "confint")]

write.csv(results, "moderated_1112_unstandardized.csv", na = "")

results_std <- printResultsTable(sem_mod_1112, parameters = "stdyx.standardized", keepCols = c("label", "est_sig", "confint", "Group"))
results_std <- results_std[!grepl("^New.Add", results_std$label), ]
results_std <- cbind(results_std[results_std$Group == "LO", 1:3], results_std[results_std$Group == "HI", 2:3])
names(results_std) <- c("Parameter", "est_un_std", "CI_un_std", "est_emp_std", "CI_emp_std")

results_std_out <- results
results_std_out[, c(2:5)] <- results_std[, c(2:5)]

write.csv(results_std_out, "moderated_1112_STDYXstandardized.csv", na = "")

write.csv(SummaryTable(list(sem_mod_1112), keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename"), "fits1112.csv")



# Evaluate on test data ---------------------------------------------------
df_test <- df[!df$idpc4i %in% train, ]

if(run_everything){
  cfa_clus_test <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = syntax_cfa_mplus(parcels$scales_list),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_test,
    usevariables = c(unlist(parcels$scales_list), "clus")
  ), modelout = "cfa_test.inp", run = 1L)$results
} else {
  cfa_clus_test <- readModels("cfa_test.out")
}

SummaryTable(cfa_clus_test, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


if(run_everything){
  sem_mod_1112_test <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dworkhome (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = model3011_mod,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_test,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dworkhome")
  ), modelout = "sem_mod_1112_test.inp", run = 1L)$results
} else {
  sem_mod_1112_test <- readModels("sem_mod_1112_test.out")
}

SummaryTable(sem_mod_1112_test, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")


# Evaluate on complete data for parameter estimation ----------------------

df_full <- df

if(run_everything){
  cfa_clus_full <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = syntax_cfa_mplus(parcels$scales_list),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_full,
    usevariables = c(unlist(parcels$scales_list), "clus")
  ), modelout = "cfa_full.inp", run = 1L)$results
} else {
  cfa_clus_full <- readModels("cfa_full.out")
}

SummaryTable(cfa_clus_full, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")


if(run_everything){
  sem_mod_1112_full <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dworkhome (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = model3011_mod,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_full,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dworkhome")
  ), modelout = "sem_mod_1112_full.inp", run = 1L)$results
} else {
  sem_mod_1112_full <- readModels("sem_mod_1112_full.out")
}

SummaryTable(sem_mod_1112_full, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")


# Prepare output for manuscript -------------------------------------------

models <- list(
  cfa_clus,
  cfa_clus_test,
  cfa_clus_full,
  sem_mod_1112,
  sem_mod_1112_test,
  sem_mod_1112_full
)

fit_tab <- SummaryTable(models, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR"), sortBy = "Filename")
fit_tab$Data <- "Train"
fit_tab$Data[grepl("test", fit_tab$Filename)] <- "Test"
fit_tab$Data[grepl("full", fit_tab$Filename)] <- "Full"
fit_tab$Model <- "CFA"
fit_tab$Model[grepl("sem_mod", fit_tab$Filename)] <- "SEM"
fit_tab$Filename <- NULL

write.csv(fit_tab[c(1,3,2,4,6,5), c("Model", "Data", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "SRMR")], "fits_23-1-2020.csv")



# Make graph --------------------------------------------------------------
layout <- read.csv("layout2020.csv", header = FALSE, stringsAsFactors = FALSE)
edg <- tidySEM:::get_edges.mplus.object(sem_mod_1112_full)
edg$group <- c("LO" = "Home", "HI" = "Away")[edg$group]
edg <- edg[!edg$arrow == "both", ]
nod <- get_nodes(sem_mod_1112_full)
nod$group <- c("LO" = "Home", "HI" = "Away")[nod$group]
nod <- nod[nod$name %in% unlist(layout),]
nod$label <- gsub("\\n.*$", "", nod$label)

labels <- c("GREEN" = "Green space",
               "BLUE" = "Blue space",
               "PLEA" = "Pleasantness",
               "NOI" = "Environmental\ndisturbance",
               "SAFE"  = "Safety",
               "COH" = "Social cohesion",
               "PHYS" = "Physical activity",
               "PSS" = "Stress",
               "PHQ" ="Depression"
)
nod$label <- labels[nod$label]
edg <- edg[edg$from %in% unlist(layout) & edg$to %in% unlist(layout), ]
#edg$curvature[edg$curvature == 60] <- -20
#edg <- edg[is.na(edg$curvature), ]

prep <- prepare_graph(edges = edg, layout = layout, nodes = nod, text_size = 3.5, ellipses_height = 1.3, spacing_y = 2.5, angle = 1)

edges(prep)$label_location <- .5
edges(prep)[edges(prep)$to == "PHYS" & edges(prep)$from %in% c("NOI", "SAFE", "COH"), ]$label_location <- .2
edges(prep)[edges(prep)$to == "PSS" & edges(prep)$from %in% c("GREEN", "BLUE", "PLEA"), ]$label_location <- .2
edges(prep) <- edges(prep)[(is.na(edges(prep)$curvature) | (edges(prep)$from =="PSS" & edges(prep)$to =="PHYS")),]
edges(prep)$curvature[!is.na(edges(prep)$curvature)] <- -75
edges(prep)$connect_from[!is.na(edges(prep)$curvature)] <- "top"
edges(prep)$connect_to[!is.na(edges(prep)$curvature)] <- "bottom"


p <- plot(prep)
p <- p + facet_grid(group~., labeller = labeller(group = c('Away'="Low strength of exposure", 'Home'="High strength of exposure")))
p
library(ggplot2)
ggsave("Standardized_regs.png", p, width = 8, height = 10)
ggsave("Standardized_regs.pdf", p, width = 8, height = 10)


# Get indirect effect specification ---------------------------------------
total_indirect <- model3011_mod
total_indirect <- total_indirect[-115] # Remove one of the indirect specifications
total_indirect[115] <- gsub(" phys", "", total_indirect[115])
if(run_everything){
  sem_mod_indirect <- mplusModeler(mplusObject(
    VARIABLE = "CLUSTER = clus;\nGROUPING = Dworkhome (1 = LO 2 = HI);\n",
    ANALYSIS = "TYPE = COMPLEX;\n",
    MODEL = total_indirect,
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_full,
    usevariables = c(unlist(parcels$scales_list), c("phys", "clus", control_vars), "Dworkhome")
  ), modelout = "sem_mod_indirect.inp", run = 1L)$results
} else {
  sem_mod_indirect <- readModels("sem_mod_indirect.out")
}


# Table results -----------------------------------------------------------

results <- tidySEM:::table_results(sem_mod_1112_full, columns = c("label", "est_sig", "confint", "group"))
eval_constr <- results[grepl("^New.Add", results$label), ]
results <- results[!grepl("^New.Add", results$label), ]
results <- cbind(results[results$group == "LO", 1:3], results[results$group == "HI", 2:3])
all(gsub("(LO|HI)$", "", results[[1]]) == gsub("(LO|HI)$", "", results[[4]]))

constraints <- readRDS("constraints.RData")
eval_constr$label <- constraints[tolower(gsub("^New\\.Additional\\.Parameters\\.(.+?)\\..*$", "\\1", eval_constr$label))]
names(results) <- c("Parameter", "est_un", "CI_un", "est_emp", "CI_emp")
results$D_CI <- results$D_est <- NA
results$Parameter <- gsub("\\.LO$", "", results$Parameter)
results[match(eval_constr$label, results$Parameter), c("D_est", "D_CI")] <- eval_constr[, c("est_sig", "confint")]

results_total_ind <- tidySEM:::table_results(sem_mod_indirect, columns = c("label", "est_sig", "confint", "group"))
results_total_ind <- results_total_ind[grepl("^Total\\.(?!indirect)", results_total_ind$label, perl = TRUE), ]
results_total_ind$label <- gsub("\\.(LO|HI)$", "", results_total_ind$label)
results_total_ind <- reshape(results_total_ind, v.names = names(results_total_ind)[2:3], timevar = "group", idvar = "label", direction = "wide")
results_total_ind <- cbind(results_total_ind[, c(1,4,5,2,3)], NA, NA)
names(results_total_ind) <- names(results)
results <- rbind(results, results_total_ind)

write.csv(results, "final_results_unstandardized.csv", na = "")

results_std <- tidySEM::table_results(sem_mod_1112_full, columns = c("label", "est_sig_std", "confint_std", "group"))
results_std <- results_std[!grepl("^New.Add", results_std$label), ]
results_std <- cbind(results_std[results_std$group == "LO", 1:3], results_std[results_std$group == "HI", 2:3])
names(results_std) <- c("Parameter", "est_un_std", "CI_un_std", "est_emp_std", "CI_emp_std")

results_total_ind <- tidySEM:::table_results(sem_mod_indirect, columns = c("label", "est_sig_std", "confint_std", "group"))
results_total_ind <- results_total_ind[grepl("^Total\\.(?!indirect)", results_total_ind$label, perl = TRUE), ]
results_total_ind$label <- gsub("\\.(LO|HI)$", "", results_total_ind$label)
results_total_ind <- reshape(results_total_ind, v.names = names(results_total_ind)[2:3], timevar = "group", idvar = "label", direction = "wide")
results_total_ind <- results_total_ind[, c(1,4,5,2,3)]
names(results_total_ind) <- names(results_std)
results_std <- rbind(results_std, results_total_ind)

results_std_out <- results
results_std_out[, c(2:5)] <- tail(results_std[, c(2:5)])

write.csv(results_std_out, "final_results_standardized.csv", na = "")



# Add correlation tables 27-05-2020 ---------------------------------------

cors <- table_cors(sem_mod_1112_full)

lapply(names(cors), function(x){
  write.csv(cors[[x]], paste0("correlation_table_", x, ".csv"))
})

library(semTools)
meas_inv <- lapply(names(scales_list), function(scalename){
  #scalename = names(scales_list)[1]
  mod <- paste0(scalename, " =~ ", paste0(scales_list[[scalename]], collapse = " + "))
  syntax_config <- measEq.syntax(configural.model = mod, 
                                 data = df, 
                                 group = "Dworkhome")
  res_config <- cfa(as.character(syntax_config), data = df, group = "Dworkhome")
  syntax_metric <- measEq.syntax(configural.model = mod, 
                                 data = df, 
                                 group = "Dworkhome",
                                 group.equal = "loadings")
  res_metric <- cfa(as.character(syntax_metric), data = df, group = "Dworkhome")
  syntax_scalar <- measEq.syntax(configural.model = mod, 
                                 data = df, 
                                 group = "Dworkhome",
                                 group.equal = c("loadings", "intercepts"))
  res_scalar <- cfa(as.character(syntax_scalar), data = df, group = "Dworkhome")
  cf <- compareFit(res_config, res_metric, res_scalar)
  
  out1 <- cf@nested
  out1$p <- tidySEM::est_sig(x = out1$`Pr(>Chisq)`, sig = out1$`Pr(>Chisq)`)
  out1$`Pr(>Chisq)` <- NULL
  out1 <- cbind(variable = scalename, model = rownames(out1), out1, cf@fit[c("rmsea", "cfi", "tli", "srmr")])
  deltas <- cf@fit[c("rmsea", "cfi", "tli", "srmr")]
  colnames(deltas) <- paste0("D", colnames(deltas))
  deltas[3,] <- deltas[3,]-deltas[2,]
  deltas[2,] <- deltas[2,]-deltas[1,]
  deltas[1, ] <- NA
  out1 <- cbind(out1, deltas)

  tidySEM:::format_numeric(out1)
})

meas_tab <- do.call(rbind, meas_inv)
meas_tab$model <- gsub("^res_", "", meas_tab$model)
#meas_tab <- meas_tab[!meas_tab$model == "res_scalar", ]
write.csv(meas_tab, "measurement_invariance.csv", row.names = FALSE)