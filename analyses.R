library(haven)
library(psych)
library(motley)
library(MplusAutomation)
source("c:/git_repositories/miniMAC/factor_syntax.R")
source("c:/git_repositories/miniMAC/rfunctions.R")
df <- read_dta("sem_data.dta")

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
scales_list$noi <- c("traffic", "noise", "inside", "outside")

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

#length(table(df$municipality))
#hist(table(df$municipality)[!table(df$municipality) == 0])
#length(table(df$idpc4i))
#table(table(df$idpc4i))

# Randomly select training sample of clusters -----------------------------

set.seed(79)
train <- sample(unique(df$idpc4i), size = .5*length(unique(df$idpc4i)), replace = FALSE)

df_train <- df[df$idpc4i %in% train, ]

df_train <- df_train[df_train$idpc4i %in% names(table(df_train$idpc4i)[table(df_train$idpc4i) > 11]), ]

df_train$clus <- df_train$idpc4i
length(unique(df_train$idpc4i))

# Measurement model -------------------------------------------------------

prepareMplusData(df_train, "try.dat", inpfile = TRUE)

if(run_everything){
  cfa_all <- mplusModeler(mplusObject(
    VARIABLE = paste0("CATEGORICAL = \n", paste0(strwrap(paste0(unlist(scales_list), collapse = " "), width = 70), collapse = "\n"), ";"),
    #ANALYSIS = "bootstrap=1000;",
    MODEL = gsub("\\[.+?\\];", "", syntax_cfa_mplus(scales_list, standardize = TRUE)),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = unlist(scales_list)
  ), modelout = "cfa_all.inp", run = 1L)$results
} else {
  cfa_all <- readModels("cfa_all.out")
}

SummaryTable(cfa_all, keepCols = c("Filename", "LL", "Parameters", "RMSEA_Estimate", "CFI", "TLI", "WRMR"), sortBy = "Filename")

if(run_everything){
  cfa_2lev <- mplusModeler(mplusObject(
    VARIABLE = paste0("CATEGORICAL = \n",
                      paste0(strwrap(paste0(unlist(scales_list), collapse = " "), width = 70), collapse = "\n"),
                      ";\n WITHIN =  \n",
                      paste0(strwrap(paste0(unlist(scales_list), collapse = " "), width = 70), collapse = "\n"),
                      ";\nCLUSTER = clus;"),
    ANALYSIS = "TYPE = TWOLEVEL;\nINTEGRATION=MONTECARLO;",
    MODEL = paste0("%WITHIN%\n", gsub("\\[.+?\\];", "", syntax_cfa_mplus(scales_list, standardize = TRUE))),
    OUTPUT = "standardized tech4 modindices;",# CINTERVAL (BCBOOTSTRAP);",
    rdata = df_train,
    usevariables = c(unlist(scales_list), "clus")
  ), modelout = "cfa_all_2lev.inp", run = 1L)$results
  cfa_2lev$errors
} else {
  cfa_all <- readModels("cfa_all.out")
}

# 1. Get measurement model working
# 2. Add regression model
# 3. Add area level control variables

# Regression model --------------------------------------------------------


