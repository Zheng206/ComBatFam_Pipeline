suppressMessages(library(argparser))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(ComBatFamily))
suppressMessages(library(ComBatFamQC))
suppressMessages(library(lme4))
suppressMessages(library(mgcv))
suppressMessages(library(parallel))

## Read in arguments
p <- arg_parser("Apply ComBatFamily harmonization method to remove bathch effect", hide.opts = FALSE)
p <- add_argument(p, "data", help = "path to the CSV or EXCEL file that contains data to be harmonized, covariates and batch information")
p <- add_argument(p, "--visualization", short = '-v', help = "a boolean variable indicating whether to run harmonization or batch effect visualization. TRUE indicates batch effect visualization.", default = TRUE)
p <- add_argument(p, "--after", short = '-a', help = "a boolean variable indicating whether the dataset is before or after harmonization. TRUE indicates the dataset is after harmonization.", default = FALSE)
p <- add_argument(p, "--family", help = "combat family to use, comfam or covfam", default = "comfam")
#p <- add_argument(p, "--cov", help = "whether covbat is considered for harmonization, TRUE or FALSE", default = FALSE)
p <- add_argument(p, "--features", short = '-f', help = "position of data to be harmonized(column numbers), eg: 1-5,9")
p <- add_argument(p, "--covariates", short = '-c', help = "position of covariates (column numbers)", default = "NULL")
p <- add_argument(p, "--batch", short = '-b', help = "position of batch column (column number)")
p <- add_argument(p, "--model", short = '-m', help = "select the model function for harmonization, eg: lm, gam", default = "lm")
p <- add_argument(p, "--smooth", short = '-s', help = "provide the variables that require a smooth function", default = "NULL")
p <- add_argument(p, "--interaction", help = "specify the interaction terms in the format of col_index1*col_index2, eg 2*3,11*12", default = "NULL")
p <- add_argument(p, "--int_type", help = "specify an interaction type for gam models, eg: linear, categorical-continuous, factor-smooth, tensor, smooth-smooth", default = "linear")
p <- add_argument(p, "--random", short = '-r', help = "specify the random intercept-effects", default = "NULL")
p <- add_argument(p, "--eb", help = "whether to use ComBat model with empirical Bayes for mean and variance", default = TRUE)
p <- add_argument(p, "--score_eb", help = "whether to use ComBat model with empirical Bayes for score mean and variance", default = FALSE)
p <- add_argument(p, "--robust.LS", help = "whether to use robust location and scale estimators for error variance and site effect parameters", default = FALSE)
p <- add_argument(p, "--ref.batch", help = "reference batch", default = "NULL")
p <- add_argument(p, "--score.model", help = "model for scores, defaults to NULL for fitting basic location and scale model without covariates on the scores", default = "NULL")
p <- add_argument(p, "--score.args", help = "list of arguments for score model, input should have the following format arg1,arg2,...", default = "NULL")
p <- add_argument(p, "--percent.var", help = "the number of harmonized principal component scores is selected to explain this proportion of the variance", default = 0.95)
p <- add_argument(p, "--n.pc", help = "if specified, this number of principal component scores is harmonized", default = "NULL")
p <- add_argument(p, "--std.var", help = "if TRUE, scales variances to be equal to 1 before PCA", default = TRUE)
p <- add_argument(p, "--predict", help = "a boolean variable indicating whether to run ComBat from scratch or apply existing model to new dataset (currently only work for original ComBat and ComBat-GAM)", default = FALSE)
p <- add_argument(p, "--object", help = "path to the ComBat model to be used for prediction, should have an extension of .rds", default = "NULL")
p <- add_argument(p, "--reference", help = "path to the CSV or EXCEL file that contains dataset to be considered as the reference group", default = "NULL")
p <- add_argument(p, "--outdir", short = '-o', help = "full path (including the file name) where harmonized data should be written")
p <- add_argument(p, "--mout", help = "full path where ComBat model to be saved")
p <- add_argument(p, "--cores", help = "number of cores used for paralleling computing, please provide a numeric value", default = "all")
argv <- parse_args(p)

# Useful Function
form_gen = function(x, c, i, random, smooth, int_type){
  if (x == "lm"){
    if(!is.null(c)){
      if (is.null(i)){form = paste0("y ~", paste(colnames(c), collapse = "+"))}else{
        i = gsub(",", "*", i)
        form = paste0("y ~", paste(colnames(c), collapse = "+"),  " + ", paste(i, collapse = " + "))
      }
    }else{form = NULL}
  }else if (x == "lmer"){
    if(!is.null(c)){
      if (is.null(i)){form = paste0("y ~", paste(colnames(c), collapse = " + "),  " + ", paste("(1 |", random, ")", collapse = " + "))}else{
        i = gsub(",", "*", i)
        form = paste0("y ~", paste(colnames(c), collapse = " + "),  " + ", paste(i, collapse = " + "), " + ", paste("(1 |", random, ")", collapse = " + "))
      }
    }else{form = paste0("y ~", paste("(1 |", random, ")", collapse = " + "))}
  }else if(x == "gam"){
    c = c[setdiff(colnames(c), smooth)]
    if(int_type == "linear"){
      i = gsub(",", "*", i)
    }else if(int_type == "factor-smooth"){
      i = paste("s(", i, ", bs = 'fs')")
      smooth_index = sapply(smooth, function(x){
        if(sum(grepl(x, i)) > 0){return(1)}else(return(0))
      }, USE.NAMES = FALSE)
      cov_index = sapply(colnames(c), function(x){
        if(sum(grepl(x, i)) > 0){return(1)}else{return(0)}
      }, USE.NAMES = FALSE)
      smooth = smooth[which(smooth_index == 0)]
      c = c[which(cov_index == 0)]
    }else if(int_type == "tensor"){
      i = paste("ti(", i, ")")}
    
    if (is.null(i)){form = paste0("y ~ ", paste(colnames(c), collapse = " + "), " + ", paste("s(", smooth, ")", collapse = " + "))}else{
      if(length(smooth) > 0){
        form = paste0("y ~ ", paste(colnames(c), collapse = " + "), " + ", paste("s(", smooth, ")", collapse = " + "), " + ", paste(i, collapse = " + "))
        }else{
        form = paste0("y ~ ", paste(colnames(c), collapse = " + "), " + ", paste(i, collapse = " + "))
      }
    }
  }
  return(form)
}

# Preprocess inputs
message('Checking inputs...')
if(is.na(argv$data)) stop("Missing input data") else {
  if(!grepl("csv$|xls$", argv$data)) stop("Input file must be a csv or an excel file") else {
    if(grepl("csv$", argv$data)) df = read_csv(argv$data) else df = read_excel(argv$data)
  }
}

df = data.frame(df)

if(is.na(argv$batch)) stop("Please identify the position of batch column") else {
  bat_col = as.numeric(argv$batch)
}

if(is.na(argv$features)) stop("Please identify the position of data to be harmonized") else {
  col = gsub("-",":",argv$features)
  col_vec = eval(parse(text = paste0("c(", col, ")")))
}

if(argv$covariates == "NULL") {
  cov_col = NULL
  covariates = NULL
} else {
  cov_col = gsub("-",":",argv$covariates)
  cov_col = eval(parse(text = paste0("c(", cov_col, ")")))
}

if(argv$model == "gam"){
  if(argv$covariates == "NULL") stop("Please provide covariates for gam model")
  if(argv$smooth == "NULL") stop("Please provide variables that require a smoothing function") else {
    smooth_col = gsub("-",":",argv$smooth)
    smooth_col = eval(parse(text = paste0("c(", smooth_col, ")")))
    smooth_var = colnames(df)[smooth_col]
    #cov_var = colnames(df)[cov_col]
    #cov_var = setdiff(cov_var, smooth_var)
    smooth = smooth_var
  }
}else{
  smooth = eval(parse(text = argv$smooth))
}

if(argv$model == "lmer"){
  if(argv$random == "NULL") stop("Please specify random intercept-effects") else {
    random_col = gsub("-",":",argv$random)
    random_col = eval(parse(text = paste0("c(", random_col, ")")))
    random_var = colnames(df)[random_col]
    random = random_var
  }
}else{
  random_col = eval(parse(text = argv$random))
  random = eval(parse(text = argv$random))
}

if(argv$reference != "NULL"){
  if(!grepl("csv$|xls$", argv$reference)) stop("Input file must be a csv or an excel file") else {
    if(grepl("csv$", argv$reference)) reference_df = read_csv(argv$reference) else reference_df = read_excel(argv$reference)
  }
  reference_df = reference_df[complete.cases(reference_df[c(col_vec, bat_col, cov_col, random_col)]),]
}

if(!argv$object == "NULL"){
  object = readRDS(argv$object)
}

obs_n = nrow(df)
df = df[complete.cases(df[c(col_vec, bat_col, cov_col, random_col)]),]
obs_new = nrow(df)
print(paste0(obs_n - obs_new, " observations are dropped due to missing values."))

if(argv$reference != "NULL"){
  untouched = df %>% semi_join(reference_df)
  new_data = df %>% anti_join(reference_df)
}

# Feature Wranggling
features_orig = df[col_vec]
n_orig = length(colnames(features_orig))
features = features_orig[, apply(features_orig, 2, function(col) { length(unique(col)) > 1 })]
n_new = length(colnames(features))
dropped_col = NULL
if (n_orig > n_new){
  dropped_col = setdiff(colnames(features_orig), colnames(features))
  print(paste0(n_orig - n_new, " univariate features dropped: ", dropped_col))
}
features_col = colnames(features)

# Batch Wranggling
batch = as.factor(df[[bat_col]])
bat_var = colnames(df[bat_col])


# Covariates Wranggling
char_var = cov_col[sapply(df[cov_col], function(col) is.character(col) || is.factor(col))]
enco_var = cov_col[sapply(df[cov_col], function(col) length(unique(col)) == 2 && all(unique(col) %in% c(0,1)))]
df[char_var] =  lapply(df[char_var], as.factor)
df[enco_var] =  lapply(df[enco_var], as.factor)
covariates = df[cov_col]
cov_var = colnames(covariates)

# Interaction Wranggling
if(argv$interaction == "NULL"){interaction = eval(parse(text = argv$interaction))}else{
  interaction_l = lapply(str_split(argv$interaction, ",")[[1]], function(x) str_split(x,  "\\*")[[1]])
  interaction = sapply(interaction_l, function(x){
    x1 = colnames(df)[as.numeric(x[1])]
    x2 = colnames(df)[as.numeric(x[2])]
    element = paste0(x1, ",", x2)
  }, USE.NAMES = FALSE)
}

used_col = c(features_col, cov_var, bat_var)
other_info = df[-which(colnames(df) %in% used_col)]

if (is.null(cov_col)){
  visual_cov = NULL
}else{
  visual_cov = colnames(df[cov_col])
}

# Batch Effect Visualization
if(argv$visualization){
  message("preparing datasets for visualization......")
  if(argv$family == "comfam"){
    result = visual_prep(features = features_col, type = argv$model, batch = colnames(df)[bat_col], covariates = visual_cov, interaction = interaction, random = random, smooth = smooth, smooth_int_type = argv$int_type, df = df, family = argv$family, 
                         eb = argv$eb, 
                         robust.LS = argv$robust.LS, 
                         ref.batch = if(argv$ref.batch == "NULL") eval(parse(text = argv$ref.batch)) else argv$ref.batch, 
                         predict = argv$predict, 
                         object = if(argv$object == "NULL") eval(parse(text = argv$object)) else object,
                         reference = if(argv$reference == "NULL") eval(parse(text = argv$reference)) else reference_df,
                         cores = if(argv$cores == "all") detectCores() else as.numeric(argv$cores))
  }else{
    result = visual_prep(features = features_col, type = argv$model, batch = colnames(df)[bat_col], covariates = visual_cov, interaction = interaction, random = random, smooth = smooth, smooth_int_type = argv$int_type, df = df, family = argv$family, 
                         eb = argv$eb, 
                         score_eb = argv$score_eb,
                         robust.LS = argv$robust.LS, 
                         ref.batch = if(argv$ref.batch == "NULL") eval(parse(text = argv$ref.batch)) else argv$ref.batch,
                         score.model = if(argv$score.model == "NULL") eval(parse(text = argv$score.model)) else argv$score.model,
                         score.args = if(argv$score.args == "NULL") eval(parse(text = argv$score.args)) else eval(parse(text = paste0("list(", argv$score.args, ")"))),
                         percent.var = as.numeric(argv$percent.var),
                         n.pc = if(argv$n.pc == "NULL") eval(parse(text = argv$n.pc)) else as.numeric(argv$n.pc),
                         std.var = argv$std.var,
                         predict = argv$predict, 
                         object = if(argv$object == "NULL") eval(parse(text = argv$object)) else object,
                         reference = if(argv$reference == "NULL") eval(parse(text = argv$reference)) else reference_df,
                         cores = if(argv$cores == "all") detectCores() else as.numeric(argv$cores))
  }
  
  
  if(!is.na(argv$outdir)){
    message("Saving harmonized data......")
    write_csv(result$harmonized_df, argv$outdir)
    message(sprintf("Results saved at %s", argv$outdir))  
  }
  
  if(!is.na(argv$mout)){
    message("Saving ComBat model......")
    saveRDS(result$combat.object, argv$mout)
    message(sprintf("ComBat model saved at %s", argv$mout))  
  }
  
  message("datasets are ready! Starting shiny app......")
  comfam_shiny(result, argv$after, argv$eb, argv$score_eb)
  #print(paste(c(argv$model, colnames(df)[bat_col], colnames(df[cov_col]), interaction, smooth, argv$int_type), collapse = ","))
}else{
# Run Combat
# customize formula
  if (is.null(covariates)){
    form_c = NULL
    combat_c = NULL
  }else{
    if(argv$model == "lmer"){
      form_c = df[cov_col]
      combat_c = cbind(df[cov_col], df[random])
    }else{
      form_c = df[cov_col]
      combat_c = df[cov_col]
    }
  }
  if (argv$reference == "NULL"){
    if (!argv$predict){
      form = form_gen(x = argv$model, c = form_c, i = interaction, random = random, smooth = smooth, int_type = argv$int_type)
      if(argv$family == "comfam"){
            ComBat_run = ComBatFamily::comfam(data = features,
                                              bat = batch , 
                                              covar = combat_c,
                                              model = eval(parse(text = argv$model)), 
                                              formula = as.formula(form),
                                              eb = argv$eb,
                                              robust.LS = argv$robust.LS, 
                                              ref.batch = if(argv$ref.batch == "NULL") eval(parse(text = argv$ref.batch)) else argv$ref.batch
            )
    comf_df = ComBat_run$dat.combat
    comf_df = cbind(other_info, df[bat_col], df[cov_col], comf_df)
    comf_df = comf_df[colnames(df)]
    print(head(comf_df))
    write_csv(comf_df, paste0(argv$outdir))
    }else if(argv$family == "covfam"){
        ComBat_run = ComBatFamily::covfam(data = features,
                                        bat = batch , 
                                        covar = combat_c,
                                        model = eval(parse(text = argv$model)), 
                                        formula = as.formula(form),
                                        eb = argv$eb,
                                        score_eb = argv$score_eb,
                                        robust.LS = argv$robust.LS, 
                                        ref.batch = if(argv$ref.batch == "NULL") eval(parse(text = argv$ref.batch)) else argv$ref.batch,
                                        score.model = if(argv$score.model == "NULL") eval(parse(text = argv$score.model)) else argv$score.model,
                                        score.args = if(argv$score.args == "NULL") eval(parse(text = argv$score.args)) else eval(parse(text = paste0("list(", argv$score.args, ")"))),
                                        percent.var = as.numeric(argv$percent.var),
                                        n.pc = if(argv$n.pc == "NULL") eval(parse(text = argv$n.pc)) else as.numeric(argv$n.pc),
                                        std.var = argv$std.var)
      covf_df = ComBat_run$dat.covbat
      covf_df = cbind(other_info, df[bat_col], df[cov_col], covf_df)
      covf_df = covf_df[colnames(df)]
      write_csv(covf_df, paste0(argv$outdir))
    }
      message(sprintf("Results saved at %s", paste0(argv$outdir)))
    }else{
      if(argv$model == "lmer"){
        ComBat_run = predict(object = object, newdata = features, newbat = batch, newcovar = cbind(covariates, df[random]), eb = argv$eb, robust.LS = argv$robust.LS)
      }else{
        ComBat_run = predict(object = object, newdata = features, newbat = batch, newcovar = covariates, eb = argv$eb, robust.LS = argv$robust.LS)
      }
      comf_df = ComBat_run$dat.combat
      comf_df = cbind(other_info, df[bat_col], df[cov_col], comf_df)
      comf_df = comf_df[colnames(df)]
      write_csv(comf_df, paste0(argv$outdir))
      message(sprintf("Results saved at %s", paste0(argv$outdir)))
    }
  }else{
    reference_df[[bat_col]] = "reference"
    df_c = rbind(reference_df, new_data)
    df_c[[bat_col]] = as.factor(df_c[[bat_col]])
    if (is.null(covariates)){
      form_c = NULL
      combat_c = NULL
    }else{
      if(argv$model == "lmer"){
        form_c = df_c[cov_col]
        combat_c = cbind(df_c[cov_col], df_c[random])
      }else{
        form_c = df_c[cov_col]
        combat_c = df_c[cov_col]
      }
    }
    form = form_gen(x = argv$model, c = form_c, i = interaction, random = random, smooth = smooth, int_type = argv$int_type)
    if(family == "comfam"){
        ComBat_run = ComBatFamily::comfam(data = df_c[features_col],
                                          bat = df_c[[bat_col]], 
                                          covar = combat_c,
                                          model = eval(parse(text = argv$model)), 
                                          formula = as.formula(form), 
                                          ref.batch = "reference",
                                          ...)
      comf_df = ComBat_run$dat.combat[(nrow(reference_df)+1):nrow(df_c),]
      comf_df = cbind(new_data[-which(colnames(df) %in% used_col)], new_data[bat_col], new_data[cov_col], comf_df)
      comf_df = comf_df[colnames(df)]
      comf_df = rbind(untouched, comf_df)
      write_csv(comf_df, paste0(argv$outdir))
      message(sprintf("Results saved at %s", paste0(argv$outdir)))
    }else{
        ComBat_run = ComBatFamily::covfam(data = df_c[features_col],
                                          bat = df_c[[bat_col]] , 
                                          covar = combat_c,
                                          model = eval(parse(text = argv$model)), 
                                          formula = as.formula(form),
                                          ref.batch = "reference",
                                          ...)
      covf_df = ComBat_run$dat.covbat[(nrow(reference_df)+1):nrow(df_c),]
      covf_df = cbind(new_data[-which(colnames(df) %in% used_col)], new_data[bat_col], new_data[cov_col], covf_df)
      covf_df = covf_df[colnames(df)]
      covf_df = rbind(untouched, covf_df)
      write_csv(covf_df, paste0(argv$outdir))
      message(sprintf("Results saved at %s", paste0(argv$outdir)))
    }
    
  }
  if(!is.na(argv$mout)){
    message("Saving ComBat model......")
    saveRDS(ComBat_run, argv$mout)
    message(sprintf("ComBat model saved at %s", argv$mout))  
  }
}