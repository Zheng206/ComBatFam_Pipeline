suppressMessages(library(argparser))
suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(ComBatFamQC))
suppressMessages(library(parallel))

## Read in arguments
p <- arg_parser("Visualizations of Lifespan age trends of neuroimaging-derived brain structures through shiny app.", hide.opts = FALSE)
p <- add_argument(p, "data", help = "path to the CSV or EXCEL file that contains data to be harmonized, covariates and batch information")
p <- add_argument(p, "--rois", short = '-r', help = "position of roi data(column numbers), eg: 1-5,9")
p <- add_argument(p, "--age", short = '-a', help = "column name of the age variable")
p <- add_argument(p, "--sex", help = "column name of the sex variable")
p <- add_argument(p, "--female", help = "female indicator, the value represent female in sex column.")
p <- add_argument(p, "--lowerquantile", short = '-l', help = "Specify a lower bound quantile. eg: 0.05, 0.25.", default = 0.25)
p <- add_argument(p, "--upperquantile", short = '-u', help = "Specify a upper bound quantile. eg: 0.75, 0.95.", default = 0.75)
p <- add_argument(p, "--mu", short = '-m', help = "An indicator of whether to smooth age variable, include it as a linear term or only include the intercept in the mu formula. smooth: y ~ pb(age), linear: y ~ age, default: y ~ 1.", default = "smooth")
p <- add_argument(p, "--sigma", help = "An indicator of whether to smooth age variable, include it as a linear term or only include the intercept in the sigma formula. smooth: ~ pb(age), linear: ~ age, default: ~ 1.", default = "smooth")
p <- add_argument(p, "--nu", short = '-n', help = "An indicator of whether to smooth age variable, include it as a linear term or only include the intercept in the nu formula. smooth: ~ pb(age), linear: ~ age, default: ~ 1.", default = "default")
p <- add_argument(p, "--tau", short = '-t', help = "An indicator of whether to smooth age variable, include it as a linear term or only include the intercept in the tau formula. smooth: ~ pb(age), linear: ~ age, default: ~ 1.", default = "default")
argv <- parse_args(p)

# Preprocess inputs
message('Checking inputs...')
if(is.na(argv$data)) stop("Missing input data") else {
  if(!grepl("csv$|xls$", argv$data)) stop("Input file must be a csv or an excel file") else {
    if(grepl("csv$", argv$data)) df = read.csv(argv$data) else df = read_excel(argv$data)
  }
}
df = data.frame(df)
if(is.na(argv$rois)) stop("Please identify the position of rois.") else {
  col = gsub("-",":",argv$rois)
  col_vec = eval(parse(text = paste0("c(", col, ")")))
}

features = colnames(df[col_vec])
age = argv$age
sex = argv$sex
df[[sex]] = as.factor(df[[sex]])
df = df %>% mutate(sex = case_when(sex == argv$female ~ "F",
                                   .default = "M"))
# Create sub_df for different features
df_var = paste0("sub_df_", 1:length(features))
for (i in 1:length(features)){
  sub_df = df[,c(features[i], age, sex)] %>% na.omit() 
  colnames(sub_df) = c("y", "age", "sex")
  assign(df_var[i], sub_df)
}

# Create age_list
age_list = mclapply(1:length(features), function(w){
  age_sub = age_list_gen (sub_df = eval(parse(text = paste0("sub_df_",w))),  lq = as.numeric(argv$lowerquantile), hq = as.numeric(argv$upperquantile), mu = argv$mu, sigma = argv$sigma, nu = argv$nu, tau = argv$tau)
  return(age_sub)
}, mc.cores = detectCores()) 
names(age_list) = features

quantile_type = c(paste0("quantile_", 100*as.numeric(argv$lowerquantile)), "median", paste0("quantile_", 100*as.numeric(argv$upperquantile)))

ComBatFamQC::age_shiny(age_list, features, quantile_type)
