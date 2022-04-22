library(optparse)
library(yaml)
library(Rcpp)

source('/app/SinkhornNNLSLinseedC.R')
source('/app/LinseedMetadata.R')
sourceCpp('/app/pipeline.cpp')


option_list <- list(
  make_option(c("-c", "--config"), action="store", default=NA, type='character',
              help="YAML with configutation path"),
  make_option(c("-i", "--init"), action="store", default=NA, type='integer',
              help="Initialization number")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$c)) {
  stop("YAML Configuration file path is not specified")
}

if (is.na(opt$i)) {
  opt$i <- 1
}

obj <- yaml.load_file(opt$c)
path_ <- file.path('/app', 'results', obj$path)
print(path_)
dir.create(path_, recursive = T, showWarnings = F)

analysis_name <- obj[['analysis_name']]
if (!is.na(opt$i)) {
  analysis_name <- paste0(analysis_name,"_",opt$i)
}

if (!is.null(obj[['coef_hinge_H']])) {
  coef_hinge_H <- as.numeric(obj[['coef_hinge_H']])
} else {
  coef_hinge_H <- 10
}

if (!is.null(obj[['coef_hinge_W']])) {
  coef_hinge_W <- as.numeric(obj[['coef_hinge_W']])
} else {
  coef_hinge_W <- 10
}

if (!is.null(obj[['coef_pos_D_h']])) {
  coef_pos_D_h <- as.numeric(obj[['coef_pos_D_h']])
} else {
  coef_pos_D_h <- 0.01
}

if (!is.null(obj[['coef_pos_D_w']])) {
  coef_pos_D_w <- as.numeric(obj[['coef_pos_D_w']])
} else {
  coef_pos_D_w <- 0.01
}

if (!is.null(obj[['coef_der_X']])) {
  coef_der_X <- as.numeric(obj[['coef_der_X']])
} else {
  coef_der_X <- 0.00001
}

if (!is.null(obj[['coef_der_Omega']])) {
  coef_der_Omega <- as.numeric(obj[['coef_der_Omega']])
} else {
  coef_der_Omega <- 0.0000001
}

if (!is.null(obj[['global_iterations']])) {
  global_iterations <- as.numeric(obj[['global_iterations']])
} else {
  global_iterations <- 5000
}

cat(paste("\n",analysis_name,":",
                  "\nNumber of iterations:",global_iterations,
                  "\nCell types:",obj[['cell_types']],
                  "\nDataset:",obj[['dataset']],
                  "\nPath:",path_,
                  "\nX derivative coefficient:",coef_der_X,
                  "\nOmega derivative coefficient:",coef_der_Omega,
                  "\nSum-to-one D_h coefficient:",coef_pos_D_h,
                  "\nSum-to-one D_w coefficient:",coef_pos_D_w,
                  "\nPositive proportions coefficient:",coef_hinge_H,
                  "\nPositive basis coefficient:",coef_hinge_W,"\n"))


if (!is.null(obj[['data']])) {
  data_ <- readRDS(obj[['data']][['path']])
  tmp_snkhrn <- SinkhornNNLSLinseed$new(dataset = obj[['dataset']], path = path_, analysis_name = analysis_name,
                                                cell_types = obj[['cell_types']],
                                                global_iterations = global_iterations,
                                                coef_der_X = coef_der_X, coef_der_Omega = coef_der_Omega,
                                                coef_hinge_H = coef_hinge_H, coef_hinge_W = coef_hinge_W,
                                                coef_pos_D_h = coef_pos_D_h, coef_pos_D_w = coef_pos_D_w,
                                                data = data_)
} else {
  tmp_snkhrn <- SinkhornNNLSLinseed$new(dataset = obj[['dataset']], path = path_, analysis_name = analysis_name,
                                                cell_types = obj[['cell_types']],
                                                global_iterations = global_iterations,
                                                coef_der_X = coef_der_X, coef_der_Omega = coef_der_Omega,
                                                coef_pos_D_h = coef_pos_D_h, coef_pos_D_w = coef_pos_D_w,
                                                coef_hinge_H = coef_hinge_H, coef_hinge_W = coef_hinge_W)
}
  

  tmp_snkhrn$selectTopGenes(obj[['top_genes']])

  if (!is.null(obj[['svd_k']])) {
    k <- max(obj[['cell_types']],obj[['svd_k']])
  } else {
    k <- obj[['cell_types']]
  }

  print("Scaling dataset")
  tmp_snkhrn$scaleDataset(iterations = 20)
  
  if (!is.null(obj[['projections']])) {
    print("Loading projections")
    tmp_snkhrn$readProjections(file = file.path(obj[['projections']][['path']],paste0("ct",tmp_snkhrn$cell_types),paste0(tmp_snkhrn$dataset,"_",opt$i,".rds")))
  } else {
    print("Get S and R projections")
    tmp_snkhrn$getSvdProjectionsNew(k=k)
  }


  if (!is.null(obj[['inits']])) {
    print("Loading initial points")
    tmp_snkhrn$readInitValues(file = file.path(obj[['inits']][['path']],paste0("ct",tmp_snkhrn$cell_types),paste0(tmp_snkhrn$dataset,"_",opt$i,".rds")))
  } else {
    print("Find initial points")
    tmp_snkhrn$selectInitOmega()
  }
  
  print("Run optimization")
  tmp_snkhrn$runOptimization()
  print("Save results")
  tmp_snkhrn$saveResults()
  print("Deconvolution complete")

  metadata_ <- LinseedMetadata$new(tmp_snkhrn)
  save(metadata_,file=paste0(metadata_$path_,"/",metadata_$analysis_name,".meta"))

