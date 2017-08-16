## ----setup, results="hide", warning=FALSE, message=FALSE, cache=FALSE----
knitr::opts_chunk$set(error = TRUE, message = FALSE, warning = FALSE, cache = TRUE)
library(kmr4toxicogenetics)

savepath <- "results_large/" # to save data locally for use later
if (!dir.exists(savepath))
  dir.create(savepath)

## ----data----------------------------------------------------------------
# Load toxicity response - sizeof ncell x nchem
data("tox", package = "kmr4toxicogenetics")
cellnames <- rownames(tox)
ncelllines <- nrow(tox)
# Number of total cell lines
ncelllines
chemnames <- colnames(tox)
nchemicals <- ncol(tox)
# Number of total chemicals
nchemicals

# Load design matrix - sizeof ncell x p
data("design", package = "kmr4toxicogenetics")
# Cell line design matrix dimension
dim(design)
# Number of dim in each feature category
table(as.vector(
  sapply(colnames(design), function(a){
    if (!grepl("gram",a)) {
      return("covariate")
    } else {
      return(paste(strsplit(a,split="")[[1]][1:7],collapse=""))
    }
  })
))

# Load RNA-seq only cell line identifiers
data("id.rna", package = "kmr4toxicogenetics")
# Number of such cell lines
length(id.rna)

## ----kernel--------------------------------------------------------------
# Load cell line kernels
data("kcell", package = "kmr4toxicogenetics")
kcellname <- names(kcell)

# Load chemical kernels
data("kchem", package = "kmr4toxicogenetics")
kchemname <- names(kchem)

# Add mean kernel when we have several RBF kernels of different bandwidth for both cell lines and chemicals
kcell <- c(kcell, 
           list("KrnaseqMean" = do.call(WGCNA::pmean, kcell[grep('rnaseq', kcellname)]), 
                "KsnpMean" = do.call(WGCNA::pmean, kcell[grep('snp', kcellname)])))
kcellname <- c(kcellname, 'KrnaseqMean', 'KsnpMean')
kchem <- c(kchem, 
           list("KcdkMean" = do.call(WGCNA::pmean, kchem[grep('cdk', kchemname)]) ,
                "KpredtargetMean" = do.call(WGCNA::pmean, kchem[grep('predtarget', kchemname)]) ,
                "KsirmsMean" = do.call(WGCNA::pmean, kchem[grep('sirms', kchemname)])))
kchemname <- c(kchemname, 'KcdkMean', 'KpredtargetMean', 'KsirmsMean')

# Add multitask kernel for chemicals (interpolation between Dirac and uniform)
ninter <- 11
kchem <- c(kchem, lapply(as.list(seq(ninter)), function(i){
  kk <- ((i-1) * matrix(1, nchemicals, nchemicals) + (ninter - i) * diag(nchemicals)) / (ninter - 1)
  colnames(kk) <- chemnames
  rownames(kk) <- chemnames
  return(kk)
}))
kchemname <- c(kchemname, sapply(seq(ninter), function(i) paste('Kmultitask', i, sep = '')))

# Add empirical correlation kernel for chemicals
kchem <- c(kchem, 
           list("Kemp" = cor(tox)))
kchemname <- c(kchemname, 'Kempirical')

# Add mean kernel of heterogeneous sources for both cell lines and chemicals
kcell <- c(kcell, 
           list("Kint" = do.call(WGCNA::pmean, kcell[match(c("Kcovariates", "KrnaseqMean", "KsnpMean"), kcellname)])))
kcellname <- c(kcellname, 'Kint')
kchem <- c(kchem, 
           list("Kint" = do.call(WGCNA::pmean, kchem[match(c("KcdkMean", "KpredtargetMean", "KsirmsMean", "Kmultitask1", "Kmultitask11", "Kempirical"), kchemname)])))
kchemname <- c(kchemname, 'Kint')

# Preview kernel names
# Cell line kernels
kcellname
# Chemical kernels
kchemname

## ----pred----------------------------------------------------------------
# CV setup
nfolds <- 5
nrepeats <- 10

# Seed
seed <- 47

# Number of cores
mc.cores <- 1

# Main loop : Try each cell line kernel with each chemical kernel
sT = sF = matrix(NA, nrow = length(kcell), ncol = length(kchem), dimnames = list(kcellname, kchemname))
for (flagrna in c(TRUE, FALSE)) {
  # whether only RNA-seq cell lines vs all cell lines are used
  for (icell in seq(length(kcell))) {
    # for each cell line
    for (ichem in seq(length(kchem))) {
      # for each chemical
      objname <- paste("cvKMR", flagrna, kcellname[icell], kchemname[ichem], sep = "_")
      filename <- paste0(savepath, objname, ".RData")
      if (file.exists(filename)) {
        # CASE I - file exists from earlier runs
        # load the results
        cvres <- get(load(filename))
      } else {
        # CASE II - file not found
        message("running... keep rna cell ", flagrna, 
                " cell ", icell," out of ", length(kcell),
                " chem ", ichem," out of ", length(kchem),
                " ", objname)
        # specify kernel data and y data
        chemicalsKernel <- kchem[[ichem]]
        if(flagrna) {
          celllinesKernel <- kcell[[icell]][id.rna, id.rna]
          toxicity <- tox[id.rna, ]
        } else {
          celllinesKernel <- kcell[[icell]]
          toxicity <- tox
        }
        # run CV evaluation with KMR
        cvres <- evaluateCV(mypredictor = "predictorKMR", 
                            celllinesKernel = celllinesKernel, 
                            chemicalsKernel = chemicalsKernel, 
                            toxicity = toxicity, 
                            nfolds = nfolds, 
                            nrepeats = nrepeats, 
                            seed = seed, 
                            mc.cores = mc.cores)
        assign(objname, cvres)
        save(list = objname, file = filename)
      }
      rm(list = objname)
      # focus on CI scores and impute NA to 0.5 (random guess)
      cvres$matrix.ci[is.na(cvres$matrix.ci)] <- 0.5
      # average CI scores (over CV folds) per chemical per cell line
      if (flagrna)
        sT[icell,ichem] <- mean(cvres$matrix.ci)
      else
        sF[icell,ichem] <- mean(cvres$matrix.ci)
    }
  }
}

## ----perf----------------------------------------------------------------
# save default graphics options as to restore later
op <- par(no.readonly = TRUE)

## For all cell lines:
# CI heatmap
gplots::heatmap.2(sF, trace = "none", margin = c(8, 9), main = "CI")
# Mean CI over cell line kernels
par(mar = c(5.1,10.1,4.1,0.1))
barplot(sort(apply(sF, 1, mean)), horiz = TRUE, las = 1, xlim = c(min(sF) - 0.005, max(sF) + 0.005), xpd = FALSE, main = "Mean CI for cell line kernels")
par(op)
# Mean CI over chemical kernels
par(mar = c(5.1,8.1,4.1,0.1))
barplot(sort(apply(sF, 2, mean)), horiz = TRUE, las = 1, xlim = c(min(sF) - 0.005, max(sF) + 0.005), xpd = FALSE, main = "Mean CI for chemicals kernels", cex.names = 0.8)
par(op)
# Table of top scoring kernel pair
s <- reshape2::melt(sF, varnames = c("kcell", "kchem"))
s <- s[order(s$value, decreasing = TRUE), ]
s <- cbind("rank" = seq(nrow(s)), s)
rownames(s) <- s$rank
knitr::kable(head(s, 30), caption = "Top scoring kernel pair")

## For RNA-seq cell lines only:
# CI heatmap
gplots::heatmap.2(sT, trace = "none", margin = c(8, 9), main = "CI (RNA-seq only)")

## ----lambda--------------------------------------------------------------
# set kernels and data
celllinesKernel <- kcell[["Kint"]]
cat("Cell line kernel dim = ", dim(celllinesKernel))
chemicalsKernel <- kchem[["Kemp"]]
cat("Chemical kernel dim = ", dim(chemicalsKernel))
toxicity <- tox
cat("Toxicity response dim = ", dim(toxicity))

# train a model on full dataset
modelKMR <- kmr::cv.kmr(x = celllinesKernel, 
                        y = toxicity, 
                        kx_type = "precomputed", 
                        kt_type = "precomputed", 
                        kt_option = list(kt = chemicalsKernel), 
                        lambda = exp(-15:25), 
                        nfolds = 5, 
                        nrepeats = 1)

# plot regularization path
par(cex.axis = 1.5, cex.lab = 1.5, font.axis = 2, font.lab = 2)
plot(modelKMR)
par(op)
rm(modelKMR)

## ----stateart------------------------------------------------------------
# set parameter for running time comparison
nfolds.run <- 2
nrepeats.run <- 1
seed.run <- seed
mc.cores.run <- 1

# 1) ElasticNet
message("ElasticNet ...")
# pred perf
filename <- paste0(savepath, "cvElasticNet.RData")
if (file.exists(filename)) {
  load(filename)
} else {
  message("running... ", filename)
  cvElasticNet <- evaluateCV(mypredictor = "predictorElasticNet", celllines = design, toxicity = tox, nfolds = nfolds, nrepeats = nrepeats, seed = seed, mc.cores = mc.cores)
  save(cvElasticNet, file = filename)
}
# running time
ptElasticNet <- system.time(evaluateCV(mypredictor = "predictorElasticNet", celllines = design, toxicity = tox, nfolds = nfolds.run, nrepeats = nrepeats.run, seed = seed.run, mc.cores = mc.cores.run))

# 2) Lasso
message("Lasso ...")
# pred perf
filename <- paste0(savepath, "cvLasso.RData")
if (file.exists(filename)) {
  load(filename)
} else {
  message("running... ", filename)
  cvLasso <- evaluateCV(mypredictor = "predictorLasso", celllines = design, toxicity = tox, nfolds = nfolds, nrepeats = nrepeats, seed = seed, mc.cores = mc.cores)
  save(cvLasso, file = filename)
}
# running time
ptLasso <- system.time(evaluateCV(mypredictor = "predictorLasso", celllines = design, toxicity = tox, nfolds = nfolds.run, nrepeats = nrepeats.run, seed = seed.run, mc.cores = mc.cores.run))

# 3) RF
message("RF ...")
# pred perf
filename <- paste0(savepath, "cvRF.RData")
if (file.exists(filename)) {
  load(filename)
} else {
  message("running... ", filename)
  cvRF <- evaluateCV(mypredictor = "predictorRF", celllines = design, toxicity = tox, nfolds = nfolds, nrepeats = nrepeats, seed = seed, mc.cores = mc.cores)
  save(cvRF, file = filename)
}
# running time
ptRF <- system.time(evaluateCV(mypredictor = "predictorRF", celllines = design, toxicity = tox, nfolds = nfolds.run, nrepeats = nrepeats.run, seed = seed.run, mc.cores = mc.cores.run, ntree = 100))

# 4) KMR (with integrated kernel on cell lines and empirical kernel on chemicals)
message("KMR ...")
# pred perf
cvKMR <- get(load(paste0(savepath, "cvKMR_FALSE_Kint_Kempirical.RData")))
# running time
ptKMR <- system.time(evaluateCV(mypredictor = "predictorKMR", celllinesKernel = kcell[["Kint"]], chemicalsKernel = kchem[["Kemp"]], toxicity = tox, nfolds = nfolds.run, nrepeats = nrepeats.run, seed = seed.run, mc.cores = mc.cores.run))

## ----compare-------------------------------------------------------------
methodlist <- c("ElasticNet", "Lasso", "RF", "KMR")
methodlist <- ordered(methodlist, levels = methodlist)

# Running time per method
times <- lapply(methodlist, function(u) get(paste0("pt",u))[1])
times <- do.call("rbind",times)
rownames(times) <- methodlist
colnames(times) <- "Time (sec)"
knitr::kable(ceiling(times/2), caption = "Running time per method")

# Averaged CI across CV experiments per method per chemical
scores <- list()
for (i in seq_along(methodlist)) {
  cvres <- get(paste0("cv",methodlist[i]))
  # focus on CI and impute NA to 0.5 (random guess)
  cvres$matrix.ci[is.na(cvres$matrix.ci)] <- 0.5
  scores[[i]] <- data.frame(
    "method" = methodlist[i],
    "chemical" = colnames(cvres$matrix.ci),
    "CI" = colMeans(cvres$matrix.ci),
    stringsAsFactors = FALSE
  )
}
scores <- do.call("rbind", scores)
rownames(scores) <- seq(nrow(scores))

# Mean CI over chemicals per method
tabscores <- tapply(scores$CI, scores$method, function(u){
  c("Mean" = mean(u), "SD" = sd(u))
})
tabscores <- do.call("rbind", tabscores)
# preview
knitr::kable(tabscores, digits = 4, caption = "Mean CI over chemicals per method")

# Mean CI over chemicals per method - Signif test by two-sided t.test
pmatrix <- matrix(NA, nrow = length(methodlist), ncol = length(methodlist),
                  dimnames = list(methodlist, methodlist))
for (i in 1:(length(methodlist) - 1)) {
  for (j in (i + 1):length(methodlist)) {
    pmatrix[i,j] <- t.test(x = scores$CI[scores$method == methodlist[i]], 
                           y = scores$CI[scores$method == methodlist[j]], 
                           alternative = 'two.sided', mu = 0, paired = TRUE)$p.value
  }
}
# correct p-value for multiple testing with Benjamini-Hochberg
pmatrix.adj <- p.adjust(pmatrix, "BH")
attributes(pmatrix.adj) <- attributes(pmatrix)
# preview
pmatrix.adj

# Mean CI over chemicals per method - Boxplot
par(cex.axis = 1.5, cex.lab = 1.5, font.axis = 2, font.lab = 2)
boxplot(CI~method, data = scores, ylab = "CI", col = c("#f4a582", "#0571b0", "#92c5de", "#ca0020"), outline = FALSE)
par(op)

## ----session-info--------------------------------------------------------
sessionInfo()

