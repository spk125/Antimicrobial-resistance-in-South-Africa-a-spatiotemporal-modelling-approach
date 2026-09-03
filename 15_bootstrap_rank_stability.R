#!/usr/bin/env Rscript

# Parametric-refit bootstrap stability of the prevalence ranks and the
# continuous surveillance-priority score displayed in Figure 4.
suppressPackageStartupMessages({
  library(DiceKriging)
  library(dplyr)
  library(writexl)
})

set.seed(20260827)
B <- 1000L
project <- "/Users/snehakotian/Documents/SouthAfrica_AMR_District Analysis"
model_root <- file.path(project, "data/models/scenario_crosswalkeddata")
joined_root <- file.path(project, "data/observed_stgpr_joined")
out_dir <- file.path(project, "outputs/bootstrap_figure4_priority_stability")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spec <- data.frame(
  slug = c("abaumanii_carbapenem_final_withcrosswalk_merged",
           "kpneumoniae_carbapenem_final_withcrosswalk_merged",
           "saureus_penicillinase_final_withcrosswalk_merged",
           "ecoli_3gc_final_withcrosswalk_merged",
           "paeruginosa_aminoglycoside_final_withcrosswalk_merged"),
  joined = c("abaumanii_carbapenem_observed_concat_stgpr.csv",
             "kpneumoniae_carbapenem_observed_concat_stgpr.csv",
             "saureus_penicillinase_observed_concat_stgpr.csv",
             "ecoli_3gc_observed_concat_stgpr.csv",
             "paeruginosa_aminoglycoside_observed_concat_stgpr.csv"),
  label = c("Carbapenem-resistant Acinetobacter baumannii",
            "Carbapenem-resistant Klebsiella pneumoniae",
            "Methicillin-resistant Staphylococcus aureus",
            "Third-generation cephalosporin-resistant Escherichia coli",
            "Aminoglycoside-resistant Pseudomonas aeruginosa"),
  stringsAsFactors = FALSE
)

annual_values <- function(joined, predicted) {
  joined$boot_pred <- predicted
  joined$value <- ifelse(!is.na(joined$resistance_per), joined$resistance_per,
                         joined$boot_pred)
  districts <- unique(joined$district_code[joined$year == 2023 &
                                             !is.na(joined$boot_pred)])
  out <- lapply(districts, function(d) {
    z <- joined[joined$district_code == d, ]
    starts <- which(!is.na(z$value) & abs(z$year - 2014) <= 1)
    ends <- which(!is.na(z$value) & abs(z$year - 2023) <= 1)
    if (!length(starts) || !length(ends)) return(NULL)
    si <- starts[which.min(abs(z$year[starts] - 2014))]
    ei <- ends[which.min(abs(z$year[ends] - 2023))]
    ny <- z$year[ei] - z$year[si]
    if (!is.finite(ny) || ny <= 0) return(NULL)
    data.frame(district_code=d, prevalence_2023=z$value[ei],
               annual_change=(z$value[ei]-z$value[si])/ny)
  })
  do.call(rbind, out)
}

run_one <- function(slug, joined_file, label, idx) {
  stgpr <- file.path(model_root, slug, "stgpr")
  stack <- file.path(model_root, slug, "stack")
  mod <- readRDS(file.path(stgpr, paste0(slug, "_stgpr_model.rds")))
  tr <- read.csv(file.path(stack, paste0(slug, "_meta_train.csv")), check.names=FALSE)
  mp <- read.csv(file.path(stack, paste0(slug, "_meta_pred.csv")), check.names=FALSE)
  pp <- read.csv(file.path(stgpr, paste0(slug, "_prediction_with_stgpr.csv")), check.names=FALSE)
  joined <- read.csv(file.path(joined_root, joined_file), check.names=FALSE)
  cols <- colnames(mod@X)
  X <- as.matrix(tr[,cols,drop=FALSE]); XP <- as.matrix(mp[,cols,drop=FALSE])
  map <- match(paste(joined$district_code,joined$year,sep="_"),
               paste(pp$district_code,pp$year,sep="_"))
  original_full <- rep(NA_real_,nrow(joined))
  original_full[!is.na(map)] <- pp$stgpr_pred[map[!is.na(map)]]
  original <- annual_values(joined,original_full)
  nr <- nrow(original)
  prevalence <- matrix(NA_real_,nr,B)
  change <- matrix(NA_real_,nr,B)
  set.seed(20260827 + idx*10000L)
  for (b in seq_len(B)) {
    ans <- tryCatch({
      ys <- as.numeric(simulate(mod,nsim=1,cond=FALSE,checkNames=FALSE))
      mb <- NULL
      invisible(capture.output(mb <- km(~.,design=X,response=ys,
        covtype=mod@covariance@name,control=list(parinit="ARD"))))
      pb <- as.numeric(predict(mb,newdata=XP,type="SK",se.compute=FALSE,
        cov.compute=FALSE,checkNames=FALSE)$mean)
      jf <- rep(NA_real_,nrow(joined)); jf[!is.na(map)] <- pb[map[!is.na(map)]]
      av <- annual_values(joined,jf)
      av[match(original$district_code,av$district_code),]
    }, error=function(e)e)
    if (!inherits(ans,"error")) {
      prevalence[,b] <- ans$prevalence_2023
      change[,b] <- ans$annual_change
    }
    if (b %% 100L == 0L) message(label, ": ", b, "/", B)
  }
  list(label=label,original=original,prevalence=prevalence,change=change)
}

results <- lapply(seq_len(nrow(spec)), function(i)
  run_one(spec$slug[i],spec$joined[i],spec$label[i],i))

# Figure 4 normalises all five displayed combinations against the same maxima.
orig_prev <- unlist(lapply(results,function(x)x$original$prevalence_2023))
orig_change <- unlist(lapply(results,function(x)x$original$annual_change))
max_prev <- max(orig_prev,na.rm=TRUE)
max_pos <- max(orig_change[orig_change>0],na.rm=TRUE)
max_neg <- max(abs(orig_change[orig_change<0]),na.rm=TRUE)
score_fun <- function(p,c,mp,mpos,mneg) p/mp*100 +
  ifelse(c>0,c/mpos*100,0) - ifelse(c<0,abs(c)/mneg*100,0)

all_prev <- do.call(rbind,lapply(results,`[[`,"prevalence"))
all_change <- do.call(rbind,lapply(results,`[[`,"change"))
ok <- which(colSums(is.finite(all_prev)&is.finite(all_change))==nrow(all_prev))
boot_max_prev <- apply(all_prev[,ok,drop=FALSE],2,max)
boot_max_pos <- apply(all_change[,ok,drop=FALSE],2,function(z)max(z[z>0]))
boot_max_neg <- apply(all_change[,ok,drop=FALSE],2,function(z)max(abs(z[z<0])))

district_tables <- vector("list",length(results))
summary_tables <- vector("list",length(results))
offset <- 0L
for (i in seq_along(results)) {
  x <- results[[i]]; nr <- nrow(x$original); rows <- offset+seq_len(nr); offset <- offset+nr
  pr <- x$prevalence[,ok,drop=FALSE]
  ch <- x$change[,ok,drop=FALSE]
  score <- sweep(pr,2,boot_max_prev,"/")*100 +
    sweep(pmax(ch,0),2,boot_max_pos,"/")*100 -
    sweep(abs(pmin(ch,0)),2,boot_max_neg,"/")*100
  score <- pmax(score,0)
  prevalence_ranks <- apply(pr,2,function(z)rank(-z,ties.method="average"))
  priority_ranks <- apply(score,2,function(z)rank(-z,ties.method="average"))
  original_score <- pmax(score_fun(x$original$prevalence_2023,
    x$original$annual_change,max_prev,max_pos,max_neg),0)
  original_prevalence_rank <- rank(-x$original$prevalence_2023)
  original_priority_rank <- rank(-original_score)
  prevalence_rho <- apply(prevalence_ranks,2,function(z)
    cor(original_prevalence_rank,z,method="spearman"))
  priority_rho <- apply(priority_ranks,2,function(z)
    cor(original_priority_rank,z,method="spearman"))
  district_tables[[i]] <- data.frame(
    pathogen_antibiotic=x$label,district_code=x$original$district_code,
    original_2023_prevalence=x$original$prevalence_2023,
    original_prevalence_rank=original_prevalence_rank,
    median_bootstrap_prevalence_rank=apply(prevalence_ranks,1,median),
    lower_95_bootstrap_prevalence_rank=apply(prevalence_ranks,1,quantile,.025,names=FALSE),
    upper_95_bootstrap_prevalence_rank=apply(prevalence_ranks,1,quantile,.975,names=FALSE),
    original_annual_change=x$original$annual_change,
    original_priority_score=original_score,
    original_priority_rank=original_priority_rank,
    median_bootstrap_priority_rank=apply(priority_ranks,1,median),
    lower_95_bootstrap_priority_rank=apply(priority_ranks,1,quantile,.025,names=FALSE),
    upper_95_bootstrap_priority_rank=apply(priority_ranks,1,quantile,.975,names=FALSE),
    stringsAsFactors=FALSE)
  summary_tables[[i]] <- data.frame(
    pathogen_antibiotic=x$label,districts=nr,successful_refits=length(ok),
    median_prevalence_rank_rho=median(prevalence_rho),
    prevalence_rank_rho_lower_95=quantile(prevalence_rho,.025,names=FALSE),
    prevalence_rank_rho_upper_95=quantile(prevalence_rho,.975,names=FALSE),
    median_priority_rank_rho=median(priority_rho),
    priority_rank_rho_lower_95=quantile(priority_rho,.025,names=FALSE),
    priority_rank_rho_upper_95=quantile(priority_rho,.975,names=FALSE))
}

district <- bind_rows(district_tables) %>%
  arrange(pathogen_antibiotic,original_priority_rank)
summary <- bind_rows(summary_tables)
write.csv(summary,file.path(out_dir,"figure4_rank_stability_summary.csv"),row.names=FALSE)
write.csv(district,file.path(out_dir,"figure4_district_rank_stability.csv"),row.names=FALSE)
write_xlsx(list(Summary=summary,District_results=district),
  file.path(out_dir,"figure4_priority_rank_stability.xlsx"))
print(summary)
