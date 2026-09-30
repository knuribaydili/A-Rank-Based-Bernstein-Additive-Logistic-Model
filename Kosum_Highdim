
rm(list = ls()); gc()
Sys.setenv(BRCM_AUTORUN = "FALSE")
stopifnot(file.exists("BRCM_simulation.R"))
source("BRCM_simulation.R")

cat("\n=== self-test ===\n")
stopifnot(brcm_self_test(verbose = TRUE))

t0 <- Sys.time()
hd <- highdim_analysis(
  num_sim   = 500,
  ps        = c(10, 20, 50),
  ns        = c(500, 1000),
  top_feats = 8L,
  verbose   = TRUE
)
cat(sprintf("\nGecen sure: %.1f saat\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))

saveRDS(hd, "02_highdim_final.rds")
write.csv(hd, "T_highdim_screening_final.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo_highdim_final.txt")

cat("\n=== dogrulama ===\n")
stopifnot(nrow(hd) == 18L)
stopifnot(setequal(unique(hd$p), c(10, 20, 50)))
stopifnot(setequal(unique(hd$n), c(500, 1000)))
stopifnot(all(hd$Attempted == 500))
stopifnot(all(hd$n_possible_pairs == hd$p * (hd$p - 1) / 2))
stopifnot(all(hd$FWER_false_selection_pct >= 0 &
              hd$FWER_false_selection_pct <= 100))

stopifnot(all(hd$FWER_false_selection_pct / 100 <=
              hd$Mean_false_surfaces + 0.002))

z <- grepl("null", hd$Process, ignore.case = TRUE)
stopifnot(all(hd$Pct_any_selection[z] == hd$FWER_false_selection_pct[z]))

cat("  butun kontroller GECTI\n\n")

print(hd[, c("p", "n", "Process", "Attempted", "Successful",
             "n_possible_pairs", "Mean_candidate_union_pairs",
             "Correct_pair_pct", "Mean_n_surfaces",
             "Mean_false_surfaces", "FWER_false_selection_pct")])

draw_screening_figure(hd, "F3_screening_by_dimension.png")

cat("\n============================================================\n")
cat("  TAMAMLANDI. Gonderilecek dosya:\n")
cat("    T_highdim_screening_final.csv\n")
cat("    F3_screening_by_dimension.png  (Figure 3)\n")
cat("============================================================\n")
