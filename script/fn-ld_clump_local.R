ld_clump_local =  function(out_dat, threshold) {
  tmp = subset(out_dat, pval.outcome < threshold) %>% TwoSampleMR::convert_outcome_to_exposure()
  snps = ieugwasr::ld_clump(
    dplyr::tibble(rsid = tmp$SNP,
                  pval = tmp$pval.exposure),
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 0.99,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = paste0(rdsf_personal, "data/1kg_eur/EUR")
  )
  return(subset(tmp,SNP%in%snps$rsid))
}