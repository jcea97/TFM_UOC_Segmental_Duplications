# ============================ Setup ===========================================

library(BSgenome.Hsapiens.NCBI.GRCh38)
library(HiCaptuRe)
library(GenomicInteractions)
library(HiContacts)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(data.table)
library(patchwork)
library(scales)


g <- BSgenome.Hsapiens.NCBI.GRCh38
ord <- c(as.character(1:22), "X", "Y", "M")

# ============================ Our Paths ===========================================
in_dir       <- "data_conf"
out_dir      <- "data_processed_liCHic"
fig_dir      <- file.path(out_dir, "figuras")
ibed_raw_dir <- file.path(out_dir, "ibed_raw")
ibed_sig_dir <- file.path(out_dir, "ibed_sig")

dir.create(out_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(ibed_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ibed_sig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================ Input original peakmatrix ==============================

# definitive_liCHiC_peakmatrix_replicates_B_cell_roadmap must be in memory first
data_filtered <- definitive_liCHiC_peakmatrix_replicates_B_cell_roadmap
rm(definitive_liCHiC_peakmatrix_replicates_B_cell_roadmap)

cols_common <- names(data_filtered)[1:11]
cols_keep <- c(
  "immtransB_fetal_1", "immtransB_fetal_2",
  "Mon_500K_1",        "Mon_500K_2",
  "nB_1M_3",           "nB_1M_4",
  "PC_WT_1",           "PC_WT_2",
  "memB_WT_3",         "memB_WT_1"
)
stopifnot(all(cols_keep %in% names(data_filtered)))

data_filtered <- dplyr::select(data_filtered, all_of(c(cols_common, cols_keep)))

pm_txt <- file.path(in_dir, "data_filtered_peakmatrix.txt")
write.table(data_filtered, pm_txt, sep = "\t", row.names = FALSE, quote = FALSE)

# ============================ Export the IBEDs ================================
pm <- load_interactions(
  file = pm_txt,
  genome = g,
  select_chr = c(as.character(1:22), "X", "Y", "MT")
)

interactions_list <- peakmatrix2list(pm)

purrr::iwalk(interactions_list, function(x, nm) {
  export_interactions(
    interactions = x,
    file         = file.path(ibed_raw_dir, paste0(nm, ".ibed")),
    format       = "ibed",
    over.write   = TRUE
  )
})

# ============================ Functions =========================================
get_ct <- function(gi, label) {
  chr1 <- sub("^chr", "", as.character(seqnames(anchorOne(gi))))
  chr2 <- sub("^chr", "", as.character(seqnames(anchorTwo(gi))))
  is_cis <- chr1 == chr2
  
  tibble(chr = c(chr1, chr2), is_cis = c(is_cis, is_cis)) |>
    group_by(chr) |>
    summarise(
      cis_n     = sum(is_cis, na.rm = TRUE),
      trans_n   = sum(!is_cis, na.rm = TRUE),
      total     = n(),
      trans_pct = trans_n / total,
      .groups   = "drop"
    ) |>
    mutate(group = label) |>
    arrange(factor(chr, levels = ord, ordered = TRUE))
}

make_plot_trans_by_chr <- function(
    ct_all,
    group_colors = c(Male = "#56B4E9", Female = "#E69F00"),
    main_title = NULL,
    save_path = NULL,
    save_width = 12, save_height = 6, save_units = "in", save_dpi = 300,
    y_limits = c(0, 1)
) {
  p <- ggplot(ct_all, aes(x = factor(chr, levels = ord), y = trans_pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.8)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    scale_y_continuous(labels = scales::percent, limits = y_limits, expand = c(0, 0)) +
    labs(x = "Chromosomes", y = "% of trans contacts", fill = "Group", title = main_title) +
    scale_fill_manual(values = group_colors)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = save_width, height = save_height, units = save_units, dpi = save_dpi)
  }
  p
}

make_plot_cistrans_panels <- function(
    gi_male, gi_female,
    titles = c("Male", "Female"),
    main_title = NULL,
    save_path = NULL,
    save_width = 10, save_height = 4.5, save_units = "in", save_dpi = 300
) {
  p_male   <- plotCisTrans(gi_male)   + ggtitle(titles[1])
  p_female <- plotCisTrans(gi_female) + ggtitle(titles[2])
  
  p <- (p_male | p_female) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = save_width, height = save_height, units = save_units, dpi = save_dpi)
  }
  p
}

make_plot_annotations_panels <- function(
    gi_male, gi_female,
    viewpoint = "B",
    titles = c("Male", "Female"),
    main_title = NULL,
    save_path = NULL,
    save_width = 10, save_height = 5.5, save_units = "in", save_dpi = 300
) {
  p_male   <- plotInteractionAnnotations(gi_male,   viewpoints = viewpoint, legend = TRUE) + ggtitle(titles[1])
  p_female <- plotInteractionAnnotations(gi_female, viewpoints = viewpoint, legend = TRUE) + ggtitle(titles[2])
  
  p <- (p_male | p_female) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = save_width, height = save_height, units = save_units, dpi = save_dpi)
  }
  p
}

# ============================ Cell type analysis function ============================
run_ct <- function(ct_name, peakmatrix_file, cs_male, cs_female) {
  
  message("==> ", ct_name)
  
  gi <- HiCaptuRe::load_interactions(
    file = peakmatrix_file,
    genome = g,
    select_chr = c(as.character(1:22), "X", "Y", "MT")
  )
  
  # Sex-specific interactions with CHiCAGO Score threshold
  gi_male   <- gi[mcols(gi)[[cs_male]]   > 5 & mcols(gi)[[cs_female]] < 5]
  gi_female <- gi[mcols(gi)[[cs_female]] > 5 & mcols(gi)[[cs_male]]   < 5]
  
  # Export filtered IBEDs
  export_interactions(gi_male,
                      file = file.path(ibed_sig_dir, paste0(ct_name, "_male_sig.ibed")),
                      format = "ibed", over.write = TRUE
  )
  export_interactions(gi_female,
                      file = file.path(ibed_sig_dir, paste0(ct_name, "_female_sig.ibed")),
                      format = "ibed", over.write = TRUE
  )
  
  # Summary table and plots
  ct_all <- bind_rows(get_ct(gi_male, "Male"), get_ct(gi_female, "Female"))
  
  p1 <- make_plot_trans_by_chr(
    ct_all,
    main_title = ct_name,
    save_path  = file.path(fig_dir, paste0("p1_", ct_name, ".png"))
  )
  p2 <- make_plot_cistrans_panels(
    gi_male, gi_female,
    main_title = ct_name,
    save_path  = file.path(fig_dir, paste0("p2_", ct_name, ".png"))
  )
  p3 <- make_plot_annotations_panels(
    gi_male, gi_female,
    viewpoint = "B",
    main_title = ct_name,
    save_path  = file.path(fig_dir, paste0("p3_", ct_name, ".png"))
  )
  
  list(
    gi_male   = gi_male,
    gi_female = gi_female,
    ct_table  = ct_all,
    plots     = list(p1 = p1, p2 = p2, p3 = p3)
  )
}

# ============================ Cell type peakmatrices ===============
cell_types <- list(
  immtransB_fetal = c("immtransB_fetal_1", "immtransB_fetal_2"),
  Mon_500K        = c("Mon_500K_1",        "Mon_500K_2"),
  nB_1M           = c("nB_1M_3",           "nB_1M_4"),
  PC_WT           = c("PC_WT_1",           "PC_WT_2"),
  memB_WT         = c("memB_WT_3",         "memB_WT_1")
)

purrr::iwalk(cell_types, function(cols, ct) {
  out_file <- file.path(out_dir, paste0(ct, "_peakmatrix.txt"))
  fwrite(data_filtered[, c(cols_common, cols), drop = FALSE], out_file, sep = "\t")
})

# ============================ Run all cell types ==============================
results <- list(
  immtransB = run_ct(
    ct_name         = "immtransB",
    peakmatrix_file = file.path(out_dir, "immtransB_fetal_peakmatrix.txt"),
    cs_male         = "CS_immtransB_fetal_1",
    cs_female       = "CS_immtransB_fetal_2"
  ),
  memB = run_ct(
    ct_name         = "memB",
    peakmatrix_file = file.path(out_dir, "memB_WT_peakmatrix.txt"),
    cs_male         = "CS_memB_WT_1",
    cs_female       = "CS_memB_WT_3"
  ),
  PC = run_ct(
    ct_name         = "PC",
    peakmatrix_file = file.path(out_dir, "PC_WT_peakmatrix.txt"),
    cs_male         = "CS_PC_WT_1",
    cs_female       = "CS_PC_WT_2"
  ),
  mon500 = run_ct(
    ct_name         = "mon500",
    peakmatrix_file = file.path(out_dir, "Mon_500K_peakmatrix.txt"),
    cs_male         = "CS_Mon_500K_2",
    cs_female       = "CS_Mon_500K_1"
  ),
  nB = run_ct(
    ct_name         = "nB",
    peakmatrix_file = file.path(out_dir, "nB_1M_peakmatrix.txt"),
    cs_male         = "CS_nB_1M_4",
    cs_female       = "CS_nB_1M_3"
  )
)

saveRDS(results, file.path(out_dir, "liCHiC_results_all.rds"))

# ============================ Figure % trans per chr ==========================
p1_all <- wrap_plots(lapply(results, function(x) x$plots$p1)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "% of trans per chromosome")

ggsave(file.path(fig_dir, "p1_all.png"), p1_all, width = 14, height = 8, dpi = 600)

# ============================ Load IBED objects (sex samples) ==================
ibed_files <- list.files(ibed_sig_dir, pattern = "\\.ibed$", full.names = TRUE)

lichic_objects_sex <- setNames(
  lapply(ibed_files, function(f) {
    load_interactions(file = f, genome = "BSgenome.Hsapiens.NCBI.GRCh38")
  }),
  tools::file_path_sans_ext(basename(ibed_files))
)

# Normalize names
names(lichic_objects_sex) <- sub("(.*_(female|male)).*", "\\1", names(lichic_objects_sex))

saveRDS(lichic_objects_sex, file.path(out_dir, "lichic_objects_sex.rds"))

# ============================ chrY top trans regions =========================
top_trans_Y <- function(gi, top_n = 20) {
  a1 <- anchorOne(gi); a2 <- anchorTwo(gi)
  chr1 <- sub("^chr", "", as.character(seqnames(a1)))
  chr2 <- sub("^chr", "", as.character(seqnames(a2)))
  
  keep <- (chr1 == "Y" | chr2 == "Y") & (chr1 != chr2)
  if (!any(keep)) return(data.frame())
  
  a1 <- a1[keep]; a2 <- a2[keep]
  chr1 <- chr1[keep]
  
  y_gr <- a1
  y_gr[chr1 != "Y"] <- a2[chr1 != "Y"]
  
  y_key <- paste0("Y:", start(y_gr), "-", end(y_gr))
  out <- sort(table(y_key), decreasing = TRUE)
  
  data.frame(
    region_Y = names(out)[seq_len(min(top_n, length(out)))],
    n_trans  = as.integer(out)[seq_len(min(top_n, length(out)))]
  )
}

topY_list <- lapply(lichic_objects_sex, top_trans_Y, top_n = 20)
saveRDS(topY_list, file.path(out_dir, "topY_list.rds"))

dfY <- bind_rows(lapply(names(topY_list), function(nm) {
  df <- topY_list[[nm]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  sex <- ifelse(grepl("_female$", nm), "Female",
                ifelse(grepl("_male$", nm), "Male", NA))
  cell_type <- sub("_(female|male)$", "", nm)
  
  df %>% mutate(sample = nm, sex = sex, cell_type = cell_type)
})) %>%
  filter(!is.na(sex))

# Extract coordinates from "Y:start-end"
parts <- do.call(rbind, strsplit(dfY$region_Y, "[:-]"))
dfY$start  <- as.numeric(parts[, 2])
dfY$end    <- as.numeric(parts[, 3])
dfY$mid_mb <- ((dfY$start + dfY$end) / 2) / 1e6

# Bin along chrY (Mb)
bin_size <- 0.1
dfY$bin_mb <- floor(dfY$mid_mb / bin_size) * bin_size

dfY_bin <- dfY %>%
  group_by(sex, cell_type, bin_mb) %>%
  summarise(n_trans = sum(n_trans), .groups = "drop")

write.table(dfY_bin, "dfY_bin_chrY_trans_contacts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(arrange(dfY, desc(n_trans)), "chrY_trans_contacts_ordered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

p_bin <- ggplot(dfY_bin, aes(x = cell_type, y = bin_mb, size = n_trans, colour = sex)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ sex, ncol = 1, scales = "free_y") +
  scale_size(range = c(2, 8)) +
  scale_colour_manual(values = c(Male = "#E69F00", Female = "#56B4E9")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title = "Trans contacts along chrY",
       x = "Cell type",
       y = "chrY position (Mb)",
       size = "Trans contacts",
       colour = "Sex")

# ============================ Viewpoints: partner regions =====================
vp1 <- GRanges("Y", IRanges(56852859, 56871834))
vp2 <- GRanges("Y", IRanges(2785864,  2937754))

partners_table <- function(gi, viewpoint, viewpoint_label, top_n = 50, score_col = "counts") {
  a1 <- anchorOne(gi); a2 <- anchorTwo(gi)
  
  h1 <- overlapsAny(a1, viewpoint)
  h2 <- overlapsAny(a2, viewpoint)
  keep <- h1 | h2
  if (!any(keep)) return(data.frame())
  
  a1k <- a1[keep]; a2k <- a2[keep]
  h2k <- h2[keep]
  
  # Partner = opposite anchor
  partner <- a2k
  partner[h2k] <- a1k[h2k]
  
  partner_key <- paste0(sub("^chr","", as.character(seqnames(partner))), ":",
                        start(partner), "-", end(partner))
  
  out <- tibble(partner_region = partner_key) %>%
    count(partner_region, name = "n_contacts") %>%
    mutate(viewpoint = viewpoint_label)
  
  sc <- tryCatch(mcols(gi)[[score_col]], error = function(e) NULL)
  if (!is.null(sc)) {
    sc <- sc[keep]
    sumtab <- tibble(partner_region = partner_key, w = sc) %>%
      group_by(partner_region) %>%
      summarise(sum_counts = sum(w, na.rm = TRUE), .groups = "drop")
    
    out <- left_join(out, sumtab, by = "partner_region")
  }
  
  out %>% arrange(desc(n_contacts)) %>% head(top_n)
}

get_sex <- function(nm) {
  if (grepl("female", nm, ignore.case = TRUE)) return("Female")
  if (grepl("male",   nm, ignore.case = TRUE)) return("Male")
  NA_character_
}

tab_all <- bind_rows(lapply(names(lichic_objects_sex), function(nm) {
  gi <- lichic_objects_sex[[nm]]
  
  bind_rows(
    partners_table(gi, vp1, "Y:56852859-56871834", top_n = 50),
    partners_table(gi, vp2, "Y:2785864-2937754",  top_n = 50)
  ) %>%
    mutate(sample = nm, sex = get_sex(nm))
}))

write.table(tab_all, "tabla_resumen_chrY.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

tab_vp_chr_ordered <- tab_all %>%
  extract(
    partner_region,
    into = c("partner_chr","partner_start","partner_end"),
    regex = "^([^:]+):(\\d+)-(\\d+)$",
    remove = FALSE
  ) %>%
  group_by(viewpoint, partner_chr, partner_region) %>%
  summarise(n_contacts_total = sum(n_contacts, na.rm = TRUE), .groups = "drop") %>%
  group_by(partner_chr) %>%
  mutate(chr_total = sum(n_contacts_total, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(chr_total), partner_chr, viewpoint, desc(n_contacts_total)) %>%
  select(viewpoint, partner_chr, partner_region, n_contacts_total)

write.table(tab_vp_chr_ordered, "tabla_resumen_chrY_ordered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

tab_top_region_per_chr <- tab_vp_chr_ordered %>%
  filter(partner_chr != "Y") %>%
  group_by(viewpoint, partner_chr) %>%
  slice_max(n_contacts_total, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(tab_top_region_per_chr, "tabla_resumen_chrY_top_regiones.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

dfp <- tab_top_region_per_chr %>%
  mutate(partner_chr = factor(partner_chr, levels = c(as.character(1:22), "X"))) %>%
  arrange(partner_chr, desc(n_contacts_total)) %>%
  mutate(partner_region = factor(partner_region, levels = unique(partner_region)))

p_top_partners <- ggplot(dfp, aes(x = n_contacts_total, y = partner_region, fill = viewpoint)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_contacts_total), hjust = -0.1, size = 3) +
  facet_wrap(~ viewpoint, scales = "free_y") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(x = "N contacts", y = "Top partner region (per chr)",
       title = "Most interactive partner region per chromosome") +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank()) +
  expand_limits(x = max(dfp$n_contacts_total) * 1.1)


