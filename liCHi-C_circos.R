# ========================= CIRCOS PLOTS ===========================

library(GenomicRanges)
library(GenomeInfoDb)
library(HiCaptuRe)
library(GenomicInteractions)
library(rtracklayer)   
library(dplyr)
library(circlize)


# ============================== SD annotations (hg38) =========================
sess <- browserSession("UCSC")
genome(sess) <- "hg38"
sd_tab_hg38 <- getTable(ucscTableQuery(sess, table = "genomicSuperDups"))

sd_gr_hg38 <- GRanges(
  seqnames = sd_tab_hg38$chrom,
  ranges   = IRanges(sd_tab_hg38$chromStart + 1, sd_tab_hg38$chromEnd)
)
seqlevelsStyle(sd_gr_hg38) <- "UCSC"

# FANK1 locus (hg38)
fank1_gr <- GRanges("chr10", IRanges(125896564, 126009592), strand = "+")

# ---- Building FANK1 SD family ----
tab <- sd_tab_hg38
a_start <- tab$chromStart + 1
a_end   <- tab$chromEnd
b_start <- tab$otherStart + 1
b_end   <- tab$otherEnd

caseA <- tab$chrom == as.character(seqnames(fank1_gr)) &
  a_end >= start(fank1_gr) & a_start <= end(fank1_gr)

caseB <- tab$otherChrom == as.character(seqnames(fank1_gr)) &
  b_end >= start(fank1_gr) & b_start <= end(fank1_gr)

sd_fank1 <- rbind(
  data.frame(chr = tab$otherChrom[caseA], start = b_start[caseA], end = b_end[caseA]),
  data.frame(chr = tab$chrom[caseB],     start = a_start[caseB], end = a_end[caseB])
)

sd_fank1_gr <- GRanges(sd_fank1$chr, IRanges(sd_fank1$start, sd_fank1$end), strand = "*")
seqlevelsStyle(sd_fank1_gr) <- "UCSC"

# ============================== R list with my input objects =================================
lichic_objects_sex <- readRDS("data_processed_liCHic/lichic_objects_sex.rds")

# ============================== Functions ===============================
to_ucsc_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

norm_chr <- function(gr) {
  old <- seqlevels(gr)
  new <- to_ucsc_chr(old)
  names(new) <- old
  seqlevels(gr) <- new
  gr
}

safe_name_spaces <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("[^A-Za-z0-9 ]+", "", x)
  trimws(x)
}

empty_links_df <- function() {
  data.frame(
    chr1 = character(), start1 = integer(), end1 = integer(),
    chr2 = character(), start2 = integer(), end2 = integer(),
    color = character(), stringsAsFactors = FALSE
  )
}

# ============================== Export links =========================
# Exports trans interactions where at least one anchor is on target_chr
# Colors legend is a = overlaps FANK1 SD family, b = overlaps any SD, c = no overlap
export_target_chr_links <- function(gi, file, target_chr = "Y", sd_fank1_gr, sd_gr_hg38) {
  
  target_chr <- to_ucsc_chr(target_chr)
  
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  keep_chr <- as.character(seqnames(a1)) == target_chr | as.character(seqnames(a2)) == target_chr
  gi <- gi[keep_chr]
  if (length(gi) == 0) {
    write.csv(empty_links_df(), file, row.names = FALSE)
    return(invisible(NULL))
  }
  
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  keep_trans <- as.character(seqnames(a1)) != as.character(seqnames(a2))
  gi <- gi[keep_trans]
  if (length(gi) == 0) {
    write.csv(empty_links_df(), file, row.names = FALSE)
    return(invisible(NULL))
  }
  
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  g1 <- GRanges(seqnames(a1), IRanges(start(a1), end(a1)))
  g2 <- GRanges(seqnames(a2), IRanges(start(a2), end(a2)))
  
  in_fank1 <- overlapsAny(g1, sd_fank1_gr, ignore.strand = TRUE) |
    overlapsAny(g2, sd_fank1_gr, ignore.strand = TRUE)
  
  in_sd <- overlapsAny(g1, sd_gr_hg38, ignore.strand = TRUE) |
    overlapsAny(g2, sd_gr_hg38, ignore.strand = TRUE)
  
  color <- ifelse(in_fank1, "a", ifelse(in_sd, "b", "c"))
  
  df <- data.frame(
    chr1   = as.character(seqnames(a1)),
    start1 = start(a1) - 1,  # 0-based for circos input
    end1   = end(a1),
    chr2   = as.character(seqnames(a2)),
    start2 = start(a2) - 1,
    end2   = end(a2),
    color  = color,
    stringsAsFactors = FALSE
  )
  
  write.csv(df, file, row.names = FALSE)
  invisible(df)
}

# Function to export trans interactions where one anchor overlaps the FANK1 locus
# Colors are  only for the partner anchor (opposite end of FANK1)
export_fank1_links <- function(
    gi, file,
    fank1_gr, sd_fank1_gr, sd_gr_hg38
) {
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  keep_fank1 <- overlapsAny(a1, fank1_gr, ignore.strand = TRUE) |
    overlapsAny(a2, fank1_gr, ignore.strand = TRUE)
  
  gi <- gi[keep_fank1]
  if (length(gi) == 0) {
    write.csv(data.frame(), file, row.names = FALSE)
    return(invisible(NULL))
  }
  
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  keep_trans <- as.character(seqnames(a1)) != as.character(seqnames(a2))
  gi <- gi[keep_trans]
  if (length(gi) == 0) {
    write.csv(data.frame(), file, row.names = FALSE)
    return(invisible(NULL))
  }
  
  a1 <- norm_chr(anchorOne(gi))
  a2 <- norm_chr(anchorTwo(gi))
  
  # Partner = opposite anchor of the one overlapping FANK1
  hit2 <- overlapsAny(a2, fank1_gr, ignore.strand = TRUE)
  partner <- a2
  partner[hit2] <- a1[hit2]
  
  partner_gr <- GRanges(seqnames(partner), IRanges(start(partner), end(partner)))
  
  in_fank1 <- overlapsAny(partner_gr, sd_fank1_gr, ignore.strand = TRUE)
  in_sd    <- overlapsAny(partner_gr, sd_gr_hg38,  ignore.strand = TRUE)
  color <- ifelse(in_fank1, "a", ifelse(in_sd, "b", "c"))
  
  df <- data.frame(
    chr1   = as.character(seqnames(a1)),
    start1 = start(a1) - 1,
    end1   = end(a1),
    chr2   = as.character(seqnames(a2)),
    start2 = start(a2) - 1,
    end2   = end(a2),
    color  = color,
    stringsAsFactors = FALSE
  )
  
  write.csv(df, file, row.names = FALSE)
  invisible(df)
}

# ============================== Regions for circos annotations ============================
# FANK1 SD family annotation
sd_gr <- sd_fank1_gr
mcols(sd_gr)$label <- paste0("FANK1_SD_", seq_along(sd_gr))
mcols(sd_gr)$col   <- "dodgerblue3"

# PARs (hg38 coordinates)
par_gr <- GRanges(
  seqnames = c("chrX","chrX","chrY","chrY"),
  ranges   = IRanges(
    start = c(10001, 155701383, 10001, 56887903),
    end   = c(2781479, 156030895, 2781479, 57217415)
  )
)
mcols(par_gr)$label <- c("PAR1","PAR2","PAR1","PAR2")
mcols(par_gr)$col   <- "goldenrod2"

# chrY regions of interest
yreg_gr <- GRanges(
  seqnames = c("chrY","chrY"),
  ranges   = IRanges(
    start = c(2785864, 56852859),
    end   = c(2937754, 56871834)
  )
)
mcols(yreg_gr)$label <- c("Y:2.8–2.9 Mb", "Y:56.85–56.87 Mb")
mcols(yreg_gr)$col   <- "firebrick3"

annot_gr <- c(par_gr, yreg_gr, sd_gr)

# Order: PAR first, chrY ROIs next, then SDs
prio <- ifelse(grepl("^PAR", mcols(annot_gr)$label), 1,
               ifelse(as.character(seqnames(annot_gr)) == "chrY", 2, 3))
annot_gr <- annot_gr[order(prio, seqnames(annot_gr), start(annot_gr))]

# ============================== Circos plot function ==================================
plot_circos <- function(
    csv_file,
    title = NULL,
    annot_gr,
    save_path = NULL,
    chromosomes = paste0("chr", c(1:22, "X", "Y")),
    link_lwd = 3,
    label_cex = 0.30
) {
  add_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  
  df <- tryCatch(read.csv(csv_file, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(df) > 0 && all(c("chr1","chr2") %in% names(df))) {
    df$chr1 <- add_chr(df$chr1)
    df$chr2 <- add_chr(df$chr2)
    df <- df[df$chr1 %in% chromosomes & df$chr2 %in% chromosomes, , drop = FALSE]
  } else {
    df <- data.frame()
  }
  
  # Link colors
  link_col <- NULL
  if (nrow(df) > 0) {
    pal <- c(
      a = adjustcolor("dodgerblue", 0.85),
      b = adjustcolor("red",       0.35),
      c = adjustcolor("grey20",    0.35)
    )
    if ("color" %in% names(df)) {
      link_col <- pal[as.character(df$color)]
      link_col[is.na(link_col)] <- pal["c"]
    } else {
      link_col <- rep(pal["c"], nrow(df))
    }
  }
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    svg(save_path, width = 7, height = 7)
    on.exit(dev.off(), add = TRUE)
  }
  
  circos.clear()
  circos.par(
    start.degree = 90, gap.degree = 2,
    cell.padding = c(0,0,0,0),
    track.margin = c(0.005, 0.005)
  )
  circos.initializeWithIdeogram(species = "hg38", chromosome.index = chromosomes)
  
  # Annotation labels
  annot_gr2 <- annot_gr[as.character(seqnames(annot_gr)) %in% chromosomes]
  if (length(annot_gr2) > 0) {
    annot_df <- data.frame(
      chr   = as.character(seqnames(annot_gr2)),
      start = start(annot_gr2),
      end   = end(annot_gr2),
      label = mcols(annot_gr2)$label,
      col   = if ("col" %in% colnames(mcols(annot_gr2))) mcols(annot_gr2)$col else "black",
      stringsAsFactors = FALSE
    )
    
    circos.par(track.height = 0.10)
    circos.genomicLabels(
      annot_df[, c("chr","start","end","label")],
      labels.column = 4,
      side = "inside",
      cex = label_cex,
      col = annot_df$col,
      line_col = annot_df$col
    )
  }
  
  # Links
  if (nrow(df) > 0) {
    for (i in seq_len(nrow(df))) {
      p1 <- (df$start1[i] + df$end1[i]) / 2
      p2 <- (df$start2[i] + df$end2[i]) / 2
      circos.link(df$chr1[i], p1, df$chr2[i], p2, col = link_col[i], border = NA, lwd = link_lwd)
    }
  }
  
  if (is.null(title)) title <- basename(csv_file)
  title(title)
  invisible(TRUE)
}

# ============================== Run on the chrY links ===============================
dir.create("Links_data_sex", showWarnings = FALSE, recursive = TRUE)

for (nm in names(lichic_objects_sex)) {
  gi  <- lichic_objects_sex[[nm]]
  out <- file.path("Links_data_sex", paste0("chrY_trans_links_", safe_name_spaces(nm), ".csv"))
  
  export_target_chr_links(
    gi,
    file         = out,
    target_chr   = "Y",
    sd_fank1_gr  = sd_fank1_gr,
    sd_gr_hg38   = sd_gr_hg38
  )
}

dir.create("Circos_Plot_Sex", showWarnings = FALSE, recursive = TRUE)
files <- list.files("Links_data_sex", pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  nm <- sub("\\.csv$", "", basename(f))
  plot_circos(
    csv_file  = f,
    title     = nm,
    annot_gr  = annot_gr,
    save_path = file.path("Circos_Plot_Sex", paste0(nm, "_circos.svg"))
  )
}

# ============================== Run on the FANK1 links ==============================
dir.create("Links_data_FANK1", showWarnings = FALSE, recursive = TRUE)

for (nm in names(lichic_objects_sex)) {
  gi  <- lichic_objects_sex[[nm]]
  out <- file.path("Links_data_FANK1", paste0("FANK1_trans_links_", safe_name_spaces(nm), ".csv"))
  
  export_fank1_links(
    gi,
    file         = out,
    fank1_gr     = fank1_gr,
    sd_fank1_gr  = sd_fank1_gr,
    sd_gr_hg38   = sd_gr_hg38
  )
}

dir.create("Circos_Plot_FANK1", showWarnings = FALSE, recursive = TRUE)
files_fank1 <- list.files("Links_data_FANK1", pattern = "\\.csv$", full.names = TRUE)

for (f in files_fank1) {
  nm <- sub("\\.csv$", "", basename(f))
  plot_circos(
    csv_file  = f,
    title     = nm,
    annot_gr  = annot_gr,
    save_path = file.path("Circos_Plot_FANK1", paste0(nm, "_circos.svg"))
  )
}
