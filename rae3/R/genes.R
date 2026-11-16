#' Cancer genes for labelling copy-number peaks
#'
#' `rae_fetch_cancer_genes()` downloads the OncoKB cancer gene list and gene coordinates for one
#' genome build, and writes a small table to the user's cache. `rae_cancer_genes()` reads it back.
#' The list is fetched on demand and is **not** shipped with rae3, because its licence does not
#' allow redistribution.
#'
#' **Gene list.** OncoKB cancer gene list, public endpoint
#' <https://www.oncokb.org/api/v1/utils/cancerGeneList.txt> (no token needed). Please cite
#' Chakravarty et al. (2017) OncoKB: a precision oncology knowledge base. *JCO Precision Oncology*
#' 1:1-16, \doi{10.1200/PO.17.00011}. OncoKB is free for academic research use only; checking that
#' your use complies with the OncoKB terms (<https://www.oncokb.org/terms>) is your responsibility.
#' `oncogene` and `tsg` come from the OncoKB "Gene Type" column (`ONCOGENE`, `TSG`,
#' `ONCOGENE_AND_TSG`); genes without a type have both `FALSE`.
#'
#' **Coordinates.** NCBI RefSeq annotation as distributed by UCSC, without any account:
#' `ncbiRefSeqSelect` (one RefSeq Select / MANE transcript per gene) for hg38 and hg19, and
#' `hs1.ncbiRefSeq.gp.gz` (curated `NM_`/`NR_` transcripts only) for T2T-CHM13 (UCSC `hs1`).
#' Genes are matched by symbol, then by the RefSeq accession OncoKB gives for the transcript, then
#' by an unambiguous OncoKB alias. One interval per gene is the union of its transcripts on
#' chromosomes 1-22, X and Y (alternate, fix and random contigs are dropped). Genes that cannot be
#' placed are left out and counted in a message.
#'
#' @param build genome build: `"hg38"`, `"hg19"` or `"chm13"` (`"hs1"` and `"t2t"` are accepted as
#'   aliases of `"chm13"`).
#' @param dir cache directory. Defaults to the per-user cache given by
#'   `tools::R_user_dir("rae3", "cache")`.
#' @param overwrite download again even if the file already exists.
#' @return `rae_fetch_cancer_genes()`: the path of `<dir>/cancer_genes_<build>.tsv`, invisibly.
#'   `rae_cancer_genes()`: a data.frame with columns `gene`, `chrom` (`"1"`-`"22"`, `"X"`, `"Y"`,
#'   no `"chr"` prefix), `start` (1-based), `end`, `oncogene`, `tsg` (logical) and `source` (URLs
#'   and download date).
#' @examples
#' \dontrun{
#' rae_fetch_cancer_genes("hg38") # once; needs internet access
#' genes <- rae_cancer_genes("hg38")
#' genes[genes$gene == "TP53", ]
#' }
#' @name rae_cancer_genes
NULL

oncokb_url <- "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt"

gene_table_urls <- c(
  hg38 = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz",
  hg19 = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqSelect.txt.gz",
  chm13 = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gp.gz"
)

gene_build <- function(build) {
  build <- tolower(build[1])
  if (build %in% c("hs1", "t2t", "t2t-chm13", "chm13v2")) build <- "chm13"
  if (!build %in% names(gene_table_urls)) {
    stop("`build` must be one of \"hg38\", \"hg19\", \"chm13\" (or \"hs1\", \"t2t\").", call. = FALSE)
  }
  build
}

gene_cache_file <- function(build, dir) file.path(dir, paste0("cancer_genes_", build, ".tsv"))

download_or_stop <- function(url, dest) {
  old <- options(HTTPUserAgent = "R (rae3 package)")
  on.exit(options(old))
  ok <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(dest) || file.size(dest) == 0) {
    stop("Could not download ", url, ".\nCheck your internet connection (or proxy settings) ",
         "and try again.", call. = FALSE)
  }
  dest
}

#' @rdname rae_cancer_genes
#' @export
rae_fetch_cancer_genes <- function(build = c("hg38", "hg19", "chm13"),
                                   dir = tools::R_user_dir("rae3", "cache"),
                                   overwrite = FALSE) {
  build <- gene_build(build)
  out <- gene_cache_file(build, dir)
  if (file.exists(out) && !overwrite) {
    message("Using cached ", out, " (set overwrite = TRUE to download again).")
    return(invisible(out))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("rae3_genes_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  ## OncoKB cancer gene list
  okb <- utils::read.delim(download_or_stop(oncokb_url, file.path(tmp, "oncokb.txt")),
                           check.names = FALSE, quote = "", colClasses = "character",
                           na.strings = character())
  need <- c("Hugo Symbol", "Gene Type", "GRCh38 RefSeq", "GRCh37 RefSeq", "Gene Aliases")
  if (!all(need %in% names(okb))) {
    stop("Unexpected OncoKB file format; missing column(s): ",
         paste(setdiff(need, names(okb)), collapse = ", "), call. = FALSE)
  }
  okb <- okb[nzchar(okb[["Hugo Symbol"]]), ]
  okb <- okb[!duplicated(okb[["Hugo Symbol"]]), ]

  ## Gene coordinates (UCSC genePred; hg38/hg19 tables carry a leading `bin` column)
  url <- gene_table_urls[[build]]
  gp <- utils::read.delim(download_or_stop(url, file.path(tmp, basename(url))), header = FALSE,
                          quote = "", colClasses = "character", comment.char = "")
  if (ncol(gp) == 16) gp <- gp[-1]
  gp <- data.frame(acc = sub("\\.[0-9]+$", "", gp[[1]]), chrom = gp[[2]],
                   start = as.numeric(gp[[4]]) + 1, end = as.numeric(gp[[5]]), sym = gp[[12]])
  gp$chrom <- sub("^chr", "", gp$chrom)
  gp <- gp[gp$chrom %in% c(1:22, "X", "Y") & grepl("^N[MR]_", gp$acc), ]

  ## Match: symbol, then OncoKB RefSeq accession, then unambiguous alias
  syms <- okb[["Hugo Symbol"]]
  hit <- ifelse(syms %in% gp$sym, syms, NA_character_)
  acc_sym <- stats::setNames(gp$sym, gp$acc)
  for (col in c("GRCh38 RefSeq", "GRCh37 RefSeq")) {
    acc <- sub("\\.[0-9]+$", "", okb[[col]])
    miss <- is.na(hit) & acc %in% names(acc_sym)
    hit[miss] <- acc_sym[acc[miss]]
  }
  aliases <- strsplit(okb[["Gene Aliases"]], ",\\s*")
  for (i in which(is.na(hit))) {
    cand <- unique(intersect(setdiff(aliases[[i]], syms), gp$sym))
    if (length(cand) == 1) hit[i] <- cand
  }
  hit[duplicated(hit) & !is.na(hit)] <- NA # never place two OncoKB genes on one locus

  ## One interval per gene: union of transcripts on the chromosome holding most of them
  gp <- gp[gp$sym %in% hit, ]
  main_chr <- tapply(gp$chrom, gp$sym, function(x) names(which.max(table(x))))
  gp <- gp[gp$chrom == main_chr[gp$sym], ]
  lo <- tapply(gp$start, gp$sym, min)
  hi <- tapply(gp$end, gp$sym, max)
  keep <- !is.na(hit)
  type <- okb[["Gene Type"]][keep]
  res <- data.frame(gene = syms[keep], chrom = unname(main_chr[hit[keep]]),
                    start = unname(lo[hit[keep]]), end = unname(hi[hit[keep]]),
                    oncogene = type %in% c("ONCOGENE", "ONCOGENE_AND_TSG"),
                    tsg = type %in% c("TSG", "ONCOGENE_AND_TSG"),
                    source = paste0("OncoKB ", oncokb_url, "; coordinates ", url, "; downloaded ",
                                    format(Sys.Date())))
  chr_rank <- match(res$chrom, c(1:22, "X", "Y"))
  res <- res[order(chr_rank, res$start), ]
  n_miss <- sum(!keep)
  if (n_miss > 0) {
    message(n_miss, " of ", length(syms), " OncoKB genes could not be placed on ", build,
            " and were left out: ", paste(utils::head(syms[!keep], 10), collapse = ", "),
            if (n_miss > 10) ", ...")
  }
  utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", nrow(res), " genes to ", out)
  invisible(out)
}

#' @rdname rae_cancer_genes
#' @export
rae_cancer_genes <- function(build = "hg38", dir = tools::R_user_dir("rae3", "cache")) {
  build <- gene_build(build)
  f <- gene_cache_file(build, dir)
  if (!file.exists(f)) {
    stop("No cancer gene table for ", build, " in ", dir, ".\n",
         "Run rae_fetch_cancer_genes(\"", build, "\") once (needs internet access).\n",
         "The list is the OncoKB cancer gene list (Chakravarty et al. 2017, JCO Precision ",
         "Oncology, doi:10.1200/PO.17.00011); it is not redistributed with rae3 and is free for ",
         "academic use under the OncoKB terms.", call. = FALSE)
  }
  utils::read.delim(f, colClasses = c("character", "character", "numeric", "numeric",
                                      "logical", "logical", "character"), quote = "")
}
