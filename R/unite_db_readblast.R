#' Reads in a UNITE taxonomy table
#'
#' This function is adapted from \code{yingtools2::read.blastn.file} . It allows reading in a taxonomic annotation file blasted to UNITE.
#' The top hit is chosen first by evalue (lowest evalue), if several hits with same evalue, the Species with the highest
#' number of hits is chosen.
#'
#'
#'
#' @param tax.file The output file from blasting the ASV table to the UNITE database. Default is \code{"asv_seqs.fasta.unite.txt"}
#' @param tax_table logical, if \code{TRUE} (default), will return a data frame of taxonomy with top hit, ready to be converted to a tax_table, if \code{FALSE}, returns all hits
#' @return Taxonomy table data frame
#' @author Thierry Rolling
#' @export

read.blastn.unite <- function(tax.file="asv_seqs.fasta.unite.txt",tax_table=TRUE) {
  requireNamespace("data.table",quietly=TRUE)
  requireNamespace("tidyr",quietly = TRUE)
  requireNamespace("stringr",quietly=TRUE)
  requireNamespace("data.table",quietly = TRUE)
  #tax.file="uparse/total.5.repset.fasta.blastn.refseq_rna.txt";tax_table=TRUE;blastn.data=FALSE
  t <- data.table::fread(tax.file,quote="") %>% tbl_df()
  colnames(t) = c("qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send", "nident", "mismatch", "gapopen", "gaps", "ppos",
                  "frames", "qframe", "sframe", "qcovs", "qcovhsp", "evalue", "bitscore", "score", "length", "pident")
  t = t %>% tidyr::separate(sseqid, into=c("accession_unite","taxonomy","species_hypothesis"),sep="[|]") %>%
    tidyr::separate(taxonomy,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep="[;]",remove=F)
  t$Kingdom=stringr::str_remove(t$Kingdom,"k__")
  t$Phylum=stringr::str_remove(t$Phylum,"p__")
  t$Class=stringr::str_remove(t$Class,"c__")
  t$Order=stringr::str_remove(t$Order,"o__")
  t$Family=stringr::str_remove(t$Family,"f__")
  t$Genus=stringr::str_remove(t$Genus,"g__")
  t$Species=stringr::str_remove(t$Species,"s__")
  t$Species=stringr::str_replace(t$Species,"_"," ")
  t=t%>%
    mutate(otu=qseqid,
           otu.number=as.numeric(stringr::str_extract(otu,"(?<=ASV_)[0-9]+"))) %>%
    group_by(otu,Species)%>%
    mutate(occurence.sp=n())%>%
    ungroup()%>%
    mutate(genuslevel=grepl(" sp",Species))%>%
    group_by(otu) %>%
    arrange(evalue,genuslevel,desc(occurence.sp)) %>%
    filter(!duplicated(taxonomy)) %>%
    mutate(evalue.rank=dense_rank(evalue)) %>%
    select(otu,Phylum,Family,Species,evalue,evalue.rank,occurence.sp,pident,length,everything())
  if (!tax_table) {
    t <- t %>% ungroup() %>% arrange(otu.number)
  } else {
    t <- t %>%
      # mutate(n.ties=sum(dense_rank(evalue)==1),blast.data=paste0(Species," (eval=",evalue,",pid=",pident,")",collapse=";")) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      arrange(otu.number) %>%
      select(otu,evalue,pident,Species, Kingdom, Class, Order, Family, Genus,Phylum)
  }
  return(t)
}



#' Creates a color palette based on fungal taxonomy
#'
#' This function is adapted from \code{yingtools2::get.yt.palette} . It allows the use of the UNITE taxonomy (which differs from NCBI taxonomy).
#'
#'
#'
#' @param tax either a data.frame, phyloseq, or tax_table
#' @param parapsilosis logical, if \code{FALSE} (default), C. parapsilosis will only have one color, if \code{TRUE}, will color according to C. parapsilosis ASV
#' @author Thierry Rolling
#' @export

get.fungal.palette.unite <- function (tax, parapsilosis=FALSE) {
  #tax=phy.bmt.fungal
  requireNamespace("data.table",quietly=TRUE)
  requireNamespace("tidyr",quietly = TRUE)
  requireNamespace("stringr",quietly=TRUE)
  requireNamespace("data.table",quietly = TRUE)
  requireNamespace("yingtools2",quietly=TRUE)
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- yingtools2::get.tax(tax)
  }
  if (parapsilosis==TRUE){
    tax=tax%>% mutate(Species=if_else(Species=="Candida parapsilosis",otu,Species))
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!all(ranks %in% names(tax))) {
    stop("YTError: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  tax.dict <- tax[, ranks] %>% distinct()
  tax.dict$color <- rep(yingtools2::shades("gray", variation=0.25),length.out = nrow(tax.dict))
  basidio <- tax.dict$Phylum == "Basidiomycota"
  tax.dict$color[basidio] <- rep(yingtools2::shades("#C48C66", variation = 0.4), length.out = sum(basidio))
  malassezia <- tax.dict$Genus == "Malassezia"
  tax.dict$color[malassezia] <- rep(yingtools2::shades("#8A3030", variation = 0.6), length.out = sum(malassezia))
  molds <- c("Arthoniomycetes","Coniocybomycetes","Dothideomycetes",
             "Eurotiomycetes","Geoglossomycetes","Laboulbeniomycetes",
             "Lecanoromycetes","Leotiomycetes","Lichinomycetes","Orbiliomycetes",
             "Pezizomycetes","Sordariomycetes","Xylonomycetes")
  mold_group <- tax.dict$Class %in% molds
  tax.dict$color[mold_group] <- rep(yingtools2::shades("#ADDADA",variation=0.4),length.out=sum(mold_group))

  aspergillus <- tax.dict$Genus == "Aspergillus"
  tax.dict$color[aspergillus] <- rep(yingtools2::shades("#3F8D3D", variation = 0.4), length.out = sum(aspergillus))

  saccharo <- tax.dict$Order == "Saccharomycetales"
  tax.dict$color[saccharo] <- rep(yingtools2::shades("#F0C3C3", variation = 0.4), length.out = sum(saccharo))
  candida_group=c("Candida albicans","Candida tropicalis","Candida dubliniensis")
  candida <- tax.dict$Species %in%candida_group
  tax.dict$color[candida] <- rep(yingtools2::shades("#DE0000", ncolor=6,variation = 0.6), length.out = sum(candida))
  parapsilosis <- grepl("ASV",tax.dict$Species)
  tax.dict$color[parapsilosis] <- rep(yingtools2::shades("#F5990F", ncolor=6, variation = 0.6),length.out=sum(parapsilosis))
  parapsilosis1 <- grepl("parapsilosis|orthopsilosis|metapsilosis",tax.dict$Species)
  tax.dict$color[parapsilosis1] <- rep(yingtools2::shades("#F5990F", ncolor=6, variation = 0.6),length.out=sum(parapsilosis1))
  sc <- tax.dict$Species == "Saccharomyces cerevisiae"
  tax.dict$color[sc] <- rep("#756EAC")
  tax.palette <- structure(tax.dict$color, names = as.character(tax.dict$Species))
  tax.palette
}
