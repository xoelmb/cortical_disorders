# Script for plots related to the intersection of Mato Blanco et al results with Jourdon et al 2023 pairwise differential expression results
library(tidyverse)
library(cowplot)

# savedir <- 
# plotdir <-

# Load. human TF list from Lambert et al
TF.db.list <- read.csv(paste0(savedir, "DatabaseExtract_v_1.01.csv")) %>% filter(Is.TF. == "Yes") %>% pull(HGNC.symbol) %>% unique()

# Load DEG from Mato Blanco et al ASD vs Ctrl (merge UP in ASD and UP in Control tabs)
MatoBlanco.de <- bind_rows( readxl::read_xlsx(paste0(savedir, "markers_per_donor.nohighmito.with_filters_and_crossed.simple.xlsx"),sheet=1),
                         readxl::read_xlsx(paste0(savedir, "markers_per_donor.nohighmito.with_filters_and_crossed.simple.xlsx"),sheet=2))

# create 4 differential expresson (de) subsets of interest
MatoBlanco.de.autgenes <- MatoBlanco.de %>% filter(!in.sex_chr) %>% pull(gene) %>% unique
MatoBlanco.de.TF <- MatoBlanco.de %>% filter(gene.is.tf & !in.sex_chr) %>% pull(gene) %>% unique
MatoBlanco.de.earlyRG <- MatoBlanco.de %>% filter( !in.sex_chr, anno_cluster_fct == "RG early")  %>% pull(gene) %>% unique
MatoBlanco.de.earlyRG.TF <- MatoBlanco.de %>% filter( !in.sex_chr, gene.is.tf, anno_cluster_fct == "RG early")  %>% pull(gene) %>% unique

# Load pairwise results from Jourdon et al 
# Exlude sex genes (chrX or chrY genes) and only RGs cells: RG-hem, RG, RG-oRG, RG-tRG
Jourdon.de.pairwise <- readRDS(paste0(savedir, "DE.Jourdon.et.al._RGs.cell.types_pairwise.rds")) 


# Counts the number of pairs with significant differential expression in Jourdon et al dataset
de.jourdon.pairwise.counts <- 
  de.jourdon.pairwise %>% 
  group_by(stage,gene,cell.type) %>% summarize(n.test=n(), n.de= sum(adj_pval< 0.01 & abs(lfc) > .25), .groups="drop") %>% 
  group_by(stage,cell.type) %>% mutate(max.pairs.tested= max(n.test)) %>% ungroup() %>% # Compute the number of pairs tested per cell type and pair
  mutate(freq.deg = round(n.de / max.pairs.tested, 4)) %>% # Since the number of pairs tested varies, report the number as frequency of pairs tested
  filter(n.de > 0) # Keep only genes that are differentially expressed in at least 1 pair


# PLOTS
# Frequency/count Dotplot for main figure 5
de.jourdon.pairwise.counts %>% filter( gene %in% MatoBlanco.de.earlyRG.TF ) %>%
  # Arrange genes by values for RG at TD0
  mutate(n.de.order= ifelse(stage=="TD0"& cell.type=="RG",n.de,0)) %>% 
  group_by(gene) %>% mutate(val=sum(n.de.order)) %>% ungroup() %>%
  arrange(desc(val)) %>%
  
  ggplot(aes(y=factor(gene,unique(gene)), x= factor(cell.type, c("hem","RG","tRG","oRG"))))+
  geom_point(aes(size=n.de, color=freq.deg))+
  facet_grid(cols=vars(stage))+
  scale_color_viridis_c(option="A",direction=-1)+
  theme_minimal()+scale_x_discrete("",drop=F)+scale_y_discrete("")+
  scale_size_continuous(range=c(.5,5), breaks=c(1,3,6,12))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+coord_cartesian(clip="off")


# Cumulative plot for Main or Supplementary figure 
curve.order <- c("de.all","de.TFs","de.intersecting.all" , "de.intersecting.TF", "de.intersecting.RGearly", "de.intersecting.TF.in.RGearly")
curve.color <- c("de.all"="grey", "de.TFs"="black",  "de.intersecting.all"="red", "de.intersecting.TF"="green3", "de.intersecting.RGearly"="orange", "de.intersecting.TF.in.RGearly"="blue")
bind_rows(
  de.jourdon.pairwise.counts %>%   mutate(curve= "de.all"),
  de.jourdon.pairwise.counts %>% filter(gene %in% TF.db.list) %>% mutate(curve= "de.TFs"),
  de.jourdon.pairwise.counts %>% filter(gene %in% MatoBlanco.de.autgenes)  %>% mutate(curve= "de.intersecting.all"),
  de.jourdon.pairwise.counts %>% filter(gene %in% MatoBlanco.de.TF) %>% mutate(curve= "de.intersecting.TF"),
  de.jourdon.pairwise.counts %>% filter( gene %in% MatoBlanco.de.earlyRG ) %>% mutate(curve= "de.intersecting.RGearly"),
  de.jourdon.pairwise.counts %>% filter( gene %in% MatoBlanco.de.earlyRG.TF ) %>% mutate(curve= "de.intersecting.TF.in.RGearly")) %>% 
  filter( !(stage=="TD60" & cell.type =="hem")) %>% # Remove a case with too limiting results and cells for distribution analysis
  
  ggplot(aes( freq.deg, colour=factor(curve, curve.order)))+
  stat_ecdf(geom="step", size=1,alpha=.7)+
  facet_grid(rows=vars(factor(cell.type, c("hem","RG","tRG","oRG"))), cols=vars(stage), scales = "free")+
  theme_minimal()+ scale_x_reverse()+ scale_color_manual("")


# Heatmap to plot log2FC for TF of interest (for Supp figure)
Jourdon.de.pairwise %>% 
  filter(gene %in% MatoBlanco.de.earlyRG.TF ) %>%
  filter(adj_pval< 0.01 , abs(lfc) > .25  ) %>%
  group_by(gene) %>% mutate(n.de=n()) %>% ungroup() %>% arrange(desc(n.de))  %>%
  
  ggplot(aes(y= factor(gene, unique(gene)), x= factor(cell.type, c("hem","RG","tRG","oRG"))))+
  geom_tile(aes(fill = lfc), color="black",size=.1)+
  facet_grid( . ~  factor(pair.id, family.order ) + factor(stage)   , scales = "free", space="free")+
  theme_bw()+
  scale_color_manual(values=c("TRUE"="black", "FALSE"=NA),na.value=NA)+
  scale_x_discrete("",expand=c(0,0))+scale_y_discrete("",expand=c(0,0))+
  scale_fill_gradient2("logFC",high="red",low="blue", na.value = "lightgrey", limits = c(-5, 5), oob = scales::squish)+
  theme(axis.text.x = element_text(size=5,angle=90, hjust=1, vjust=0.5),
        strip.text.x = element_text(angle=90,hjust=0),
        panel.spacing = unit(0,"mm"),panel.grid = element_blank())+coord_cartesian(clip="off")


# Heatmap to plot logFC for significant differentially expressed genes in RG at TD0 for any gene set of interest
geneset.toplot <- readRDS(paste0(savedir,"gene.of.interest.rds")) %>% unique
family.order <- c("f7938","f07","fi03","f10530","f11175","f9230","fS8270","f10789","fU10999","f11251","fACE1575","fS1123") # Ordered by hierachical clustering 
Jourdon.de.pairwise %>% 
  filter(gene %in% geneset.toplot ) %>%
  filter(stage=="TD0",cell.type == "RG", adj_pval< 0.01 , abs(lfc) > .25  ) %>% 
  
  ggplot(aes( y= factor(gene, geneset.toplot), x= factor(pair.id, family.order) )) +
  geom_tile(aes(fill = lfc), color="black")+
  scale_color_manual(values=c("TRUE"="black", "FALSE"=NA),na.value=NA)+
  scale_x_discrete("",expand=c(0,0))+ scale_y_discrete("",expand=c(0,0))+
  scale_fill_gradient2("logFC",high="red",low="blue", 
                       na.value = "lightgrey",limits = c(-5, 5), oob = scales::squish)+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=0.5), axis.text.y = element_text(size=7),
        strip.text.x = element_text(angle=90,hjust=0), panel.grid = element_blank(), panel.border = element_rect(color="black",fill = NA))


# Venn Diagram for any intersection
library(ggvenn)
ggvenn(list("MatoBlanco.aut.deg"= MatoBlanco.de.autgenes , 
            "jourdon.aut.deg.inRG"= Jourdon.de.pairwise %>% filter( abs(lfc) > 0.25 & adj_pval < 0.01) %>% pull(gene) %>% unique(),
            "TF.db" = TF.db.list), show_percentage = F,text_size = 7)


# Compare frequency across pairs distribution between all degs and the set of interest
# Compute Kolmogorov-Smirnov test
imprinted.genes <- read.table(paste0(savedir , "imprinted_human_genes.tsv"),sep="\t", header=T) %>% filter(Status=="Imprinted") %>% pull(Gene) %>% unique()
stage.OI <- "TD0"
celltype.OI <- "RG"
ks.res <- ks.test(de.jourdon.pairwise.counts %>% filter(stage==stage.OI, cell.type==celltype.OI, gene %in%imprinted.genes) %>% pull(freq.deg) ,
                  de.jourdon.pairwise.counts  %>% filter(stage==stage.OI, cell.type==celltype.OI) %>% pull(freq.deg) , alternative="less")