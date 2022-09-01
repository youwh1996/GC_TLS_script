GC[["percent.mt"]] <- PercentageFeatureSet(GC, pattern ="^MT-")
plot1 <- FeatureScatter(GC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
GC <- subset(GC, subset = nFeature_RNA > 200 & nFeature_RNA < 3000  & percent.mt<25 )


GC <- NormalizeData(GC, normalization.method = "LogNormalize", scale.factor = 10000)
GC <- FindVariableFeatures(GC, selection.method = "vst", nfeatures = 3000)
GC <- ScaleData(GC, features = rownames(GC))
GC <- RunPCA(GC, features = VariableFeatures(object = GC))
RunHarmony(GC,"orig.ident")->GC
GC <- FindNeighbors(GC, dims = 1:30,reduction = "harmony")
GC <- FindClusters(GC, resolution = 0.3)
GC <- RunUMAP(GC, dims = 1:30,reduction = "harmony")

DimPlot(GC,reduction='umap',label=T)
subset(GC,idents=c(0:6))->Tcell

library(foreach)
library(pheatmap)
library(metacell)
library(tgconfig)
library(tgstat)
#
set_param("mc_plot_device",'pdf', "metacell")
set_param("scm_spike_regexp","^ERCC-","metacell")
set_param("scm_mc_mark_k_per_clust",100,"metacell") #default: 5
set_param("scm_mc_mark_min_gene_cov",0.3,"metacell") # default: 0.25
set_param("scm_mc_mark_min_gene_fold",2,"metacell") # default: 1.5

set_param("mcell_mc2d_K",30,"metacell") # default: 20
set_param("mcell_mc2d_T_edge",0.02,"metacell") # default: 0.05
set_param("mcell_mc2d_max_confu_deg",4,"metacell") # default: 5
set_param("mcell_mc2d_edge_asym",FALSE,"metacell") # default: TRUE
set_param("mcell_mc2d_proj_blur",0.02,"metacell") # default: 0.02
grDevices::pdf.options(useDingbats = FALSE)


if(!dir.exists("GCdb")) dir.create("GCdb/")
scdb_init("GCdb/", force_reinit=T)
if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="GZMB")

scdb_add_mat('GC',mat)


mcell_plot_umis_per_cell("GC")

####
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat)) 
nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
pre_nr_term <- c("^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("MALAT1", "XIST", "XIST_intron")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
pre_bad_genes
##
genes_RBC <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
cells_RBC <- names(which((apply(mat@mat[intersect(genes_RBC,rownames(mat@mat)),],2,sum))>=1))
mcell_mat_ignore_genes(new_mat_id='GC', mat_id='GC', pre_bad_genes, reverse=F)
mcell_mat_ignore_cells(new_mat_id='GC', mat_id='GC', ig_cells = genes_RBC, reverse = F)

gene_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','MX1','RSAD2','TOP2A','MKI67','STMN1')
tab_fn = "./lateral_gmods.txt"
mcell_mat_rpt_cor_anchors(mat_id='GC', gene_anchors = gene_anchors, cor_thresh = 0.1,
                          gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table('./lateral_gmods.txt', header=T)
foc_genes = apply(gcor_mat[, intersect(colnames(gcor_mat),gene_anchors)], 1, which.max)
#
mat = scdb_mat("GC")
mcell_add_gene_stat(gstat_id="GC", mat_id="GC", force=T)
mcell_gset_filter_varmean(gset_id="GC_feats", gstat_id="GC", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "GC_feats", gstat_id="GC", T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id="GC", gset_id="GC_feats")
gset <- scdb_gset("GC_feats")
pst_genes <- names(gset@gene_set)
pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                 "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                 "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
pst_ex_genes <- c()
pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
pst_add_genes <- c()
final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
final_genes


gset = gset_new_gset(sets = gset@gene_set[final_genes], desc = "final genes")
scdb_add_gset("GC_feats", gset)

###一般K为总细胞数的平方根的数目，若细胞数很少，K=20-40

mcell_add_cgraph_from_mat_bknn(mat_id="GC",
                gset_id = "GC_feats",
                graph_id="GC_graph",
                K=277,
                dsamp=T)


mcell_coclust_from_graph_resamp(
                coc_id="GC_coc1000",
                graph_id="GC_graph",
                min_mc_size=90,
                p_resamp=0.75, n_resamp=1000)



mcell_mc_from_coclust_balanced(
                coc_id="GC_coc1000",
                mat_id= "GC",
                mc_id= "GC_mc",
                K=100, min_mc_size=90, alpha=2) 


              

mc = scdb_mc("GC_mc")
mc@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc@mc_fp))
scdb_add_mc("GC_mc",mc)
mcell_mc2d_force_knn(mc2d_id="GC_2dproj",mc_id="GC_mc", graph_id="GC_graph")
mcell_mc2d_plot(mc2d_id="GC_2dproj")#

lfp <- round(log2(mc@mc_fp),2)
write.table(lfp,file='./lfp.txt',sep='\t',quote=F)

mcell_mc2d_plot_by_factor(
  mc2d_id='GC_2dproj',
  mat_id='GC',
  meta_field='stage',
  meta_data_vals = NULL,
  single_plot =T,
  filter_values = NULL,
  filter_name = NULL,
  ncols = NULL,
  neto_points = F,
  colors = mc@colors
)
##

lat_genes=names(scdb_gset("GC_feats")@gene_set)
T_genes = intersect(names(which(apply(abs(lfp),1,max)>log2(1.5))),lat_genes)

mat_b_ds = scm_downsamp(mat@mat, 500)
library(tgstat)
T_cor_c = tgs_cor(t(as.matrix(mat_b_ds[T_genes, ])), spearman=T)
diag(T_cor_c) = NA
T_cor_mc=cor(t(lfp[T_genes,]))
diag(T_cor_mc) = NA
blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
pdf(file='./figs/submc_T_cor.pdf', width=max(700, 300 + length(T_genes) * 12)/72, height=max(700, 300 + length(T_genes) * 12)/72)
pheatmap(pmin(pmax(T_cor_mc, -0.7), 0.7), clustering_method="ward.D2", cutree_rows=15, cutree_cols=15, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols)
dev.off()

e_gc = mc@e_gc[T_genes, ] * 1000
pdf(file='./figs/submc_T_exp.pdf', width=max(700, 300 + length(T_genes) * 12)/72, height=max(700, 300 + length(T_genes) * 12)/72)
pheatmap(e_gc,clustering_method='ward.D2',border_color=NA,col=blwtrd_cols,scale='row',cellwidth=10, cellheight=10,fontsize=12,breaks=unique(c(seq(-5,5, length=101))))
dev.off()



mc@colors->color
names(color)=mc@annots
color[unique(mc@annots)]->color
pdf("./TCELL_annot.pdf",useDingbats=F,width=8,height=6)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=2)+
scale_fill_manual(values=color)
dev.off()


 mc=scdb_mc("mc_annot")
 mc2d=scdb_mc2d("GC_2dproj")
 sc_info=data.frame(sc_x=mc2d@sc_x,sc_y=mc2d@sc_y,row.names=names(mc2d@sc_y))
sc_info$type=factor(mc@mc)
levels(sc_info$type)=mc@annots


pdf("./2dproject.pdf",useDingbats=F,width=10,height=7)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=2)+
scale_fill_manual(values=c('CD4-CXCL13-TNFRSF18'='#44CEE0','CD4-Treg'='#6591CF','Th17-like'='#7FBE8C','CD4-FOS'='#FF0095','CD4-CXCL13-TCF7'='#FFD600','CD4-Naive'='#4E6D25','CD8-CTL'='#BA6FB5','CD8-TRM-IFNG'='#3551AA','CD8-Cycling'='#B41322','CD8-MT2A'='#F89373','CD8-CXCL13'='#215050','CD8-Naive'='#A1DBEA','CD8-ISG15'='#FF2D3C','CD8-IL17A'='#915D63','CD8-TRM-XCL2'='#96362B','CD4-MAIT'='#E26500'))
dev.off()

minnumbcells <- min(sapply(c(as.character(unique(sc_info$S))), function(x) length(rownames(sc_info[sc_info$S==x,]))))
set.seed(22)
plotcells <-  as.vector(sapply(c(as.character(unique(sc_info$S))), function(x) sample(rownames(sc_info[sc_info$S==x,]), size=minnumbcells, replace=F, prob=NULL )))


pdf('./figs/density_2d_S.pdf',useDingbats=F)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$S=='pos'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$S=='pos'),],aes(x=sc_x,y=sc_y),color="black")


ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$S=='nes'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$S=='nes'),],aes(x=sc_x,y=sc_y),color="black")
dev.off()

df = data.frame(type =mc@annots[ord_by_id[["mc_annot"]]])
ha = HeatmapAnnotation(df = df,col = list(type =c('outlier'='white','CD4-CXCL13-TNFRSF18'='#44CEE0','CD4-Treg'='#6591CF','Th17-like'='#7FBE8C','CD4-FOS'='#FF0095','CD4-CXCL13-TCF7'='#FFD600','CD4-Naive'='#4E6D25','CD8-CTL'='#BA6FB5','CD8-TRM-IFNG'='#3551AA','Cycling-CD8'='#B41322','CD8-MT2A'='#F89373','CD8-CXCL13'='#215050','CD8-Naive'='#A1DBEA','CD8-ISG15'='#FF2D3C','CD8-IL17A'='#915D63','CD8-TRM-XCL2'='#96362B','CD4-MAIT'='#E26500')))
Heatmap(mat2[,ord_by_id[["mc_annot"]]],top_annotation = ha, name = "Rel. E",cluster_rows = F, cluster_columns = F, show_row_names = T, column_names_gp = gpar(fontsize = fontsize))


res_scDC_noClust1 <- scDC::scDC_noClustering(as.character(sc_info$type), as.character(sc_info$TLS), calCI = TRUE, 
                                     calCI_method = c("BCa"),
                                     nboot = 1000,ncores=10)
res1=res_scDC_noClust1$result
res1$mean=apply(res_scDC_noClust1$thetastar,1,mean)
res1<-data.frame(type=as.character(res1$cellTypes[1:16]),high=res1$mean[1:16],low=res1$mean[17:32])
row.names(res1)=res1$type
res1$type<-NULL
res1$fc<-c()
res1$fc=log2(res1$high/res1$low)

for(i in 1:nrow(res1)){
    if(res1$fc[i]>1){
        res1$reg[i]='UP'
    }else if(res1$fc[i]<(-1)){
         res1$reg[i]='DOWN'
    }else{
        res1$reg[i]='NAN'
    }
}
ggplot(res1,aes(x=group,y=type))+
geom_point(aes(size=abs(fc),color=reg))+
scale_color_manual(values=c('UP'='#ED2A28','DOWN'='#2C9685','NAN'='#B5B6B8'))+
theme_bw()+
theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


plotgeneexp(c('CD4','CD8A','ITGAE','FOXP3','IL17A','PDCD1','HAVCR2','LAG3','CTLA4','TCF7','CXCL13'),'GC','mc_annot','GC_2dproj',fat=sc_info,figwidth = 40, figheight=50, figtextsize = 8, dotsize = 0.2, stroke=0.05, spec = "Tcell")



TRM <- NormalizeData(TRM, normalization.method = "LogNormalize", scale.factor = 10000)
TRM <- ScaleData(TRM, features = rownames(TRM))
library(clusterProfiler)
library(enrichplot)
geneList<-

TRM <- irGSEA.score(object = TRM, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 10,
                             min.cells = 0, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens",  
                              geneid = "symbol",
                             method = c("AUCell"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')
ridgeplot1 <- irGSEA.ridgeplot(object = TRM,group.by='mc',
                              method = "AUCell",
                              show.geneset = c("HALLMARK-PI3K-AKT-MTOR-SIGNALING"))
ridgeplot2<-irGSEA.ridgeplot(object = TRM,group.by='mc',
                              method = "AUCell",
                              show.geneset = c("HALLMARK-MTORC1-SIGNALING"))


ridgeplot3<-irGSEA.ridgeplot(object = TRM,group.by='mc',
                              method = "AUCell",
                              show.geneset = c("HALLMARK-TNFA-SIGNALING-VIA-NFKB"))

gene_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','MX1','RSAD2','TOP2A','MKI67','STMN1')
tab_fn = "./lateral_gmods_trm.txt"
mcell_mat_rpt_cor_anchors(mat_id='TRM', gene_anchors = gene_anchors, cor_thresh = 0.1,
                          gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table('./lateral_gmods_trm.txt', header=T)
foc_genes = apply(gcor_mat[, intersect(colnames(gcor_mat),gene_anchors)], 1, which.max)

mcell_mat_ignore_genes(new_mat_id='TRM', mat_id='TRM', pre_bad_genes, reverse=F)
mat = scdb_mat("TRM")

mcell_add_gene_stat(gstat_id="TRM", mat_id="TRM", force=T)
mcell_gset_filter_varmean(gset_id="TRM_feats", gstat_id="TRM", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "TRM_feats", gstat_id="TRM", T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id="TRM", gset_id="TRM_feats")
gset <- scdb_gset("TRM_feats")
pst_genes <- names(gset@gene_set)
pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                 "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                 "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
pst_ex_genes <- c()
pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
pst_add_genes <- c()
final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
final_genes


gset = gset_new_gset(sets = gset@gene_set[final_genes], desc = "final genes")
scdb_add_gset("TRM_feats", gset)

mcell_add_cgraph_from_mat_bknn(mat_id="TRM",
                gset_id = "TRM_feats",
                graph_id="TRM_graph",
                K=85,
                dsamp=T)###一般K为总细胞数的平方根的数目，若细胞数很少，K=20-40

mcell_coclust_from_graph_resamp(
                coc_id="TRM_coc1000",
                graph_id="TRM_graph",
                min_mc_size=20,
                p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
                coc_id="TRM_coc1000",
                mat_id= "TRM",
                mc_id= "TRM_mc",
                K=30, min_mc_size=30, alpha=2) 

mc = scdb_mc("TRM_mc")
mc@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc@mc_fp))
scdb_add_mc("TRM_mc",mc)
mcell_mc2d_force_knn(mc2d_id="TRM_2dproj",mc_id="TRM_mc", graph_id="TRM_graph")
mcell_mc2d_plot(mc2d_id="TRM_2dproj")#



GC <- NormalizeData(GC, normalization.method = "LogNormalize", scale.factor = 10000)
GC <- FindVariableFeatures(GC, selection.method = "vst", nfeatures = 3000)
GC <- ScaleData(GC, features = rownames(GC))
GC <- RunPCA(GC, features = VariableFeatures(object = GC))
RunHarmony(GC,"orig.ident")->GC
GC <- FindNeighbors(GC, dims = 1:30,reduction = "harmony")
GC <- FindClusters(GC, resolution = 0.3)
GC <- RunUMAP(GC, dims = 1:30,reduction = "harmony")
DimPlot(GC,reduction='umap',label=T)


data.frame(GC@reductions$umap@cell.embeddings)->sc_info
sc_info$type=GC$seurat_clusters
ggplot()+
geom_point(data=sc_info,aes(x=UMAP_1,y=UMAP_2,fill=type),color='black',pch=21,size=2)+
scale_fill_chameleon()

FeaturePlot(GC,features='ITGAE',reduction='umap',order=T)+
    scale_colour_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000))
#####
countexp.Seurat<-sc.metabolism.Seurat(obj = TRM, method = "VISION", imputation = F, ncores = 12, metabolism.type = "KEGG")

datamean=group_by(score,group) %>% summarize_each(funs(mean))

base_mean = rowMeans(mean)
mat_scaled = t(scale(t(mean)))


col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
circos.heatmap(mat_scaled,  col = col_fun,dend.side = "inside",rownames.side = "outside",rownames.cex=4)
circos.clear()


subset(TRM,TLS=='low')->low
FeaturePlot(low,features='CXCL13',reduction='umap',order=T)+
    scale_colour_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000))

theme_publa <- function(base_size = 8, base_family = "sans", legend_position = "right",
                     title_size=8){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
  theme(legend.position = legend_position, legend.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000"), axis.ticks = element_line(colour = "#000000"), 
        legend.key = element_blank(), 
        axis.text = element_text(size = base_size, face="plain"), plot.title=element_text(face="plain", size = title_size),
        axis.title = element_text(face="plain", size = base_size), legend.text=element_text(size = base_size),
        legend.title=element_text(face="plain", size = base_size), strip.text=element_text(face="plain", size = base_size)
        )
}




ggplot()+
geom_point(data=data.frame(sc_info),aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
geom_point(data=data.frame(sc_info[which(sc_info$S=='pos'),]),aes(x=sc_x,y=sc_y,color=CXCL13,fill=CXCL13),size=0.7)+
scale_colour_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)) +
scale_fill_gradientn(colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                             c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                               "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)) +
       theme_publa(base_size = 8) +
      theme(legend.position = "bottom", legend.background = element_blank(), 
            legend.text = element_text(size=8), 
            legend.title = element_text(size=8), plot.title = element_text(size=8),
            axis.line.x = element_blank(), axis.line.y = element_blank(),
            axis.ticks = element_blank(), axis.text = element_blank(),
            legend.key.size = unit(1.5,"mm"))      


res.cut <- surv_cutpoint(OS, time = "OS.time", event = "OS",variables = c("score"))
  res.cat <- surv_categorize(res.cut)
  fit<-survfit(Surv(OS.time/30,OS)~score,data=res.cat)

  ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'),legend='none')


dat2[["percent.mt"]] <- PercentageFeatureSet(dat2, pattern ="^MT-")
plot1 <- FeatureScatter(H07T, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(H07T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
H07T <- subset(H07T, subset = nFeature_RNA > 200 & nFeature_RNA < 6000  & percent.mt<25 )


mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="CD3D")



SpatialFeaturePlot(object,features='TLS',images = NULL,image.alpha = 0)+
scale_fill_gradientn(limits = c(0,1),colours = colorRampPalette(#c("white", "orange", "tomato","mediumorchid4", "midnightblue"))
                               c("white","white","white", "#FEF8E5","#FFEABC","#FDBC52","#F68523",
                                 "#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)) 

pdf("./ig_genes_3L.pdf",useDingbats=F) 
SpatialFeaturePlot(data,features='IGHG1',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='IGHA1',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='IGHM',images = NULL,image.alpha = 0)
dev.off()
pdf("./COL1A1_3L.pdf",useDingbats=F) 
SpatialFeaturePlot(data,features='COL1A1',images = NULL,image.alpha = 0)
dev.off()

pdf("./HCC_3T_CELL.pdf",useDingbats=F)
SpatialFeaturePlot(data,features='Fibroblasts',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='T.cells',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='B.cells',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='NK.cells',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='Myeloid.cells',images = NULL,image.alpha = 0)
SpatialFeaturePlot(data,features='Endothelial.cells',images = NULL,image.alpha = 0)
dev.off()

SpatialFeaturePlot(data,features='LTB',images = NULL,image.alpha = 0)


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern ="^MT-")
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000  & percent.mt<25 )


data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
data <- ScaleData(data, features = rownames(data))
data <- RunPCA(data, features = VariableFeatures(object = data))
RunHarmony(data,"orig.ident")->data
data <- FindNeighbors(data, dims = 1:30,reduction = "harmony")
data <- FindClusters(data, resolution = 0.3)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony")
DimPlot(data,reduction='umap',label=T)
FeaturePlot(data,reduction='umap',features='CD79A')





minnumbcells <- min(sapply(c(as.character(unique(sc_info$pid))), function(x) length(rownames(sc_info[sc_info$pid==x,]))))
set.seed(22)
plotcells <-  as.vector(sapply(c(as.character(unique(sc_info$pid))), function(x) sample(rownames(sc_info[sc_info$pid==x,]), size=minnumbcells, replace=F, prob=NULL )))

pdf('./density_2d_group.pdf',useDingbats=F,width=10,height=8)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$pid=='yes'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$pid=='yes'),],aes(x=sc_x,y=sc_y),color="black")+
theme_classic()


ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y),size=0.7,color='#E1DFDD')+
stat_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$pid=='no'),],aes(x=sc_x,y=sc_y,fill=stat(density)),geom = "raster", contour = FALSE,alpha=0.5)+
scale_fill_gradientn(colours = colorRampPalette(c('white','#FFF8F7','#FFEBE4','#FFE2D6','#FFBCA7','#FF997F','#FF7A5B','#FF6344','#FF5737','#FF3314','#FF0011'))(1000))+
geom_density_2d(data=sc_info[plotcells,][which(sc_info[plotcells,]$pid=='no'),],aes(x=sc_x,y=sc_y),color="black")+
theme_classic()
dev.off()


ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=2)+
scale_fill_manual(values=c('CD4-CXCL13-TNFRSF18'='#44CEE0','CD4-Treg'='#6591CF','Th17-like'='#7FBE8C','CD4-FOS'='#FF0095','CD4-CXCL13-TCF7'='#FFD600','CD4-Naive'='#4E6D25','CD8-CTL'='#BA6FB5','CD8-TRM-IFNG'='#3551AA','CD8-Cycling'='#B41322','CD8-MT2A'='#F89373','CD8-CXCL13'='#215050','CD8-Naive'='#A1DBEA','CD8-ISG15'='#FF2D3C','CD8-IL17A'='#915D63','CD8-TRM-XCL2'='#96362B','CD4-MAIT'='#E26500'))
dev.off()

ggplot(aa,aes(x=Var2,y=per,fill=Var1))+
geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=c('CD4-CXCL13-TNFRSF18'='#44CEE0','CD4-Treg'='#6591CF','Th17-like'='#7FBE8C','CD4-FOS'='#FF0095','CD4-CXCL13-TCF7'='#FFD600','CD4-Naive'='#4E6D25','CD8-CTL'='#BA6FB5','CD8-TRM-IFNG'='#3551AA','CD8-Cycling'='#B41322','CD8-MT2A'='#F89373','CD8-CXCL13'='#215050','CD8-Naive'='#A1DBEA','CD8-ISG15'='#FF2D3C','CD8-IL17A'='#915D63','CD8-TRM-XCL2'='#96362B','CD4-MAIT'='#E26500'))+
theme_classic()

pdf("./sc2d.pdf",useDingbats=F,width=10,height=8)
ggplot()+
geom_point(data=sc_info,aes(x=sc_x,y=sc_y,fill=type),color='black',pch=21,size=2)+
ggsci::scale_fill_d3('category20')
dev.off()


ggplot(aa,aes(x=Var2,y=percentage,fill=Var1))+
geom_bar(stat = 'identity', position = 'fill')+
scale_fill_manual(values=color)+
theme_classic()

p<-list()
for(i in unique(bb$Var1)){
p[[i]]<-ggplot(bb[which(bb$Var1==i),],aes(x=group,y=per,color=group))+
geom_boxplot(lwd=1)+
geom_point(aes(color=group),size=1,position='jitter')+
scale_colour_manual(values=c('grey','black'))+
scale_fill_manual(values=c('grey','black'))+
theme(axis.text.x = element_text(angle =45, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title =element_blank(),legend.position = 'none')+
        stat_compare_means(method='t.test')+
ylab("percentage of CD3+ T cells")+ ggtitle(i)
}

mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="IL17A",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="CD4",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="CD8A",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="ITGAE",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="CXCL13",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="FOXP3",min_lfp=-1.5,max_lfp=1.5)
mcell_mc2d_plot_gene(mc2d_id="GC_2dproj",gene="GZMB",min_lfp=-1.5,max_lfp=1.5)

pdf("./survival.pdf")
for(i in 20:60){
  sig1<-sig[1:i]
  gsva(as.matrix(stad),gset.idx.list=list(sig1),method='ssgsea',ssgsea.norm=T)->score
  survival$score=score[1,rownames(survival)]
  res.cut <- surv_cutpoint(survival, time = "OS.time", event = "OS",
   variables = c("score"))
  res.cat <- surv_categorize(res.cut)
  fit<-survfit(Surv(OS.time/30,OS)~score,data=res.cat)
  p<-ggsurvplot(fit,pval=T,palette=c('#E7B800','#38A2DF'))
  print(p)

}
dev.off()

res_scDC_noClust1 <- scDC::scDC_noClustering(as.character(sc_info$type), as.character(sc_info$Sp), calCI = TRUE, 
                                     calCI_method = c("BCa"),
                                     nboot = 1000,ncores=10)


res1=res_scDC_noClust1$result
res1$mean=apply(res_scDC_noClust1$thetastar,1,mean)

res1[,c('cellTypes','subject','mean')]->res1
res1<-data.frame(type=as.character(res1$cellTypes[1:15]),high=res1$mean[1:15],low=res1$mean[16:30])
res1$fc=res1$low/res1$high
res1$logfc=log2(res1$fc)
res_scDC_noClust1$thetastar->sample_matrix
rownames(sample_matrix)=paste0(res_scDC_noClust1$result$cellTypes,'_',res_scDC_noClust1$result$subject)
for(i in unique(res1$type)){
  res1[i,'pvalue']=wilcox.test(sample_matrix[paste0(i,'_low'),],sample_matrix[paste0(i,'_high'),])$p.value
}

for(i in 1:nrow(res1)){
    if(res1$logfc[i]>0.5){
        res1$reg[i]='UP'
    }else if(res1$logfc[i]<(-0.5)){
         res1$reg[i]='DOWN'
    }else{
        res1$reg[i]='NAN'
    }
}




per2<-data.frame(compare=c(rep('lowvshigh',15)))
per2$celltype=rep(res1$type[1:15],1)
per2$logfc=res1$logfc
per2$reg=res1$reg
ggplot(per2,aes(x=compare,y=celltype))+
geom_point(aes(size=abs(logfc),color=reg))+
scale_color_manual(values=c('UP'='#ED2A28','DOWN'='#2C9685','NAN'='#B5B6B8'))+
theme_bw()+
theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
