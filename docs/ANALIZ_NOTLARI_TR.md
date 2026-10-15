> **Not:** Bu belge, orijinal Türkçe Word dokümanındaki metinlerin GitHub üzerinden okunabilen sürümüdür. Görseller ve özgün sayfa düzeni için [Word belgesini](Atopic_Dermatitis_Analysis_Notes_TR.docx) inceleyin.

# ATOPİK DERMATİT

Atopik dermatit, diğer adıyla atopik egzama, ciltte kuruluk, yoğun kaşıntı ve tekrarlayan inflamatuvar lezyonlarla seyreden kronik bir deri hastalığıdır. Genellikle çocukluk döneminde başlasa da ergenlikte veya yetişkinlikte de görülebilir. Hastalık bulaşıcı değildir ve alevlenme-remisyon dönemleriyle ilerler (Carr et al., 2024).

Hastalığın gelişiminde genetik yatkınlık, deri bariyerinin bozulması, bağışıklık sisteminin düzensiz aktivasyonu ve çevresel faktörler birlikte rol oynar. Sağlıklı deride epidermal bariyer su kaybını sınırlar ve dış etkenlerin deri içine girişini engeller. Atopik dermatitte bu bariyer zayıfladığı için deri kurur, irritanlara ve alerjenlere karşı daha geçirgen hâle gelir ve inflamasyon kolaylaşır (Afshari et al., 2024; Yue et al., 2024).

Atopik dermatitin en belirgin belirtisi kaşıntıdır. Kaşıma, deri bariyerini daha fazla bozarak kızarıklık, çatlama, kabuklanma ve inflamasyonu artırır. Böylece kaşıntı-kaşıma döngüsü oluşur. Lezyonların yerleşimi yaşa göre değişebilir; çocuklarda dirsek ve diz içleri, yetişkinlerde ise el, boyun, yüz ve fleksural bölgeler sık etkilenir (Carr et al., 2024).

Sabunlar, deterjanlar, parfümlü ürünler, terleme, sıcaklık değişimleri, stres ve bazı çevresel alerjenler hastalığı alevlendirebilir. Bununla birlikte atopik dermatit yalnızca alerjik bir hastalık değildir ve her hastada belirli bir alerjen saptanmaz (National Institute for Health and Care Excellence \[NICE\], 2025).

Tanı genellikle klinik muayene, kaşıntı, lezyonların görünümü ve hastalığın tekrarlayıcı seyri değerlendirilerek konur. Tedavide temel amaç deri bariyerini güçlendirmek, kaşıntı ve inflamasyonu azaltmak ve alevlenmeleri kontrol altına almaktır. Nemlendiriciler, topikal antiinflamatuvar tedaviler ve orta-şiddetli olgularda sistemik veya biyolojik tedaviler kullanılabilir (American Academy of Dermatology Association, 2025).

## Bu çalışma açısından önemi

Bu çalışmada atopik dermatitli bireylerden alınan lesional ve non-lesional deri örnekleri karşılaştırılmıştır. Lesional deri aktif lezyonlu bölgeleri, non-lesional deri ise aynı hastalarda görünür lezyon bulunmayan bölgeleri ifade eder. Non-lesional deri tamamen sağlıklı deriyle aynı kabul edilmemelidir; görünür lezyon olmasa da moleküler ve immün değişiklikler bulunabilir.

Tek hücreli RNA dizileme analizi sayesinde keratinosit, fibroblast, endotelyal ve bağışıklık hücreleri ayrı ayrı değerlendirilmiş; bu hücrelerdeki gen ekspresyon farklılıkları, biyolojik yolaklar ve hücre-hücre iletişimi incelenmiştir.

# KOD AÇIKLAMALARI

## 1. Başlangıç ayarları

############################################################

# GSE147424 - Atopic dermatitis scRNA-seq

setwd("~/Desktop/GSE147424_RAW")

set.seed(1234)

############################################################

setwd() analizde kullanılacak çalışma klasörünü belirler. Kod içinde göreli olarak belirtilen tüm çıktı yolları bu klasörün altında değerlendirilir. set.seed(1234), rastgele başlangıçlardan etkilenebilen analizlerin yeniden çalıştırıldığında mümkün olduğunca aynı sonucu vermesini sağlar. 1234 sayısının biyolojik bir anlamı yoktur; önemli olan seed değerinin sabit tutulmasıdır.

## 2. Paketlerin çağrılması

## 0. Libraries

############################################################

library(Seurat)

library(dplyr)

library(ggplot2)

library(patchwork)

library(Matrix)

library(SingleR)

library(celldex)

library(SingleCellExperiment)

library(scDblFinder)

library(clusterProfiler)

library(org.Hs.eg.db)

library(enrichplot)

library(AnnotationDbi)

library(CellChat)

library(NMF)

library(ggalluvial)

library(ComplexHeatmap)

library(harmony)

############################################################

library() daha önce kurulmuş R paketlerini aktif oturuma yükler. Seurat tek hücreli RNA-seq analizinin ana paketidir. dplyr veri düzenleme; ggplot2 grafik oluşturma; patchwork grafikleri birleştirme; SingleR ve celldex referans tabanlı hücre tipi anotasyonu; scDblFinder doublet tespiti; clusterProfiler ve ilgili anotasyon paketleri GO/KEGG analizi; CellChat hücre-hücre iletişimi; ComplexHeatmap ısı haritaları; harmony ise örnek kaynaklı batch etkisinin azaltılması için kullanılır.

## 3. Çıktı klasörlerinin oluşturulması

## 1. Output folders

############################################################

dir.create("plots", showWarnings = FALSE)

dir.create("Doublet_results", showWarnings = FALSE)

dir.create("DEG_GO_KEGG_by_celltype", showWarnings = FALSE)

dir.create("Global_DEG_GO_KEGG", showWarnings = FALSE)

dir.create("CellChat_results", showWarnings = FALSE)

dir.create("Metascape_inputs", showWarnings = FALSE)

dir.create("Metascape_inputs/By_celltype", showWarnings = FALSE)

dir.create("Metascape_inputs/Global", showWarnings = FALSE)

############################################################

dir.create() analiz çıktılarının düzenli biçimde kaydedileceği klasörleri oluşturur. showWarnings = FALSE, klasör zaten varsa uyarı verilmesini engeller. Grafikler, doublet tabloları, hücre tipi bazlı ve global DEG/GO/KEGG çıktıları, CellChat sonuçları ve Metascape listeleri ayrı klasörlerde tutulur.

## 4. Yardımcı fonksiyonların tanımlanması

## 2. Helper functions

############################################################

save_gg \<- function(plot, filename, width = 8, height = 6) {

print(plot)

ggsave(

filename = filename,

plot = plot,

width = width,

height = height

)

}

save_heatmap \<- function(ht, filename, width = 10, height = 8) {

ComplexHeatmap::draw(ht)

pdf(filename, width = width, height = height)

ComplexHeatmap::draw(ht)

dev.off()

}

safe_name \<- function(x) {

gsub("\[^A-Za-z0-9\_-\]", "\_", x)

}

run_enrichment \<- function(

deg_table,

prefix,

outdir,

showCategory = 15

) {



deg_sig \<- deg_table %\>%

filter(p_val_adj \< 0.05 & abs(avg_log2FC) \> 0.25)



message(prefix, " significant DEG count: ", nrow(deg_sig))



if (nrow(deg_sig) \< 5) {

message(prefix, ": fewer than 5 significant DEGs, skipping GO/KEGG")

return(NULL)

}



entrez \<- bitr(

deg_sig\$gene,

fromType = "SYMBOL",

toType = "ENTREZID",

OrgDb = org.Hs.eg.db

)



entrez_ids \<- unique(entrez\$ENTREZID)



message(prefix, " ENTREZ count: ", length(entrez_ids))



if (length(entrez_ids) \< 5) {

message(prefix, ": too few ENTREZ IDs, skipping GO/KEGG")

return(NULL)

}



##########################################################

Bu bölüm tekrar eden işlemleri tek bir fonksiyon çağrısıyla yürütmek için yardımcı fonksiyonlar tanımlar. save_gg() ggplot grafiklerini hem ekrana basar hem dosyaya kaydeder. save_heatmap() ComplexHeatmap nesnelerini PDF aygıtı üzerinden güvenli biçimde kaydeder. safe_name() hücre tipi adlarındaki boşluk ve özel karakterleri dosya adına uygun alt çizgilere dönüştürür. run_enrichment() ise anlamlı DEG seçimi, gen kimliği dönüşümü ve GO/KEGG analizini otomatikleştirir.

# 4.1 GO Biological Process analizi

# GO Biological Process

##########################################################



ego \<- enrichGO(

gene = entrez_ids,

OrgDb = org.Hs.eg.db,

keyType = "ENTREZID",

ont = "BP",

pAdjustMethod = "BH",

qvalueCutoff = 0.05

)



if (!is.null(ego) && nrow(as.data.frame(ego)) \> 0) {



write.csv(

as.data.frame(ego),

file = file.path(outdir, paste0(prefix, "\_GO_BP.csv")),

row.names = FALSE

)



p_go \<- dotplot(ego, showCategory = showCategory) +

ggtitle(paste("GO BP -", prefix))



save_gg(

p_go,

filename = file.path(outdir, paste0(prefix, "\_GO_BP.pdf")),

width = 9,

height = 7

)



} else {

message(prefix, ": no significant GO terms")

}



##########################################################

Önce düzeltilmiş p-değeri 0,05’in altında ve mutlak avg_log2FC değeri 0,25’ten büyük genler seçilir. Beşten az anlamlı gen varsa zenginleştirme güvenilir olmayabileceği için analiz atlanır. bitr() gen sembollerini ENTREZ ID’ye çevirir. enrichGO() yalnızca Biological Process kategorisini inceler; Benjamini-Hochberg yöntemi çoklu test düzeltmesi için kullanılır. Sonuç varsa ayrıntılı tablo CSV, dotplot ise PDF olarak kaydedilir.

# 4.2 KEGG zenginleştirme analizi

# KEGG

##########################################################



ekegg \<- enrichKEGG(

gene = entrez_ids,

organism = "hsa",

pvalueCutoff = 0.05

)



if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) \> 0) {



write.csv(

as.data.frame(ekegg),

file = file.path(outdir, paste0(prefix, "\_KEGG.csv")),

row.names = FALSE

)



p_kegg \<- dotplot(ekegg, showCategory = showCategory) +

ggtitle(paste("KEGG -", prefix))



save_gg(

p_kegg,

filename = file.path(outdir, paste0(prefix, "\_KEGG.pdf")),

width = 9,

height = 7

)



} else {

message(prefix, ": no significant KEGG terms")

}

}

############################################################

enrichKEGG(), anlamlı genlerin KEGG yolaklarında beklenenden fazla temsil edilip edilmediğini test eder. organism = "hsa" insan türünü belirtir. Sonuçlar tablo ve dotplot olarak kaydedilir. Hastalık veya enfeksiyon adı taşıyan KEGG yolakları, örnekte doğrudan o hastalığın ya da enfeksiyonun bulunduğunu kanıtlamaz; ortak mitokondriyal, inflamatuvar veya hücre iskeleti genlerinden kaynaklanabilir.

## 5. Seurat nesnesinin yüklenmesi

## 3. Load object

############################################################

combined \<- readRDS("~/Desktop/GSE147424_RAW/combined.rds")

dim(combined)

table(combined\$sample_id)

table(combined\$condition)

head(combined@meta.data)

############################################################

readRDS() daha önce oluşturulmuş combined.rds Seurat nesnesini R ortamına yükler. dim() nesnenin gen ve hücre boyutunu verir. table(sample_id) her örnekteki hücre sayısını; table(condition) lesional ve non-lesional hücre sayılarını gösterir. head(meta.data), analizde kullanılacak metadata sütunlarının doğru biçimde bulunduğunu kontrol etmeyi sağlar.

## 6. Kalite kontrolünün hazırlanması

## 4. QC

############################################################

# Bu veri setinde gerçek MT genleri görünmüyor.

grep("^MT-", rownames(combined), value = TRUE, ignore.case = TRUE)

# percent.mt güvenilir olmadığı için kaldırıyoruz.

if ("percent.mt" %in% colnames(combined@meta.data)) {

combined\$percent.mt \<- NULL

}

p_qc_sample \<- VlnPlot(

combined,

features = c("nFeature_RNA", "nCount_RNA"),

group.by = "sample_id",

ncol = 2,

pt.size = 0.2

)

save_gg(

p_qc_sample,

"plots/QC_before_filter_by_sample.pdf",

width = 12,

height = 6

)

p_qc_condition \<- VlnPlot(

combined,

features = c("nFeature_RNA", "nCount_RNA"),

group.by = "condition",

ncol = 2,

pt.size = 0.2

)

save_gg(

p_qc_condition,

"plots/QC_before_filter_by_condition.pdf",

width = 10,

height = 6

)

summary(combined\$nFeature_RNA)

summary(combined\$nCount_RNA)

cat("Before filtering:", ncol(combined), "cells\n")

############################################################

grep("^MT-") gen adları arasında insan mitokondriyal genlerini arar. Bu nesnede güvenilir MT genleri bulunmadığı için percent.mt yeniden hesaplanamaz ve mevcut sütun kaldırılır. Ardından nFeature_RNA ve nCount_RNA dağılımları önce örnek, sonra koşul bazında violin grafiklerle incelenir. summary() temel dağılım istatistiklerini; ncol() ise filtreleme öncesi toplam hücre sayısını verir.

## 7. Kalite kontrol filtrelemesi

## 5. QC filtering

############################################################

combined_filtered \<- subset(

combined,

subset = nFeature_RNA \> 200 &

nFeature_RNA \< 6000

)

cat("After nFeature filtering:", ncol(combined_filtered), "cells\n")

table(combined_filtered\$sample_id)

table(combined_filtered\$condition)

p_qc_after \<- VlnPlot(

combined_filtered,

features = c("nFeature_RNA", "nCount_RNA"),

group.by = "sample_id",

ncol = 2,

pt.size = 0.2

)

save_gg(

p_qc_after,

"plots/QC_after_filter_by_sample.pdf",

width = 12,

height = 6

)

p_qc_after_condition \<- VlnPlot(

combined_filtered,

features = c("nFeature_RNA", "nCount_RNA"),

group.by = "condition",

ncol = 2,

pt.size = 0.2

)

save_gg(

p_qc_after_condition,

"plots/QC_after_filter_by_condition.pdf",

width = 10,

height = 6

)

p_scatter \<- FeatureScatter(

combined_filtered,

feature1 = "nFeature_RNA",

feature2 = "nCount_RNA"

)

save_gg(

p_scatter,

"plots/QC_scatter_nFeature_vs_nCount.pdf",

width = 7,

height = 6

)

############################################################

subset() ile nFeature_RNA değeri 200’den büyük ve 6000’den küçük hücreler tutulur. Alt eşik çok az gen içeren düşük kaliteli hücreleri, üst eşik ise aşırı yüksek gen sayısına sahip olası multipletleri sınırlar. Filtreleme sonrasında hücre sayıları örnek ve koşul bazında yeniden kontrol edilir. Violin ve FeatureScatter grafikleri filtrelemenin dağılımlara etkisini ve nFeature_RNA ile nCount_RNA arasındaki ilişkiyi gösterir.

## 8. scDblFinder ile doublet tespiti

## 6. Doublet detection with scDblFinder

############################################################

set.seed(1234)

cat(

"Doublet analizi öncesi toplam hücre:",

ncol(combined_filtered),

"\n"

)

# Not:

# Bu veri clean.data / processed matrix gibi görünüyor.

# Yine de pipeline gereği scDblFinder çalıştırıyoruz.

# Sonuçları dikkatli yorumlamak gerekir.

sce_doublet \<- as.SingleCellExperiment(

combined_filtered,

assay = "RNA"

)

sce_doublet \<- scDblFinder(

sce_doublet,

samples = "sample_id",

verbose = TRUE

)

combined_filtered\$doublet_score \<- sce_doublet\$scDblFinder.score

combined_filtered\$doublet_class \<- sce_doublet\$scDblFinder.class

table(combined_filtered\$doublet_class)

table(

combined_filtered\$sample_id,

combined_filtered\$doublet_class

)

table(

combined_filtered\$condition,

combined_filtered\$doublet_class

)

############################################################

Seurat nesnesi scDblFinder’ın çalışabildiği SingleCellExperiment biçimine dönüştürülür. samples = "sample_id" parametresi doublet modelinin her örneği ayrı değerlendirmesini sağlar. scDblFinder her hücre için bir skor ve singlet/doublet sınıfı üretir. Bu bilgiler Seurat metadata’sına eklenir ve toplam, örnek bazlı ve koşul bazlı tablolarla kontrol edilir. İşlenmiş bir matris kullanıldığı için doublet sonuçlarının dikkatli yorumlanması gerektiği kodda ayrıca belirtilmiştir.

# 8.1 Doublet özet tabloları

# Doublet summary tables

############################################################

doublet_summary_sample \<- as.data.frame(

table(

sample_id = combined_filtered\$sample_id,

doublet_class = combined_filtered\$doublet_class

)

)

write.csv(

doublet_summary_sample,

file = "Doublet_results/Doublet_summary_by_sample.csv",

row.names = FALSE

)

doublet_summary_condition \<- as.data.frame(

table(

condition = combined_filtered\$condition,

doublet_class = combined_filtered\$doublet_class

)

)

write.csv(

doublet_summary_condition,

file = "Doublet_results/Doublet_summary_by_condition.csv",

row.names = FALSE

)

############################################################

table() ile örnek ve doublet sınıfının çapraz sayımı oluşturulur, as.data.frame() sonucu standart tabloya çevirir ve write.csv() ile kaydeder. Aynı işlem koşul bazında da yapılır. Böylece hangi örnek veya grupta daha fazla doublet tahmin edildiği sayısal olarak izlenebilir.

# 8.2 Doublet skoru violin grafiği

# Doublet score violin plot

############################################################

p_doublet_score \<- VlnPlot(

combined_filtered,

features = "doublet_score",

group.by = "sample_id",

pt.size = 0.1

)

save_gg(

p_doublet_score,

"plots/Doublet_score_by_sample.pdf",

width = 12,

height = 6

)

############################################################

VlnPlot(), her örnekteki doublet_score dağılımını gösterir. Düşük skorlar singlet olasılığının daha yüksek, yüksek skorlar doublet olasılığının daha yüksek olduğunu düşündürür. pt.size parametresi tek tek hücre noktalarının boyutunu belirler. Grafik hem ekranda gösterilir hem PDF olarak kaydedilir.

# 8.3 Doublet ve singlet oran grafiği

# Doublet / singlet fraction barplot

############################################################

doublet_plot_data \<- combined_filtered@meta.data %\>%

dplyr::count(

sample_id,

doublet_class,

name = "n"

) %\>%

dplyr::group_by(sample_id) %\>%

dplyr::mutate(

fraction = n / sum(n)

) %\>%

dplyr::ungroup()

write.csv(

doublet_plot_data,

file = "Doublet_results/Doublet_fraction_by_sample.csv",

row.names = FALSE

)

p_doublet_fraction \<- ggplot(

doublet_plot_data,

aes(

x = sample_id,

y = fraction,

fill = doublet_class

)

) +

geom_col() +

theme_classic() +

labs(

x = "Sample",

y = "Cell fraction",

fill = "Doublet class",

title = "Doublet / singlet fraction by sample"

) +

theme(

axis.text.x = element_text(

angle = 45,

hjust = 1

)

)

save_gg(

p_doublet_fraction,

"plots/Doublet_fraction_by_sample.pdf",

width = 10,

height = 6

)

############################################################

dplyr::count() her örnekteki singlet ve doublet sayılarını hesaplar. group_by(sample_id) ve mutate(fraction = n / sum(n)) sınıfların örnek içindeki oranını oluşturur. ggplot ve geom_col() bu oranları yığılmış sütun grafiği olarak gösterir. X ekseni etiketleri okunabilirlik için 45 derece döndürülür.

# 8.4 Doubletlerin çıkarılması

# Keep singlets only

############################################################

combined_filtered \<- subset(

combined_filtered,

subset = doublet_class == "singlet"

)

cat(

"Doublet çıkarıldıktan sonra toplam hücre:",

ncol(combined_filtered),

"\n"

)

table(combined_filtered\$sample_id)

table(combined_filtered\$condition)

saveRDS(

combined_filtered,

file = "GSE147424_AD_QC_singlets.rds"

)

############################################################

subset(doublet_class == "singlet") yalnızca singlet olarak sınıflandırılan hücreleri tutar. Doubletlerin çıkarılması, iki farklı hücre tipinin markerlarını aynı anda taşıyan yapay hücrelerin kümeleme ve anotasyonu bozmasını azaltır. Kalan hücre sayıları yeniden kontrol edilir ve kalite kontrolü tamamlanmış nesne RDS olarak kaydedilir.

## 9. Normalizasyon, değişken gen seçimi ve ölçekleme

## 7. Normalize / Variable features / Scale

############################################################

combined_filtered \<- NormalizeData(

combined_filtered,

normalization.method = "LogNormalize",

scale.factor = 10000

)

combined_filtered \<- FindVariableFeatures(

combined_filtered,

selection.method = "vst",

nfeatures = 2500

)

top20 \<- head(VariableFeatures(combined_filtered), 20)

p_var \<- VariableFeaturePlot(combined_filtered)

p_var_labeled \<- LabelPoints(

plot = p_var,

points = top20,

repel = TRUE,

xnudge = 0,

ynudge = 0

)

save_gg(

p_var,

"plots/VariableFeaturePlot.pdf",

width = 7,

height = 6

)

save_gg(

p_var_labeled,

"plots/VariableFeaturePlot_top20_labeled.pdf",

width = 7,

height = 6

)

save_gg(

p_var + p_var_labeled,

"plots/VariableFeaturePlot_side_by_side.pdf",

width = 14,

height = 6

)

combined_filtered \<- ScaleData(

combined_filtered,

features = VariableFeatures(combined_filtered),

vars.to.regress = "nCount_RNA"

)

############################################################

NormalizeData(LogNormalize), her hücredeki gen sayımlarını toplam kütüphane büyüklüğüne göre normalize eder, 10.000 ile ölçekler ve log dönüşümü uygular. FindVariableFeatures(vst, nfeatures = 2500), hücreler arasındaki biyolojik farklılığı en iyi taşıyan 2500 geni seçer. VariableFeaturePlot ve LabelPoints seçilen genleri görselleştirir. ScaleData genleri standardize eder ve nCount_RNA etkisini regresyonla azaltır.

## 10. PCA, Harmony, UMAP, t-SNE ve kümeleme

## 8. PCA / Harmony / UMAP / tSNE / Clustering

############################################################

combined_filtered \<- RunPCA(

combined_filtered,

features = VariableFeatures(combined_filtered),

verbose = TRUE

)

p_elbow \<- ElbowPlot(

combined_filtered,

ndims = 40

)

save_gg(

p_elbow,

"plots/ElbowPlot_40PCs.pdf",

width = 7,

height = 6

)

############################################################

RunPCA yüksek boyutlu gen ekspresyon verisini temel bileşenlere indirger. ElbowPlot ilk 40 bileşenin katkısını göstererek kullanılacak PC sayısının seçimine yardım eder. Sonraki analizlerde ilk 20 boyut dims_use olarak tanımlanmıştır. Bu ana bölüm Harmony düzeltmesi, komşuluk grafiği, kümeleme, UMAP ve t-SNE hesaplamalarının temelini oluşturur.

# 10.1 PCA ısı haritaları

# PCA heatmaps

############################################################

p_heat_1_10 \<- DimHeatmap(

combined_filtered,

dims = 1:10,

cells = 500,

balanced = TRUE,

ncol = 5

)

save_gg(

p_heat_1_10,

"plots/PCA_DimHeatmap_PC1_to_PC10.pdf",

width = 15,

height = 11

)

############################################################

DimHeatmap ilk 10 temel bileşene en fazla katkı sağlayan genlerin 500 hücredeki ekspresyon örüntülerini gösterir. balanced = TRUE, pozitif ve negatif yük taşıyan genlerin dengeli seçilmesini sağlar. ncol = 5, on PC panelini beşer sütunlu biçimde düzenler.

# 10.2 PCA grafiklerinin oluşturulması

# PCA plots

############################################################

p_pca_sample \<- DimPlot(

combined_filtered,

reduction = "pca",

group.by = "sample_id",

pt.size = 0.3

)

p_pca_condition \<- DimPlot(

combined_filtered,

reduction = "pca",

group.by = "condition",

pt.size = 0.3

)

save_gg(

p_pca_sample + p_pca_condition,

"plots/PCA_sample_condition_side_by_side.pdf",

width = 16,

height = 6

)

############################################################

DimPlot(reduction = "pca") hücreleri PCA uzayında gösterir. İlk grafik örnek kimliğine, ikinci grafik koşula göre renklendirilir. Bu iki grafik birlikte incelenerek örnek kaynaklı ayrışma ile lesional/non-lesional dağılımının PCA’daki görünümü değerlendirilir.

# 10.3 Harmony batch düzeltmesi ve kümeleme

# Harmony + clustering

############################################################

combined_filtered \<- RunHarmony(

object = combined_filtered,

group.by.vars = "sample_id"

)

dims_use \<- 1:20

combined_filtered \<- FindNeighbors(

combined_filtered,

reduction = "harmony",

dims = dims_use

)

combined_filtered \<- FindClusters(

combined_filtered,

resolution = 0.5

)

combined_filtered \<- RunUMAP(

combined_filtered,

reduction = "harmony",

dims = dims_use

)

combined_filtered \<- RunTSNE(

combined_filtered,

reduction = "harmony",

dims = dims_use,

check_duplicates = FALSE

)

############################################################

RunHarmony, sample_id değişkenine göre örnek kaynaklı batch etkisini PCA uzayında azaltır. FindNeighbors Harmony boyutlarında benzer hücrelerin komşuluk ağını oluşturur. FindClusters resolution = 0.5 ile bu ağdaki hücreleri kümelere ayırır. RunUMAP ve RunTSNE batch düzeltilmiş ilk 20 Harmony boyutu üzerinden iki boyutlu gösterimler üretir.

# 10.4 UMAP ve t-SNE görselleştirmeleri

# UMAP / tSNE

############################################################

p_umap_condition \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "condition",

pt.size = 0.3

)

p_umap_sample \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "sample_id",

pt.size = 0.3

)

p_umap_cluster \<- DimPlot(

combined_filtered,

reduction = "umap",

label = TRUE,

pt.size = 0.3

)

save_gg(

p_umap_condition + p_umap_sample + p_umap_cluster,

"plots/UMAP_condition_sample_cluster_combined.pdf",

width = 21,

height = 6

)

p_tsne_condition \<- DimPlot(

combined_filtered,

reduction = "tsne",

group.by = "condition",

pt.size = 0.3

)

p_tsne_cluster \<- DimPlot(

combined_filtered,

reduction = "tsne",

label = TRUE,

pt.size = 0.3

)

save_gg(

p_tsne_condition + p_tsne_cluster,

"plots/tSNE_condition_cluster_side_by_side.pdf",

width = 16,

height = 6

)

table(combined_filtered\$seurat_clusters, combined_filtered\$sample_id)

table(combined_filtered\$seurat_clusters, combined_filtered\$condition)

saveRDS(

combined_filtered,

file = "GSE147424_AD_filtered_clustered.rds"

)

############################################################

UMAP grafikleri koşul, örnek ve Seurat küme numarasına göre ayrı ayrı hazırlanıp yan yana birleştirilir. t-SNE grafikleri koşul ve küme bazında oluşturulur. Farklı örneklerin aynı kümelerde karışması batch düzeltmesinin başarılı olduğuna dair görsel destek sağlar; kümelerin iki yöntemde benzer görünmesi ise küme yapısının tutarlılığını destekler. Çapraz tablolar her kümenin örnek ve koşul dağılımını gösterir.

## 11. Küme markerlarının belirlenmesi

## 9. Cluster markers

############################################################

cluster_markers \<- FindAllMarkers(

combined_filtered,

only.pos = TRUE,

logfc.threshold = 0.25,

min.pct = 0.25

)

write.csv(

cluster_markers,

file = "GSE147424_cluster_markers.csv",

row.names = FALSE

)

top_markers \<- cluster_markers %\>%

group_by(cluster) %\>%

slice_max(order_by = avg_log2FC, n = 10)

write.csv(

top_markers,

file = "GSE147424_top10_markers_per_cluster.csv",

row.names = FALSE

)

############################################################

FindAllMarkers her kümeyi diğer tüm kümelerle karşılaştırarak pozitif marker genleri belirler. only.pos = TRUE yalnızca ilgili kümede daha yüksek genleri tutar. logfc.threshold ve min.pct, minimum etki büyüklüğü ile ifade edilen hücre oranı eşikleridir. Her küme için avg_log2FC değeri en yüksek ilk 10 marker ayrıca seçilip kaydedilir.

## 12. Manuel deri markerlarının kontrolü

## 10. Manual skin marker check

############################################################

skin_markers \<- c(

"KRT14", "KRT5", "KRT1", "KRT10", "KRT16", "KRT6A",

"COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",

"CD3D", "CD3E", "TRAC", "IL7R", "CD4", "CD8A",

"NKG7", "GNLY", "GZMB", "PRF1",

"MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN",

"LYZ", "CD68", "C1QA", "C1QB", "LST1",

"FCER1A", "CD1A", "CD207", "CLEC10A", "LAMP3",

"PECAM1", "VWF", "CLDN5",

"RGS5", "ACTA2", "TAGLN"

)

skin_markers_present \<- skin_markers\[

skin_markers %in% rownames(combined_filtered)

\]

p_dot \<- DotPlot(

combined_filtered,

features = skin_markers_present

) +

RotatedAxis()

save_gg(

p_dot,

"plots/Skin_marker_DotPlot.pdf",

width = 14,

height = 7

)

p_feature \<- FeaturePlot(

combined_filtered,

features = skin_markers_present\[1:min(12, length(skin_markers_present))\],

ncol = 3

)

print(p_feature)

ggsave(

filename = "plots/Skin_marker_FeaturePlot_first12.png",

plot = p_feature,

width = 12,

height = 10,

dpi = 300

)

############################################################

Deri dokusunda beklenen keratinosit, fibroblast, lenfosit, myeloid, dendritik, endotelyal ve perivasküler markerlardan bir panel oluşturulur. Yalnızca veri setinde bulunan genler seçilir. DotPlot nokta büyüklüğüyle ifade eden hücre oranını, renk yoğunluğuyla ortalama ekspresyonu gösterir. FeaturePlot ilk 12 markerın UMAP üzerindeki hücre bazlı dağılımını gösterir.

## 13. HPCA referansıyla SingleR anotasyonu

## 11. SingleR annotation - HPCA

############################################################

refHPCA \<- HumanPrimaryCellAtlasData()

sce_test \<- GetAssayData(

object = combined_filtered,

assay = "RNA",

layer = "data"

)

prediction_HPCA_main \<- SingleR(

test = sce_test,

ref = refHPCA,

clusters = Idents(combined_filtered),

labels = refHPCA\$label.main

)

prediction_HPCA_fine \<- SingleR(

test = sce_test,

ref = refHPCA,

clusters = Idents(combined_filtered),

labels = refHPCA\$label.fine

)

combined_filtered\$HPCA_main \<- prediction_HPCA_main\$labels\[

Idents(combined_filtered)

\]

combined_filtered\$HPCA_fine \<- prediction_HPCA_fine\$labels\[

Idents(combined_filtered)

\]

p_hpca_main \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "HPCA_main",

label = TRUE,

pt.size = 0.4

) +

NoLegend()

p_hpca_fine \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "HPCA_fine",

label = TRUE,

pt.size = 0.4

) +

NoLegend()

save_gg(

p_hpca_main + p_hpca_fine,

"plots/UMAP_HPCA_main_fine_side_by_side.pdf",

width = 18,

height = 7

)

table(combined_filtered\$HPCA_main)

table(combined_filtered\$HPCA_fine)

############################################################

HumanPrimaryCellAtlasData HPCA referans ekspresyon atlasını sağlar. GetAssayData(layer = "data") normalize edilmiş RNA matrisini alır. SingleR, Seurat kümelerini referans profillerle karşılaştırarak ana ve ayrıntılı hücre tipi etiketleri üretir. Etiketler küme kimlikleri üzerinden her hücrenin metadata’sına yazılır ve UMAP üzerinde yan yana gösterilir. Otomatik anotasyonlar manuel markerlarla doğrulanmalıdır.

## 14. Monaco referansıyla SingleR anotasyonu

## 12. SingleR annotation - Monaco immune

############################################################

refMonaco \<- MonacoImmuneData()

prediction_Monaco_main \<- SingleR(

test = sce_test,

ref = refMonaco,

clusters = Idents(combined_filtered),

labels = refMonaco\$label.main

)

prediction_Monaco_fine \<- SingleR(

test = sce_test,

ref = refMonaco,

clusters = Idents(combined_filtered),

labels = refMonaco\$label.fine

)

combined_filtered\$Monaco_main \<- prediction_Monaco_main\$labels\[

Idents(combined_filtered)

\]

combined_filtered\$Monaco_fine \<- prediction_Monaco_fine\$labels\[

Idents(combined_filtered)

\]

p_monaco_main \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "Monaco_main",

label = TRUE,

pt.size = 0.4

) +

NoLegend()

p_monaco_fine \<- DimPlot(

combined_filtered,

reduction = "umap",

group.by = "Monaco_fine",

label = TRUE,

pt.size = 0.4

) +

NoLegend()

save_gg(

p_monaco_main + p_monaco_fine,

"plots/UMAP_Monaco_main_fine_side_by_side.pdf",

width = 18,

height = 7

)

saveRDS(

combined_filtered,

file = "GSE147424_AD_with_celltypes.rds"

)

############################################################

MonacoImmuneData özellikle bağışıklık hücrelerinin daha ayrıntılı sınıflandırılmasına yönelik bir referanstır. HPCA ile aynı SingleR yaklaşımı kullanılarak ana ve ayrıntılı etiketler hesaplanır, metadata’ya eklenir ve UMAP üzerinde gösterilir. Monaco sonuçları özellikle immün popülasyonları doğrulayan tamamlayıcı bir anotasyon olarak değerlendirilmelidir.

## 15. Hücre tipine özgü DEG ve GO/KEGG analizi

## 13. Cell-type-specific DEG + GO/KEGG

############################################################

seurat_obj \<- combined_filtered

celltype_col \<- "HPCA_main"

condition_col \<- "condition"

meta \<- seurat_obj@meta.data

celltypes \<- na.omit(unique(meta\[\[celltype_col\]\]))

for (ct in celltypes) {



message("Processing cell type: ", ct)



cells_ct \<- rownames(meta)\[

meta\[\[celltype_col\]\] == ct &

!is.na(meta\[\[celltype_col\]\])

\]



obj_ct \<- subset(seurat_obj, cells = cells_ct)



tab \<- table(obj_ct@meta.data\[\[condition_col\]\])

print(tab)



if (!all(c("Lesional", "Non_lesional") %in% names(tab))) {

message(ct, ": missing Lesional or Non_lesional")

next

}



if (any(tab\[c("Lesional", "Non_lesional")\] \< 10)) {

message(ct, ": too few cells")

next

}



Idents(obj_ct) \<- obj_ct@meta.data\[\[condition_col\]\]



deg \<- FindMarkers(

obj_ct,

ident.1 = "Lesional",

ident.2 = "Non_lesional",

logfc.threshold = 0.25,

min.pct = 0.1,

test.use = "wilcox"

)



deg\$gene \<- rownames(deg)



prefix \<- paste0(

"Lesional_vs_Nonlesional\_",

safe_name(ct)

)



write.csv(

deg,

file = file.path(

"DEG_GO_KEGG_by_celltype",

paste0(prefix, "\_DEG.csv")

),

row.names = FALSE

)



############################################################

Bu döngü her HPCA_main hücre tipini ayrı seçer ve kendi içinde Lesional ile Non_lesional gruplarını karşılaştırır. Her iki koşulun bulunması ve her grupta en az 10 hücre olması zorunludur. FindMarkers Wilcoxon testiyle DEG’leri hesaplar. Pozitif avg_log2FC lesional, negatif değer non-lesional grupta yüksek ekspresyonu gösterir. Her hücre tipinin DEG tablosu kaydedilir ve run_enrichment fonksiyonu ile GO/KEGG analizi yürütülür.

# 15.1 Hücre tipi bazlı Metascape girdileri

# Metascape input - cell-type specific

############################################################



deg_sig \<- deg %\>%

filter(p_val_adj \< 0.05 & abs(avg_log2FC) \> 0.25)



if (nrow(deg_sig) \> 0) {



genes_all \<- unique(deg_sig\$gene)



genes_lesional_up \<- deg_sig %\>%

filter(avg_log2FC \> 0) %\>%

pull(gene) %\>%

unique()



genes_nonlesional_up \<- deg_sig %\>%

filter(avg_log2FC \< 0) %\>%

pull(gene) %\>%

unique()



writeLines(

genes_all,

con = file.path(

"Metascape_inputs/By_celltype",

paste0("Metascape\_", prefix, "\_all_DEGs.txt")

)

)



writeLines(

genes_lesional_up,

con = file.path(

"Metascape_inputs/By_celltype",

paste0("Metascape\_", prefix, "\_Lesional_up.txt")

)

)



writeLines(

genes_nonlesional_up,

con = file.path(

"Metascape_inputs/By_celltype",

paste0("Metascape\_", prefix, "\_Nonlesional_up.txt")

)

)

}



run_enrichment(

deg_table = deg,

prefix = prefix,

outdir = "DEG_GO_KEGG_by_celltype",

showCategory = 15

)

}

############################################################

Anlamlı hücre tipi DEG’leri tümü, lesional grupta yüksek olanlar ve non-lesional grupta yüksek olanlar şeklinde üç listeye ayrılır. writeLines() her geni ayrı satıra yazarak Metascape’e yüklenebilecek TXT dosyaları oluşturur. Yönlere göre ayrı listeler kullanmak, zenginleşen süreçlerin hangi koşulla ilişkili olduğunu daha doğru yorumlamayı sağlar.

## 16. Global DEG ve GO/KEGG analizi

## 14. Global DEG + GO/KEGG

############################################################

Idents(seurat_obj) \<- "condition"

deg_global \<- FindMarkers(

seurat_obj,

ident.1 = "Lesional",

ident.2 = "Non_lesional",

logfc.threshold = 0.25,

min.pct = 0.1,

test.use = "wilcox"

)

deg_global\$gene \<- rownames(deg_global)

write.csv(

deg_global,

file = "Global_DEG_GO_KEGG/Global_DEG_Lesional_vs_Nonlesional.csv",

row.names = FALSE

)

############################################################

Global analiz bütün hücreleri condition kimliğine göre Lesional ve Non_lesional olarak karşılaştırır. Hücre tipi ayrımı yapılmadığı için genel hastalık sinyalini verir ancak hücre tipi kompozisyonundan etkilenebilir. DEG tablosu kaydedilir ve showCategory = 20 ile global GO/KEGG analizi yürütülür. Sonuçlar hücre tipi bazlı analizlerle birlikte değerlendirilmelidir.

# 16.1 Global Metascape girdileri

# Metascape input - global DEG

############################################################

deg_global_sig \<- deg_global %\>%

filter(p_val_adj \< 0.05 & abs(avg_log2FC) \> 0.25)

if (nrow(deg_global_sig) \> 0) {



global_genes_all \<- unique(deg_global_sig\$gene)



global_lesional_up \<- deg_global_sig %\>%

filter(avg_log2FC \> 0) %\>%

pull(gene) %\>%

unique()



global_nonlesional_up \<- deg_global_sig %\>%

filter(avg_log2FC \< 0) %\>%

pull(gene) %\>%

unique()



writeLines(

global_genes_all,

con = "Metascape_inputs/Global/Metascape_Global_all_DEGs.txt"

)



writeLines(

global_lesional_up,

con = "Metascape_inputs/Global/Metascape_Global_Lesional_up.txt"

)



writeLines(

global_nonlesional_up,

con = "Metascape_inputs/Global/Metascape_Global_Nonlesional_up.txt"

)

}

run_enrichment(

deg_table = deg_global,

prefix = "Global_Lesional_vs_Nonlesional",

outdir = "Global_DEG_GO_KEGG",

showCategory = 20

)

\############################################################

Global anlamlı DEG’ler tümü, lesional-up ve non-lesional-up olarak ayrılır ve Metascape_inputs/Global klasörüne yazılır. Böylece lesional dokuda artan süreçler ile non-lesional dokuda daha yüksek olan homeostatik veya yapısal süreçler ayrı analiz edilebilir.

## 17. CellChat analizine hazırlık

## 15. CellChat

############################################################

seurat_lesional \<- subset(seurat_obj, subset = condition == "Lesional")

seurat_nonlesional \<- subset(seurat_obj, subset = condition == "Non_lesional")

cat("Lesional hücre sayısı:", ncol(seurat_lesional), "\n")

table(seurat_lesional\$HPCA_main)

cat("Non-lesional hücre sayısı:", ncol(seurat_nonlesional), "\n")

table(seurat_nonlesional\$HPCA_main)

############################################################

Seurat nesnesi condition değişkenine göre lesional ve non-lesional olmak üzere iki ayrı nesneye bölünür. Her koşuldaki toplam hücre sayısı ve HPCA_main hücre tipi dağılımı kontrol edilir. CellChat hücre sayısı ve grup kompozisyonundan etkilenebildiği için bu başlangıç kontrolleri yorumlama açısından gereklidir.

# 17.1 Lesional CellChat analizi

# CellChat - Lesional

############################################################

set.seed(1234)

cellchat_lesional \<- createCellChat(

object = seurat_lesional,

group.by = "HPCA_main",

assay = "RNA"

)

CellChatDB \<- CellChatDB.human

cellchat_lesional@DB \<- CellChatDB

cellchat_lesional \<- subsetData(cellchat_lesional)

cellchat_lesional \<- identifyOverExpressedGenes(cellchat_lesional)

cellchat_lesional \<- identifyOverExpressedInteractions(cellchat_lesional)

cellchat_lesional \<- computeCommunProb(cellchat_lesional, type = "triMean")

cellchat_lesional \<- filterCommunication(cellchat_lesional, min.cells = 10)

df_lesional \<- subsetCommunication(cellchat_lesional)

write.csv(

df_lesional,

file = "CellChat_results/CellChat_Lesional_communications.csv",

row.names = FALSE

)

cellchat_lesional \<- computeCommunProbPathway(cellchat_lesional)

cellchat_lesional \<- aggregateNet(cellchat_lesional)

groupSize_lesional \<- as.numeric(table(cellchat_lesional@idents))

netVisual_circle(

cellchat_lesional@net\$count,

vertex.weight = groupSize_lesional,

weight.scale = TRUE,

label.edge = FALSE,

title.name = "Lesional - Number of interactions"

)

pdf("CellChat_results/CellChat_Lesional_circle_count.pdf", width = 8, height = 8)

netVisual_circle(

cellchat_lesional@net\$count,

vertex.weight = groupSize_lesional,

weight.scale = TRUE,

label.edge = FALSE,

title.name = "Lesional - Number of interactions"

)

dev.off()

cellchat_lesional \<- netAnalysis_computeCentrality(

cellchat_lesional,

slot.name = "netP"

)

p1 \<- netAnalysis_signalingRole_heatmap(

cellchat_lesional,

pattern = "outgoing"

)

p2 \<- netAnalysis_signalingRole_heatmap(

cellchat_lesional,

pattern = "incoming"

)

p1 + p2

pdf(

"CellChat_results/CellChat_Lesional_outgoing_incoming.pdf",

width = 18,

height = 8

)

ComplexHeatmap::draw(p1 + p2)

dev.off()

saveRDS(

cellchat_lesional,

file = "CellChat_results/cellchat_lesional.rds"

)

############################################################

createCellChat lesional nesneyi HPCA_main hücre tiplerine göre iletişim ağına dönüştürür. İnsan ligand-reseptör veritabanı atanır. subsetData sinyal genlerini seçer; identifyOverExpressedGenes ve identifyOverExpressedInteractions ilgili genleri ve ligand-reseptör çiftlerini belirler. computeCommunProb(triMean) iletişim olasılıklarını hesaplar; min.cells = 10 küçük grupları filtreler. İletişimler CSV’ye çıkarılır, yolak düzeyinde toplanır ve dairesel ağ grafiği oluşturulur. Merkeziyet hesabından sonra outgoing ve incoming sinyal rolleri yan yana ısı haritası olarak kaydedilir.

# 17.2 Non-lesional CellChat analizi

# CellChat - Non-lesional

############################################################

set.seed(1234)

cellchat_nonlesional \<- createCellChat(

object = seurat_nonlesional,

group.by = "HPCA_main",

assay = "RNA"

)

cellchat_nonlesional@DB \<- CellChatDB.human

cellchat_nonlesional \<- subsetData(cellchat_nonlesional)

cellchat_nonlesional \<- identifyOverExpressedGenes(cellchat_nonlesional)

cellchat_nonlesional \<- identifyOverExpressedInteractions(cellchat_nonlesional)

cellchat_nonlesional \<- computeCommunProb(cellchat_nonlesional, type = "triMean")

cellchat_nonlesional \<- filterCommunication(cellchat_nonlesional, min.cells = 10)

df_nonlesional \<- subsetCommunication(cellchat_nonlesional)

write.csv(

df_nonlesional,

file = "CellChat_results/CellChat_Nonlesional_communications.csv",

row.names = FALSE

)

cellchat_nonlesional \<- computeCommunProbPathway(cellchat_nonlesional)

cellchat_nonlesional \<- aggregateNet(cellchat_nonlesional)

groupSize_nonlesional \<- as.numeric(table(cellchat_nonlesional@idents))

netVisual_circle(

cellchat_nonlesional@net\$count,

vertex.weight = groupSize_nonlesional,

weight.scale = TRUE,

label.edge = FALSE,

title.name = "Non-lesional - Number of interactions"

)

pdf("CellChat_results/CellChat_Nonlesional_circle_count.pdf", width = 8, height = 8)

netVisual_circle(

cellchat_nonlesional@net\$count,

vertex.weight = groupSize_nonlesional,

weight.scale = TRUE,

label.edge = FALSE,

title.name = "Non-lesional - Number of interactions"

)

dev.off()

cellchat_nonlesional \<- netAnalysis_computeCentrality(

cellchat_nonlesional,

slot.name = "netP"

)

p3 \<- netAnalysis_signalingRole_heatmap(

cellchat_nonlesional,

pattern = "outgoing"

)

p4 \<- netAnalysis_signalingRole_heatmap(

cellchat_nonlesional,

pattern = "incoming"

)

p3 + p4

pdf(

"CellChat_results/CellChat_Nonlesional_outgoing_incoming.pdf",

width = 18,

height = 8

)

ComplexHeatmap::draw(p3 + p4)

dev.off()

saveRDS(

cellchat_nonlesional,

file = "CellChat_results/cellchat_nonlesional.rds"

)

############################################################

Non-lesional veri için lesional bölümdeki CellChat adımlarının aynısı uygulanır. Aynı veritabanı, hücre tipi etiketi, triMean yöntemi ve minimum hücre eşiği kullanılması iki koşulun karşılaştırılabilirliğini sağlar. İletişim tablosu, dairesel ağ, merkeziyet ölçümleri ve outgoing/incoming heatmap ayrı dosyalara kaydedilir.

## 18. CellChat koşul karşılaştırması

## 16. CellChat comparison

############################################################

object.list \<- list(

Non_lesional = cellchat_nonlesional,

Lesional = cellchat_lesional

)

cellchat_merged \<- mergeCellChat(

object.list,

add.names = names(object.list)

)

compareInteractions(

cellchat_merged,

show.legend = FALSE,

group = c(1, 2),

measure = "count"

)

compareInteractions(

cellchat_merged,

show.legend = FALSE,

group = c(1, 2),

measure = "weight"

)

p_count \<- compareInteractions(

cellchat_merged,

show.legend = FALSE,

group = c(1, 2),

measure = "count"

)

p_weight \<- compareInteractions(

cellchat_merged,

show.legend = FALSE,

group = c(1, 2),

measure = "weight"

)

ggsave(

filename = "CellChat_results/CellChat_compare_count_weight.pdf",

plot = p_count + p_weight,

width = 12,

height = 6

)

netVisual_diffInteraction(

cellchat_merged,

weight.scale = TRUE,

measure = "count"

)

pdf("CellChat_results/CellChat_difference_count.pdf", width = 8, height = 8)

netVisual_diffInteraction(

cellchat_merged,

weight.scale = TRUE,

measure = "count"

)

dev.off()

netVisual_diffInteraction(

cellchat_merged,

weight.scale = TRUE,

measure = "weight"

)

pdf("CellChat_results/CellChat_difference_weight.pdf", width = 8, height = 8)

netVisual_diffInteraction(

cellchat_merged,

weight.scale = TRUE,

measure = "weight"

)

dev.off()

saveRDS(

cellchat_merged,

file = "CellChat_results/cellchat_merged_Lesional_Nonlesional.rds"

)

############################################################

İki CellChat nesnesi listede Non_lesional ve Lesional adlarıyla birleştirilir. compareInteractions(measure = "count") toplam etkileşim sayısını, measure = "weight" toplam iletişim gücünü karşılaştırır. netVisual_diffInteraction hücre tipi çiftleri arasındaki fark ağlarını sayı ve ağırlık açısından gösterir. Yön yorumlanırken object.list sırasının birinci grup Non_lesional, ikinci grup Lesional olduğu dikkate alınmalıdır.

## 19. Son özet tabloları ve nesnenin kaydedilmesi

## 17. Summary tables and final object

############################################################

cat("\n===== CONDITION =====\n")

table(seurat_obj\$condition)

cat("\n===== SAMPLE =====\n")

table(seurat_obj\$sample_id)

cat("\n===== CLUSTERS x CONDITION =====\n")

table(seurat_obj\$seurat_clusters, seurat_obj\$condition)

cat("\n===== HPCA MAIN x CONDITION =====\n")

table(seurat_obj\$HPCA_main, seurat_obj\$condition)

cat("\n===== HPCA FINE =====\n")

table(seurat_obj\$HPCA_fine)

saveRDS(

seurat_obj,

file = "GSE147424_AD_final_analysis_object.rds"

)

############################################################

Analiz sonunda koşul, örnek, küme × koşul, HPCA_main × koşul ve HPCA_fine dağılımları konsola yazdırılır. Bu tablolar son hücre kompozisyonunu kontrol etmeyi sağlar. Tüm normalizasyon, Harmony, kümeleme ve anotasyon bilgilerini içeren son Seurat nesnesi RDS olarak kaydedilir.

# 19.1 Metascape girdi özetinin hazırlanması

# Metascape input summary

############################################################

metascape_files \<- list.files(

"Metascape_inputs",

pattern = "\\.txt\$",

recursive = TRUE,

full.names = TRUE

)

metascape_summary \<- data.frame(

file = metascape_files,

gene_count = sapply(

metascape_files,

function(x) length(readLines(x))

)

)

write.csv(

metascape_summary,

file = "Metascape_inputs/Metascape_input_summary.csv",

row.names = FALSE

)

print(metascape_summary)

message("Analysis completed successfully.")

list.files() Metascape_inputs klasörü ve alt klasörlerindeki bütün TXT dosyalarını bulur. sapply() her dosyadaki satır sayısını hesaplar; her satır bir gen olduğu için bu değer gen sayısını verir. Dosya yolları ve gen sayıları Metascape_input_summary.csv dosyasına kaydedilir. Son mesaj kodun son satıra hata vermeden ulaştığını bildirir; tek başına biyolojik sonuçların doğruluğunu garanti etmez.

# ANALİZLER

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 1- QC before filter by sample

Şekil 1’de kalite kontrol filtreleri uygulanmadan önce her bir örnekte bulunan hücrelerin **nFeature_RNA** ve **nCount_RNA** değerleri violin plot ile gösterilmiştir. nFeature_RNA, her hücrede tespit edilen farklı gen sayısını; nCount_RNA ise hücrede ölçülen toplam transkript veya UMI sayısını ifade etmektedir. Violin grafiklerinin genişliği, ilgili değer aralığında kaç hücre bulunduğunu; grafik üzerine eklenen siyah noktalar ise tek tek hücreleri temsil etmektedir.

Örnekler arasında hem tespit edilen gen sayısı hem de toplam RNA miktarı bakımından belirli farklılıklar bulunduğu görülmektedir. Hücrelerin büyük bölümü orta düzeyde nFeature_RNA ve nCount_RNA değerlerine sahipken, bazı hücreler oldukça düşük veya oldukça yüksek değerler göstermektedir. Çok düşük nFeature_RNA ve nCount_RNA değerleri, yeterli RNA içermeyen, hasar görmüş veya düşük kaliteli hücrelerle ilişkili olabilir. Buna karşılık aşırı yüksek değerler, aynı damlacık içerisinde iki veya daha fazla hücrenin birlikte yakalanmış olabileceğini, yani olası doublet hücrelerini düşündürebilir.

Örnekler arasındaki dağılımların tamamen aynı olmaması; hücre tipi kompozisyonu, örnek kalitesi, hücre yakalama verimliliği ve sekanslama derinliği gibi hem biyolojik hem de teknik faktörlerden kaynaklanabilir. Bu nedenle Şekil 1, yalnızca veri kalitesini görselleştirmek için değil, aynı zamanda hücrelerin veri setine dâhil edilmesinde kullanılacak alt ve üst filtreleme eşiklerinin belirlenmesi açısından da önemlidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 2- GQ before filter by condition

Şekil 2’de filtreleme öncesindeki hücreler, lesional ve non-lesional koşullarına göre birleştirilerek karşılaştırılmıştır. Sol panelde hücre başına tespit edilen gen sayısı, sağ panelde ise hücre başına toplam RNA/UMI miktarı gösterilmektedir.

Her iki koşulda da geniş bir nFeature_RNA ve nCount_RNA dağılımı bulunmaktadır. Lesional grupta dağılımın üst bölgelere doğru kısmen daha geniş olduğu görülmektedir. Bu durum, lesional örneklerde bazı hücrelerin daha yüksek RNA içeriğine veya daha fazla tespit edilebilir gene sahip olabileceğini düşündürmektedir. Ancak bu farklılık tek başına biyolojik aktivite artışı olarak yorumlanmamalıdır. Hücre tipi kompozisyonundaki farklılıklar, örnek başına sekanslama derinliği veya doublet oranları da benzer bir görünüme neden olabilir.

İki koşulun dağılımlarının büyük ölçüde birbiriyle örtüşmesi, lesional ve non-lesional gruplar arasında kalite kontrol metrikleri açısından tamamen ayrışmış bir yapı bulunmadığını göstermektedir. Bununla birlikte, her iki grupta da düşük kaliteli ve aşırı yüksek RNA içeriğine sahip hücrelerin bulunduğu görülmüş ve bu hücrelerin sonraki analizlerden önce filtrelenmesi gerekli olmuştur.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 3 QC after filter by sample

Şekil 3’te kalite kontrol eşikleri uygulandıktan sonra veri setinde bırakılan hücrelerin örnek bazındaki nFeature_RNA ve nCount_RNA dağılımları gösterilmektedir. Filtreleme öncesindeki grafikle karşılaştırıldığında, aşırı düşük gen ve RNA sayısına sahip hücrelerin önemli ölçüde azaldığı görülmektedir. Benzer şekilde, belirlenen üst sınırların dışında kalan aşırı yüksek değerli hücreler de veri setinden uzaklaştırılmıştır.

Filtreleme sonrasında dağılımların daha kontrollü ve karşılaştırılabilir bir aralığa geldiği görülmektedir. Bununla birlikte örneklerin tamamen aynı dağılıma sahip olmaması beklenen bir durumdur. Çünkü farklı örneklerdeki hücre tipi oranları ve hücrelerin transkripsiyonel aktiviteleri birbirinden farklı olabilir. Örneğin metabolik olarak aktif veya büyük hücreler, doğal olarak daha yüksek RNA ve gen sayısına sahip olabilir.

Bu grafik, kalite kontrol işleminin veri setindeki uç değerleri azalttığını ve sonraki normalizasyon, PCA, UMAP, kümeleme ve diferansiyel gen ekspresyon analizleri için daha güvenilir bir hücre popülasyonu oluşturduğunu göstermektedir. Ancak filtreleme sonucunun yalnızca violin grafiklerine değil, hücre sayılarındaki değişime ve doublet analizine de birlikte bakılarak değerlendirilmesi gerekir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 4 QC after filter by condition

Şekil 4’te kalite kontrol sonrasında veri setinde kalan hücrelerin nFeature_RNA ve nCount_RNA değerleri lesional ve non-lesional koşulları arasında karşılaştırılmıştır. Filtreleme sonrasında her iki grupta da çok düşük ve aşırı yüksek değerlere sahip hücrelerin azaldığı ve dağılımların daha belirgin bir aralıkta toplandığı görülmektedir.

Lesional gruptaki nFeature_RNA ve nCount_RNA değerlerinin merkezi dağılımı non-lesional gruba göre kısmen daha yüksek görünmektedir. Bu farklılık, lezyonlu dokudaki hücrelerin daha yüksek transkripsiyonel aktivite göstermesiyle ilişkili olabilir. Bununla birlikte, lesional ve non-lesional grupların farklı hücre tipi kompozisyonlarına sahip olması da bu sonucu etkileyebilir. Örneğin yüksek RNA içeriğine sahip hücre tiplerinin lesional grupta daha fazla bulunması, koşullar arasında genel bir farklılık varmış gibi görünmesine yol açabilir.

Dolayısıyla bu grafik, gruplar arasında kalite metrikleri bakımından bazı farklılıklar olduğunu göstermekle birlikte, bu farklılıkların doğrudan hastalık biyolojisine bağlanması uygun değildir. Koşul etkisinin hücre tipi bazında ve örnek düzeyinde ayrıca değerlendirilmesi gerekmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 5 QC scatter nFeature vs nCount

Şekil 5’te hücre başına tespit edilen gen sayısı ile toplam RNA/UMI sayısı arasındaki ilişki saçılım grafiği ile gösterilmiştir. Her nokta bir hücreyi, renkler ise hücrelerin ait olduğu örnekleri temsil etmektedir. Grafikte hesaplanan korelasyon katsayısı **0,94** olup, nFeature_RNA ile nCount_RNA arasında oldukça güçlü ve pozitif bir ilişki bulunduğunu göstermektedir.

Toplam RNA sayısı arttıkça hücrede tespit edilen farklı gen sayısının da artması beklenen bir sonuçtur. Daha fazla RNA molekülü ölçülen hücrelerde daha fazla genin tespit edilebilmesi, veri setinin teknik olarak tutarlı olduğunu desteklemektedir. Noktaların büyük bölümünün ortak bir eğri üzerinde bulunması, örneklerin benzer bir gen sayısı–RNA miktarı ilişkisi sergilediğini göstermektedir.

Bununla birlikte, eğrinin genel yapısından belirgin biçimde uzaklaşan hücreler düşük kaliteli hücre, teknik uç değer veya doublet adayı olabilir. Özellikle benzer nFeature_RNA değerine rağmen olağan dışı yüksek nCount_RNA gösteren hücreler, az sayıda genin çok yoğun ekspresyonu veya olası doublet oluşumu açısından değerlendirilmelidir. Genel olarak yüksek korelasyon, RNA sayımı ile tespit edilen gen sayısının birbiriyle uyumlu olduğunu göstermektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 6 doublet score by sample

Şekil 6’da scDblFinder tarafından her hücre için hesaplanan **doublet_score** değerlerinin örneklere göre dağılımı gösterilmiştir. Doublet skoru, bir hücre barkodunun tek bir hücre yerine iki veya daha fazla hücreden kaynaklanma olasılığına ilişkin sayısal bir değerlendirmedir. Skor yükseldikçe hücrenin doublet olma olasılığı artmaktadır.

Örneklerin çoğunda hücrelerin büyük bölümü sıfıra yakın doublet skorlarına sahiptir. Bu durum, veri setindeki hücrelerin çoğunluğunun singlet özellik gösterdiğini desteklemektedir. Bununla birlikte bazı örneklerde, özellikle grafikte violin dağılımı daha geniş görünen örnekte, orta ve yüksek doublet skorlarına sahip hücrelerin diğer örneklere göre daha fazla olduğu görülmektedir. Bu örnekteki farklılık, hücre yükleme yoğunluğu, hücre konsantrasyonu veya örnek hazırlanması sırasında iki hücrenin aynı damlacık içerisinde yakalanma olasılığının daha yüksek olmasıyla ilişkili olabilir.

Grafikte bazı hücrelerin 1’e yakın doublet skorlarına sahip olduğu görülse de bu hücreler toplam popülasyonun küçük bir kısmını oluşturmaktadır. Doublet skorları tek başına kesin sınıflandırma değildir; scDblFinder tarafından oluşturulan singlet/doublet sınıflandırmasıyla birlikte değerlendirilmiştir. Yüksek skorlu ve doublet olarak sınıflandırılan hücrelerin çıkarılması, gerçekte bulunmayan ara hücre tiplerinin veya yapay kümelerin oluşmasını önlemek açısından önemlidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 7 doublet fraction by sample

Şekil 7’de her bir örnekteki singlet ve doublet hücrelerin oransal dağılımı yığılmış sütun grafiği ile gösterilmiştir. Y ekseni hücrelerin toplam içerisindeki oranını, x ekseni ise örnekleri temsil etmektedir. Kırmızı alanlar singlet, turkuaz alanlar doublet olarak sınıflandırılan hücreleri göstermektedir.

Tüm örneklerde hücrelerin çok büyük çoğunluğunun singlet olarak sınıflandırıldığı görülmektedir. Doublet oranı örnekler arasında küçük değişiklikler göstermekle birlikte genel olarak düşük düzeydedir. Bazı örneklerde doublet alanının diğerlerine kıyasla daha görünür olması, Şekil 6’da gözlenen doublet skoru dağılımlarıyla uyumludur.

Doublet oranlarının düşük olması, veri setinin genel kalitesi açısından olumlu bir bulgudur. Ancak düşük oranlarda bile doublet hücreler iki farklı hücre tipine ait belirteçleri aynı anda taşıyabilir. Bu durum hücre tipi anotasyonunu, marker gen analizini ve hücre-hücre iletişimi sonuçlarını etkileyebilir. Bu nedenle doublet olarak sınıflandırılan hücrelerin analizden çıkarılması, daha güvenilir kümeleme ve anotasyon sonuçları elde edilmesini sağlamıştır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 8 Variable Feature

Şekil 8’de genlerin ortalama ekspresyon değerleri ile standartlaştırılmış varyansları arasındaki ilişki gösterilmektedir. Her nokta bir geni temsil etmektedir. Siyah noktalar normal varyasyon gösteren genleri, kırmızı noktalar ise analiz için seçilen yüksek değişken genleri göstermektedir. Analizde toplam **2500 yüksek değişken gen** belirlenmiştir.

Tek hücreli RNA dizileme analizlerinde tüm genler PCA ve kümeleme için eşit derecede bilgi taşımaz. Birçok gen tüm hücrelerde benzer düzeyde eksprese edilirken, bazı genler hücre tipleri veya hücresel durumlar arasında daha fazla değişkenlik gösterir. Bu yüksek değişken genler, hücreler arasındaki biyolojik farklılıkları yakalamada daha fazla bilgi sağlar.

Sağ panelde en yüksek değişkenlik gösteren genlerden bazıları etiketlenmiştir. Bu genler arasında epitel farklılaşması, hücre iskeleti, bağışıklık yanıtı veya dokuya özgü işlevlerle ilişkili genlerin bulunması, veri setinde farklı hücresel popülasyonların bulunduğunu desteklemektedir. Ancak yüksek değişkenlik her zaman yalnızca biyolojik farklılıktan kaynaklanmaz. Hücre döngüsü, stres yanıtı, mitokondriyal aktivite veya teknik varyasyon da bazı genlerin değişkenliğini artırabilir. Bu nedenle yüksek değişken gen seçimi, sonraki PCA ve kümeleme analizleriyle birlikte değerlendirilmiştir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 9 Elbow plot

Şekil 9’da ilk 40 temel bileşenin açıkladığı varyasyon, her bir temel bileşene ait standart sapma değerleri üzerinden gösterilmiştir. X ekseni temel bileşen numarasını, y ekseni ise ilgili bileşenin standart sapmasını temsil etmektedir.

İlk temel bileşenlerde standart sapmanın oldukça yüksek olduğu ve özellikle ilk birkaç bileşende hızlı bir düşüş meydana geldiği görülmektedir. Bu durum, veri setindeki varyasyonun önemli bir bölümünün ilk temel bileşenler tarafından açıklandığını göstermektedir. Daha sonraki bileşenlerde eğri giderek yataylaşmakta ve her yeni bileşenin açıkladığı ek bilgi azalmaktadır.

Grafikte kesin ve tek bir kırılma noktası bulunmamakla birlikte, yaklaşık **15–20. temel bileşenden sonra** eğimin belirgin biçimde azaldığı görülmektedir. Bu nedenle UMAP, t-SNE ve komşuluk grafiği oluşturma aşamalarında ilk yaklaşık 20 temel bileşenin kullanılması, biyolojik sinyali korurken düşük bilgi içeren veya gürültü ağırlıklı bileşenlerin etkisini sınırlandırmak açısından uygun bir tercih olarak değerlendirilmiştir.

Elbow plot tek başına kesin bir istatistiksel test değildir. Bu nedenle seçilen temel bileşen sayısının PCA ısı haritaları, küme stabilitesi ve biyolojik marker dağılımlarıyla birlikte değerlendirilmesi gerekir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 10 PCA Heatmap

Şekil 10’da PCA analizinde elde edilen ilk 10 temel bileşene en fazla katkı sağlayan genlerin ekspresyon desenleri ısı haritaları şeklinde gösterilmiştir. Her panel farklı bir temel bileşeni temsil etmektedir. Sütunlar hücreleri, satırlar ise ilgili temel bileşene yüksek katkı sağlayan genleri göstermektedir.

Renk farklılıkları genlerin hücreler arasındaki göreceli ekspresyon düzeylerini temsil etmektedir. Sarı ve mor bölgeler, ilgili gen gruplarının farklı hücre gruplarında zıt ekspresyon örüntüleri sergilediğini göstermektedir. Bir temel bileşende belirgin gen bloklarının oluşması, bu bileşenin veri setindeki anlamlı bir biyolojik varyasyonu yakaladığını düşündürmektedir.

Örneğin bazı temel bileşenler epitel belirteçlerinin, bazıları stromal veya bağışıklık hücresi belirteçlerinin ekspresyon farklarıyla ilişkili olabilir. İlk 10 temel bileşende belirgin ve düzenli gen ekspresyon bloklarının bulunması, PCA’nın rastgele gürültü yerine hücre popülasyonları arasındaki gerçek transkripsiyonel farklılıkları yakaladığını desteklemektedir. Bununla birlikte her temel bileşenin biyolojik anlamı, o bileşene en fazla pozitif ve negatif katkı sağlayan genlerin ayrı olarak incelenmesiyle doğrulanmalıdır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 11 PCA sample and condition

Şekil 11’de hücrelerin ilk iki temel bileşen, yani PC1 ve PC2 üzerindeki dağılımı gösterilmektedir. Sol panelde hücreler örnek kimliğine, sağ panelde ise lesional ve non-lesional koşullarına göre renklendirilmiştir.

Örnek bazındaki grafikte farklı örneklerden gelen hücrelerin önemli ölçüde aynı alanlarda bir arada bulunduğu görülmektedir. Her örneğin tamamen ayrı bir bölge oluşturmadığı göz önüne alındığında, veri setinde tüm hücreleri örnek kimliğine göre ayıran çok güçlü ve tek yönlü bir batch etkisi bulunmadığı söylenebilir. Bununla birlikte bazı PCA kollarında belirli örneklerin daha yoğun bulunması, örnekler arasında kısmi kompozisyon veya teknik farklılıklar olabileceğini göstermektedir.

Koşula göre renklendirilen panelde lesional ve non-lesional hücrelerin büyük ölçüde örtüştüğü görülmektedir. Bu durum, hücrelerin yalnızca hastalık koşuluna göre iki ayrı gruba ayrılmadığını göstermektedir. Tek hücreli verilerde ana varyasyon kaynağının çoğu zaman koşuldan önce hücre tipi olması beklenir. Dolayısıyla lesional ve non-lesional hücrelerin aynı PCA alanlarında bulunması, iki koşulun benzer temel hücre tiplerini içerdiğini; koşula bağlı değişimlerin ise hücre tipleri içindeki daha ince transkripsiyonel farklılıklar şeklinde ortaya çıkabileceğini düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 12 UMAP condition sample and cluster

Şekil 12’de hücrelerin yüksek boyutlu gen ekspresyon profilleri UMAP yöntemiyle iki boyuta indirgenmiştir. Soldaki panelde hücreler koşula, ortadaki panelde örnek kimliğine, sağdaki panelde ise Seurat küme numaralarına göre renklendirilmiştir.

Koşula göre renklendirilen grafikte lesional ve non-lesional hücrelerin birçok UMAP bölgesinde birlikte yer aldığı görülmektedir. Bu durum, kümelerin yalnızca koşul farklılığı nedeniyle oluşmadığını ve iki grupta ortak hücre popülasyonlarının bulunduğunu göstermektedir. Bazı kümelerde koşullardan birinin daha yoğun bulunması ise ilgili hücre tipinin lesional veya non-lesional dokuda sayısal olarak daha fazla olabileceğini ya da koşula özgü transkripsiyonel durumlar geliştirmiş olabileceğini düşündürmektedir.

Örnek bazındaki UMAP grafiğinde farklı örneklerden gelen hücrelerin aynı kümelerde karışması olumlu bir sonuçtur. Bu görünüm, belirlenen kümelerin yalnızca tek bir örneğe ait hücrelerden oluşmadığını ve örnekler arası teknik farklılıkların küme yapısını tamamen belirlemediğini desteklemektedir. Ancak bazı kümelerde belirli örneklerin baskın olması, o hücre popülasyonunun yalnızca az sayıdaki bireyde yüksek bulunabileceğini göstermektedir. Bu durum diferansiyel ekspresyon analizlerinde örnek düzeyinin dikkate alınmasını önemli hâle getirir.

Küme numaralarına göre renklendirilen panelde yaklaşık 15 farklı transkripsiyonel hücre kümesi tanımlanmıştır. Birbirinden ayrılan kümeler farklı hücre tiplerini veya aynı hücre tipinin farklı aktivasyon ve farklılaşma durumlarını temsil edebilir. Küme numaraları tek başına hücre tipi anlamına gelmediğinden, sonraki aşamada marker gen ekspresyonları ve referans tabanlı anotasyon sonuçları kullanılarak biyolojik kimlikler atanmıştır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 13 t_SNE condition and cluster

Şekil 13’te hücrelerin t-SNE yöntemiyle oluşturulan iki boyutlu dağılımları gösterilmiştir. Sol panelde hücreler koşula, sağ panelde ise Seurat kümelerine göre renklendirilmiştir.

Koşul grafiğinde lesional ve non-lesional hücrelerin birçok bölgede birlikte bulunması, UMAP analizinde gözlenen sonuçlarla uyumludur. Bu durum, hücrelerin temel transkripsiyonel kimliğinin koşuldan çok hücre tipi tarafından belirlendiğini desteklemektedir. Bununla birlikte bazı t-SNE adacıklarında lesional veya non-lesional hücrelerin daha yoğun bulunması, belirli hücre popülasyonlarının koşula bağlı dağılım farklılıkları gösterebileceğine işaret etmektedir.

Küme bazındaki grafikte farklı hücre popülasyonlarının belirgin adacıklar hâlinde ayrıldığı görülmektedir. UMAP ve t-SNE yöntemleri farklı algoritmik yaklaşımlar kullandığından, iki yöntemde benzer kümelerin gözlenmesi kümeleme sonucunun tutarlılığını desteklemektedir. Ancak iki boyutlu haritalardaki kümeler arasındaki fiziksel uzaklığın doğrudan biyolojik uzaklık olarak yorumlanmaması gerekir. Özellikle t-SNE grafiğinde kümeler arası boşlukların büyüklüğü nicel bir benzerlik ölçüsü değildir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 14 skin marker dot plots

Şekil 14’te seçilmiş hücre tipi belirteçlerinin Seurat kümelerindeki ekspresyonları dot plot ile gösterilmiştir. X ekseni marker genleri, y ekseni ise küme numaralarını temsil etmektedir. Nokta büyüklüğü, ilgili kümedeki hücrelerin yüzde kaçının geni eksprese ettiğini; nokta renginin yoğunluğu ise o kümedeki ortalama ekspresyon düzeyini göstermektedir.

Büyük ve koyu renkli bir nokta, ilgili genin hem çok sayıda hücrede eksprese edildiğini hem de ortalama ekspresyon düzeyinin yüksek olduğunu gösterir. Küçük fakat koyu bir nokta, genin az sayıdaki hücrede güçlü şekilde eksprese edildiğini; büyük fakat açık renkli bir nokta ise çok sayıda hücrede düşük veya orta düzeyde ekspresyon bulunduğunu gösterebilir.

Keratin genlerinin belirli kümelerde birlikte yüksek ekspresyon göstermesi, bu kümelerin epitel veya keratinosit kökenli olduğunu desteklemektedir. Kollajenler, DCN ve LUM gibi ekstraselüler matriks ilişkili genlerin farklı kümelerde yoğunlaşması fibroblast veya stromal hücre popülasyonlarına işaret etmektedir. CD3D gibi lenfosit belirteçlerinin sınırlı kümelerde görülmesi ise T hücrelerinin tanımlanmasına katkı sağlamaktadır.

Marker genlerin tek tek değil, bir panel hâlinde birlikte değerlendirilmesi önemlidir. Çünkü tek bir gen birden fazla hücre tipinde düşük düzeyde eksprese edilebilir. Güvenilir anotasyon için bir hücre tipine özgü birden fazla pozitif markerın birlikte bulunması ve diğer hücre tiplerine ait markerların düşük olması beklenir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 15 skin marker feature plot

Şekil 15’te KRT14, KRT5, KRT1, KRT10, KRT16, KRT6A, COL1A1, COL1A2, DCN, LUM, COL3A1 ve CD3D genlerinin UMAP üzerindeki ekspresyon dağılımları gösterilmiştir. Gri renk düşük veya tespit edilemeyen ekspresyonu, mavi rengin yoğunlaşması ise daha yüksek ekspresyon seviyelerini temsil etmektedir.

KRT14 ve KRT5, bazal keratinositlerle ilişkili markerlardır ve belirli epitel kümelerinde yoğunlaşmaları bazal keratinosit popülasyonlarının varlığını desteklemektedir. KRT1 ve KRT10 daha farklılaşmış suprabazal keratinositlerle ilişkilidir. Bu genlerin farklı UMAP bölgelerinde yoğunlaşması, keratinositlerin farklılaşma durumlarına göre alt popülasyonlara ayrıldığını düşündürmektedir.

KRT16 ve KRT6A genleri proliferatif, stres altındaki veya aktive keratinositlerde artabilen markerlardır. Bu genlerin belirli kümelerde yoğunlaşması, lezyonla ilişkili aktive keratinosit durumlarının bulunabileceğine işaret edebilir. Bununla birlikte bu yorumun koşul bazındaki ekspresyon karşılaştırmasıyla doğrulanması gerekir.

COL1A1, COL1A2 ve COL3A1 kollajen üretimiyle; DCN ve LUM ise ekstraselüler matriks organizasyonuyla ilişkili stromal markerlardır. Bu genlerin ortak UMAP bölgelerinde eksprese edilmesi fibroblast popülasyonunun tanımlanmasını desteklemektedir. CD3D ekspresyonunun küçük ve ayrı bir bölgede yoğunlaşması ise T lenfositlerin daha sınırlı bir hücre popülasyonu olarak bulunduğunu göstermektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 16 HPCA main and fine

Şekil 16’da hücre tipleri HPCA referans veri seti kullanılarak ana ve ayrıntılı sınıflar şeklinde otomatik olarak anotasyonlanmıştır. Sol panel daha geniş hücre sınıflarını, sağ panel ise daha ayrıntılı hücre tipi atamalarını göstermektedir.

Ana anotasyon, hücreleri epitel, fibroblast, bağışıklık veya diğer geniş biyolojik kategorilere ayırmaktadır. Ayrıntılı anotasyon ise aynı geniş sınıf içerisindeki daha spesifik alt tipleri ayırt etmeye çalışmaktadır. Aynı UMAP bölgelerindeki hücrelerin benzer referans etiketleri alması, bu bölgelerin ortak bir transkripsiyonel kimliğe sahip olduğunu desteklemektedir.

Bununla birlikte referans tabanlı anotasyon sonuçları doğrudan kesin hücre kimliği olarak kabul edilmemelidir. HPCA genel amaçlı bir insan hücre atlasına dayandığı için deri dokusuna özgü bazı hücre durumlarını yeterince ayrıntılı ayıramayabilir. Ayrıca hastalık koşullarında değişen gen ekspresyonları, hücrelerin referanstaki en yakın fakat biyolojik olarak tam karşılığı olmayan bir etikete atanmasına yol açabilir.

Bu nedenle HPCA sonuçları, Şekil 14 ve Şekil 15’te gösterilen marker genlerle birlikte değerlendirilmiştir. Özellikle deri dokusunda beklenmeyen “kondrosit”, “nöron” veya benzeri anotasyonlar elde edildiğinde bu hücrelerin marker profilleri manuel olarak kontrol edilmelidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 17 monaco main and fine

Şekil 17’de hücrelerin Monaco referans veri setine göre yapılan ana ve ayrıntılı anotasyonları gösterilmiştir. Monaco referansı özellikle insan bağışıklık hücrelerinin alt tiplerinin belirlenmesi amacıyla geliştirilmiş bir referans veri setidir.

Ana anotasyon grafiğinde hücreler daha geniş bağışıklık kategorilerine ayrılırken, ayrıntılı anotasyon hücrelerin daha spesifik immün alt sınıflara atanmasını sağlamaktadır. UMAP üzerindeki belirli hücre bölgelerinin farklı Monaco etiketleri alması, bağışıklık hücrelerinin kendi içerisinde transkripsiyonel alt gruplara ayrıldığını düşündürmektedir.

Ancak Monaco referansının temel olarak bağışıklık hücrelerine odaklanması nedeniyle epitel ve stromal hücrelerin anotasyonunda HPCA veya dokuya özgü markerlar kadar uygun olmayabilir. Bu nedenle Monaco sonuçları özellikle T hücreleri, dendritik hücreler ve diğer immün popülasyonların doğrulanmasında tamamlayıcı olarak kullanılmıştır. HPCA, Monaco ve manuel marker incelemesinin birlikte değerlendirilmesi, tek bir referans veri setine dayanmanın oluşturabileceği yanlış sınıflandırma riskini azaltmaktadır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 18 Go Bp Lesional Vs Nonlesional Keratinocytes

Şekil 18’de lesional ve non-lesional keratinositler arasında farklı eksprese edilen genlerin ilişkili olduğu GO Biological Process terimleri gösterilmektedir. X eksenindeki GeneRatio, ilgili biyolojik sürece atanan genlerin analizde kullanılan toplam genler içerisindeki oranını göstermektedir. Nokta büyüklüğü ilgili süreçle eşleşen gen sayısını, renk ise düzeltilmiş p-değerini, yani p.adjust değerini temsil etmektedir. Daha düşük p.adjust değerine sahip süreçler daha güçlü istatistiksel destek göstermektedir.

Öne çıkan biyolojik süreçler arasında büyümenin düzenlenmesi, apoptotik sürecin pozitif düzenlenmesi, epidermis gelişimi, epidermal hücre farklılaşması, keratinosit farklılaşması ve yara yanıtı bulunmaktadır. Bu süreçlerin zenginleşmesi, lesional ve non-lesional keratinositler arasındaki temel farklılıkların epidermal yenilenme, farklılaşma ve hücre yaşam döngüsü ile ilişkili olduğunu göstermektedir.

Bakır, çinko, kadmiyum ve diğer metal iyonlarına karşı hücresel yanıtlarla ilişkili terimlerin görülmesi, keratinositlerde hücresel stres ve detoksifikasyon mekanizmalarının farklılaşmış olabileceğini düşündürmektedir. Bu terimler genellikle metallotiyonein ailesi gibi ortak stres yanıtı genlerinden kaynaklanabilir. Bu nedenle birden fazla benzer GO teriminin görünmesi, tamamen bağımsız süreçlerin aktive olduğu anlamına gelmeyebilir; aynı gen grubunun farklı fakat ilişkili GO kategorilerinde yer alması söz konusu olabilir.

Genel olarak sonuçlar, keratinositlerin iki koşul arasında yalnızca yapısal farklılaşma bakımından değil, stres yanıtı, apoptoz ve doku onarımı bakımından da farklı transkripsiyonel programlara sahip olduğunu göstermektedir. Ancak bu süreçlerin hangi koşulda arttığını belirlemek için sürece katkı sağlayan genlerin log2 kat değişim yönleri ayrıca incelenmelidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 19 KEGG Lesional vs Nonlesional Keratinocytes

Şekil 19’da keratinositlere ait diferansiyel genlerin KEGG yolaklarıyla ilişkisi gösterilmektedir. Analizde iki temel yolak öne çıkmıştır: Cornified envelope formation **ve** IL-17 signaling pathway**.**

Cornified envelope oluşumu, epidermisin en dış bariyer tabakasının meydana gelmesi ve keratinositlerin terminal farklılaşmasıyla ilişkilidir. Bu yolun zenginleşmesi, lesional ve non-lesional dokular arasında keratinosit olgunlaşması ve epidermal bariyer oluşumu açısından farklılık bulunduğunu göstermektedir. Atopik dermatit gibi epidermal bariyer bozukluğu ile karakterize durumlarda bu yolun farklılaşması biyolojik açıdan anlamlıdır.

IL-17 sinyal yolu, inflamatuvar yanıt, sitokin üretimi ve keratinosit aktivasyonu ile ilişkilidir. Keratinositler yalnızca yapısal epidermal hücreler değildir; aynı zamanda sitokin ve kemokinlere yanıt vererek bağışıklık hücrelerinin dokuya yönlendirilmesine katkı sağlayabilir. Bu nedenle IL-17 sinyal yoluna ait genlerin zenginleşmesi, keratinositlerin inflamatuvar mikroçevreyle etkileşiminin değiştiğini düşündürmektedir.

Bununla birlikte grafikte yalnızca iki yol bulunması ve gen sayılarının görece sınırlı olması nedeniyle sonuçlar dikkatli yorumlanmalıdır. KEGG analizinin gen listesi büyüklüğüne ve kullanılan eşiklere duyarlı olduğu unutulmamalıdır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 20 Go Bp Lesional vs Nonlesional Endothelial Cells

Şekil 20’de endotelyal hücrelere ait farklı eksprese edilen genlerin zenginleştiği biyolojik süreçler gösterilmektedir. En yüksek GeneRatio ve gen sayısına sahip süreçler arasında aktin filament organizasyonu, proteinlerin organellere yönlendirilmesi, protein katlanması, mitokondriyal gen ekspresyonu ve mitokondriyal translasyon yer almaktadır.

Aktin filament organizasyonunun zenginleşmesi, endotelyal hücrelerin hücre şekli, göç, bariyer bütünlüğü ve damar yapısının düzenlenmesiyle ilişkili değişiklikler gösterdiğini düşündürmektedir. Endotelyal hücreler inflamasyon sırasında hücre-hücre bağlantılarını ve geçirgenliklerini değiştirebilir. Bu nedenle sitoskeleton organizasyonundaki farklılıklar, inflamatuvar hücrelerin dokuya geçişi açısından önem taşıyabilir.

Mitokondriyal gen ekspresyonu, mitokondriyal translasyon ve solunumla ilişkili süreçlerin zenginleşmesi, endotelyal hücrelerin enerji üretimi ve metabolik durumunda koşula bağlı farklılıklar bulunduğunu göstermektedir. Protein katlanması ve protein-RNA komplekslerinin düzenlenmesiyle ilişkili terimler ise hücrelerde protein sentezi, stres yanıtı ve RNA işleme mekanizmalarının farklılaştığını düşündürmektedir.

Endotelyal gelişim ve hücre farklılaşmasıyla ilişkili terimlerin bulunması, damar yapısının korunması veya yeniden şekillenmesiyle ilişkili transkripsiyonel değişikliklere işaret etmektedir. Bununla birlikte farklı terimlerin aynı mitokondriyal veya ribozomal genlerden kaynaklanabileceği dikkate alınmalıdır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 21 KEGG Lesional vs Nonlesional Endothelial Cells

Şekil 21’de endotelyal hücrelerde farklılaşan genlerin KEGG yolakları gösterilmektedir. Parkinson hastalığı, Huntington hastalığı, amyotrofik lateral skleroz, prion hastalığı ve nörodejenerasyonla ilişkili yolakların üst sıralarda bulunması, bu hücrelerde söz konusu nörolojik hastalıkların geliştiği anlamına gelmemektedir.

Bu KEGG yolaklarının önemli bir bölümü mitokondriyal elektron taşıma zinciri, oksidatif fosforilasyon, ribozomal proteinler, proteostaz ve hücresel stresle ilişkili ortak genleri içermektedir. Dolayısıyla bu başlıkların birlikte zenginleşmesi, endotelyal hücrelerde esas olarak enerji metabolizması ve protein kalite kontrolüyle ilişkili ortak genlerin farklılaştığını göstermektedir.

Ribozom ve endoplazmik retikulumda protein işlenmesi yolaklarının belirlenmesi, protein sentezi ve katlanma mekanizmalarının değiştiğini desteklemektedir. Fagositoz, antijen işlenmesi ve sunumu, viral enfeksiyon yolakları gibi bağışıklıkla ilişkili başlıkların görülmesi ise endotelyal hücrelerin inflamatuvar sinyallere yanıt verdiğini veya bağışıklık hücreleriyle etkileşiminin değiştiğini düşündürmektedir.

“Fluid shear stress and atherosclerosis” yolunun zenginleşmesi, endotelyal hücrelerin mekanik stres, oksidatif stres ve vasküler sinyallemeyle ilişkili genlerinin farklılaştığını göstermektedir. Bu yolak adı doğrudan ateroskleroz bulunduğunu göstermez; yol içerisindeki ortak endotelyal stres genlerinin zenginleştiğini ifade eder.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 22 Go Bp Lesional vs Nonlesional Fibroblasts

Şekil 22’de fibroblastlarda lesional ve non-lesional koşullar arasında farklılaşan genlerin GO Biological Process sonuçları gösterilmektedir. En belirgin süreçler arasında yaralanmaya yanıt, protein katlanması, oksidatif strese yanıt, ekstraselüler matriks organizasyonu ve hücre dışı yapı organizasyonu yer almaktadır.

Fibroblastlar kollajen ve diğer matriks bileşenlerinin üretiminden sorumlu temel stromal hücrelerdir. Ekstraselüler matriks organizasyonu ve hücre dışı yapı organizasyonunun zenginleşmesi, iki koşul arasında doku mimarisi, kollajen düzenlenmesi ve matriks yeniden şekillenmesi açısından farklılık bulunduğunu göstermektedir.

Yaralanmaya yanıt teriminin yüksek GeneRatio göstermesi, fibroblastların lezyonla ilişkili doku hasarı veya onarım süreçlerinde farklı bir transkripsiyonel duruma geçtiğini düşündürmektedir. Oksidatif strese yanıt ve protein katlanmasıyla ilişkili süreçler, inflamatuvar mikroçevrenin fibroblastlarda hücresel stres oluşturabileceğine işaret etmektedir.

Hücresel solunum, elektron taşıma zinciri, ATP sentezi ve mitokondriyal işlevlerle ilişkili terimlerin birlikte görülmesi, fibroblastların enerji metabolizmasının da koşullar arasında farklılaştığını göstermektedir. Bu sonuçlar fibroblastlarda hem yapısal hem de metabolik yeniden programlanma bulunduğunu düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 23 KEGG Lesional vs Nonlesional Fibroblasts

Şekil 23’te fibroblastlara ait diferansiyel genlerin KEGG yolakları gösterilmektedir. Nörodejenerasyon, amyotrofik lateral skleroz, Alzheimer, Parkinson, Huntington ve prion hastalığı gibi başlıkların üst sıralarda bulunması, çoğunlukla ortak mitokondriyal ve protein homeostazı genlerinin zenginleşmesinden kaynaklanmaktadır.

Endoplazmik retikulumda protein işlenmesi yolunun zenginleşmesi, fibroblastlarda sentezlenen yüksek miktardaki matriks proteinlerinin katlanması ve taşınmasıyla ilişkili olabilir. Fibroblastlar yoğun kollajen ve ekstraselüler matriks proteini ürettikleri için endoplazmik retikulum işlevlerindeki değişiklikler biyolojik olarak önemlidir.

Oksidatif fosforilasyon yolunun belirlenmesi enerji üretim mekanizmalarında farklılık olduğunu; proteazom yolunun belirlenmesi ise hasarlı veya yanlış katlanmış proteinlerin uzaklaştırılmasıyla ilişkili süreçlerin değiştiğini göstermektedir. Sitoskeleton ve kas hücresi yapısı ile ilişkili yollar, fibroblastların aktifleştirilmiş veya kontraktil özellik gösteren bir duruma geçmesiyle ilişkili olabilir.

Bununla birlikte bu sonuçların yönünü belirlemek için ilgili yolaklara katkı sağlayan genlerin lesional veya non-lesional grupta hangi yönde değiştiğinin ayrıca incelenmesi gerekir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 24 Go Bp Lesional vs Nonlesional T Cells

Şekil 24’te T hücrelerinde iki koşul arasında farklılaşan genlerle ilişkili GO süreçleri gösterilmektedir. Grafikte B hücre proliferasyonunun pozitif düzenlenmesi, B hücre aktivasyonu, immün yanıtın düzenlenmesi, lenfosit proliferasyonu ve plazminojen aktivasyonu gibi terimler öne çıkmaktadır.

T hücrelerinde B hücresiyle ilişkili terimlerin bulunması ilk bakışta beklenmedik görünebilir. Ancak GO terimleri çoğu zaman yalnızca tek bir hücre tipine özgü değildir. Sitokinler, kostimülatör moleküller ve ortak lenfosit aktivasyon genleri hem T hücrelerinin işlevinde hem de B hücre yanıtının düzenlenmesinde rol alabilir. Bu nedenle sonuçlar, T hücrelerinin diğer bağışıklık hücrelerini düzenleyen işlevlerinde farklılık bulunduğunu düşündürebilir.

Bununla birlikte grafikte nokta büyüklüklerinin küçük ve Count değerlerinin yaklaşık üç-dört gen düzeyinde olduğu görülmektedir. Az sayıdaki genle elde edilen zenginleştirme sonuçları, tek bir genin eklenmesi veya çıkarılmasıyla önemli ölçüde değişebilir. Dolayısıyla bu sonuçların keşifsel olarak değerlendirilmesi ve gen düzeyinde doğrulanması gerekir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 25 KEGG Lesional vs Nonlesional T Cells

Şekil 25’te T hücrelerine ait KEGG analizinde yalnızca **Salmonella infection** yolunun istatistiksel eşikleri karşıladığı görülmektedir. Yaklaşık yedi genin bu yolakla eşleştiği ve p.adjust değerinin anlamlılık sınırı içerisinde olduğu görülmektedir.

Bu yolak doğrudan örneklerde Salmonella enfeksiyonu bulunduğunu göstermez. KEGG’de “Salmonella infection” başlığı altında NF-κB sinyallemesi, aktin sitoskeletonu, sitokin yanıtı, hücre içi sinyal iletimi ve doğuştan bağışıklık mekanizmalarıyla ilişkili çok sayıda gen bulunmaktadır. Aynı genler farklı inflamatuvar durumlarda da değişebilir.

Yalnızca tek bir yolun anlamlı bulunması, T hücrelerine ait DEG listesinin görece küçük olması veya genlerin sınırlı sayıda KEGG yoluna eşleşmesiyle ilişkili olabilir. Bu nedenle sonuç, genel bir T hücresi inflamatuvar yanıt farklılığına işaret eden sınırlı ve keşifsel bir bulgu olarak değerlendirilmelidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 26 KEGG Lesional vs Nonlesional DC

Şekil 26’da dendritik hücrelere ait diferansiyel genlerin KEGG analizinde yalnızca Salmonella infection yolunun anlamlı bulunduğu görülmektedir. Yaklaşık beş genin bu yolağa katkı sağladığı anlaşılmaktadır.

Dendritik hücreler patojen tanıma, antijen işleme ve bağışıklık yanıtının başlatılmasında önemli rol oynadığından, enfeksiyonla ilişkili KEGG yollarının zenginleşmesi biyolojik olarak şaşırtıcı değildir. Bununla birlikte bu sonuç spesifik bir bakteriyel enfeksiyon kanıtı değildir. Yolağın içerdiği Toll-like reseptör, sitokin, NF-κB veya hücre iskeleti ile ilişkili ortak genlerden kaynaklanabilir.

Gen sayısının düşük olması nedeniyle sonuç dikkatli yorumlanmalıdır. Dendritik hücre sayısının sınırlı olması, diferansiyel gen analizinin gücünü ve zenginleştirme sonucunda elde edilen yol sayısını azaltmış olabilir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 27 Go Bp Lesional vs Nonlesional Tissue Stem Cells

Şekil 27’de doku kök hücreleri olarak anotasyonlanan popülasyonda farklılaşan genlerin GO sonuçları gösterilmektedir. Ekstraselüler matriks organizasyonu, hücre dışı yapı organizasyonu ve dış kapsayıcı yapı organizasyonu en yüksek GeneRatio değerine sahip süreçler arasında yer almaktadır.

Bu sonuçlar, ilgili hücre popülasyonunun doku yapısının korunması, matriks üretimi ve hücre-matriks etkileşimleriyle yakından ilişkili olduğunu göstermektedir. Kas sistemi süreçleri, kas kontraksiyonu, düz kas hücre proliferasyonu ve ossifikasyonla ilişkili terimlerin bulunması, bu popülasyonun mezenkimal veya kontraktil özellikler taşıyabileceğini düşündürmektedir.

Hipoksiye yanıt ve apoptotik sürecin pozitif düzenlenmesi gibi terimler, hücrelerin lesional mikroçevrede oksijen düzeyi, stres ve hücre yaşam döngüsü bakımından farklılaştığını gösterebilir. Ancak “tissue stem cells” anotasyonunun referans tabanlı olduğu dikkate alınmalı; bu hücrelerin gerçekten kök hücre olup olmadığı COL1A1, COL3A1, ACTA2, TAGLN, RGS5 ve benzeri markerlarla ayrıca doğrulanmalıdır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 28 KEGG Lesional vs Nonlesional Tissue Stem Cells

Şekil 28’de doku kök hücresi olarak anotasyonlanan popülasyonda farklılaşan genlerin KEGG sonuçları gösterilmektedir. En yüksek GeneRatio değerine sahip yol “cytoskeleton in muscle cells” olup, bunu fokal adezyon, PI3K-Akt sinyallemesi, integrin sinyallemesi ve aktin sitoskeletonunun düzenlenmesi izlemektedir.

Fokal adezyon, integrin sinyallemesi ve ECM-reseptör etkileşimi hücrelerin ekstraselüler matrikse bağlanması, mekanik sinyalleri algılaması, göç etmesi ve farklılaşması açısından önemlidir. Bu yolların birlikte zenginleşmesi, ilgili hücre popülasyonunda hücre-matriks etkileşiminin temel bir farklılık oluşturduğunu göstermektedir.

PI3K-Akt sinyal yolu hücre yaşamı, proliferasyon ve büyümenin düzenlenmesinde rol oynar. TGF-β sinyal yolu ise fibroblast aktivasyonu, matriks üretimi ve hücre farklılaşmasıyla ilişkilidir. Sonuçlar, lezyonlu ve lezyonsuz dokular arasında bu hücrelerin yapısal yeniden düzenlenme ve doku onarım süreçlerinde farklı roller üstlenebileceğini düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 29 Go Bp Lesional vs Nonlesional Chondrocytes

Şekil 29’da kondrosit olarak anotasyonlanan hücrelerde hücre göçü, hücre hareketliliği ve lokomosyonun negatif düzenlenmesiyle ilişkili GO terimleri öne çıkmaktadır. Düz kas hücresi proliferasyonu ve göçünün negatif düzenlenmesine ilişkin süreçler de zenginleşmiştir.

Bu süreçler, ilgili hücre popülasyonunun hücre hareketliliğini ve mezenkimal hücre davranışlarını düzenleyen genler taşıdığını göstermektedir. Ancak analiz edilen dokunun deri olması nedeniyle “chondrocytes” anotasyonunun biyolojik uygunluğu ayrıca sorgulanmalıdır. Referans tabanlı anotasyon algoritması, ekstraselüler matriks ve kollajen genleri yüksek olan bazı fibroblast veya mezenkimal hücreleri kondrositlere benzetmiş olabilir.

Bu nedenle bu sonuçlar “deride gerçek kondrosit popülasyonu bulunduğu” şeklinde yorumlanmamalıdır. İlgili kümede COL2A1, ACAN, SOX9 ve COMP gibi gerçek kondrosit markerları kontrol edilmeli; bu markerlar bulunmuyorsa popülasyonun fibroblast, perisit veya başka bir mezenkimal hücre tipi olarak yeniden adlandırılması değerlendirilmelidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 30 Go Bp Lesional vs Nonlesional Smooth Muscle Cells

Şekil 30’da düz kas hücresi olarak anotasyonlanan popülasyonda kemotaksi ve bağışıklık hücresi yönlendirilmesiyle ilişkili süreçler öne çıkmaktadır. Kemotaksi, taksis, hücre kemotaksisi, humoral immün yanıt, eozinofil kemotaksisi ve eozinofil göçü en belirgin terimler arasında yer almaktadır.

Eozinofil kemotaksisi ve kemokin aracılı sinyalleme, atopik dermatit gibi tip 2 inflamasyonla ilişkili hastalıklarda biyolojik açıdan önemli olabilir. Düz kas veya perivasküler stromal hücreler, kemokin üreterek inflamatuvar hücrelerin dokuya çekilmesine katkı sağlayabilir. Bununla birlikte bu biyolojik yorum, ilgili sürece katkı sağlayan kemokin genlerinin incelenmesiyle doğrulanmalıdır.

Kompleman aktivasyonu ve T hücresi polaritesinin oluşturulması gibi terimlerin görülmesi, bu popülasyonun yalnızca yapısal bir hücre grubu olmadığını, aynı zamanda bağışıklık mikroçevresinin düzenlenmesinde rol oynayabileceğini düşündürmektedir. Ancak “smooth muscle cells” anotasyonunun ACTA2, TAGLN, MYL9, MYH11 ve RGS5 gibi markerlarla doğrulanması gerekmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 31 KEGG Lesional vs Nonlesional Smooth Muscle Cells

Şekil 31’de düz kas hücrelerine ait diferansiyel genlerin KEGG sonuçları gösterilmektedir. Sitokin-sitokin reseptör etkileşimi en yüksek GeneRatio değerine sahip yolak olarak öne çıkmaktadır. AGE-RAGE sinyal yolu, lipid ve ateroskleroz, kompleman ve koagülasyon kaskatları, viral protein-sitokin etkileşimleri ve TNF sinyal yolu da anlamlı yollar arasında bulunmaktadır.

Sitokin-sitokin reseptör etkileşimi ve TNF sinyal yolunun zenginleşmesi, bu hücrelerde inflamatuvar sinyal alışverişinin koşullar arasında farklılaştığını göstermektedir. Kompleman ve koagülasyon yolakları vasküler yanıt, inflamasyon ve doku hasarıyla ilişkilidir.

AGE-RAGE ve ateroskleroz başlıkları doğrudan diyabetik komplikasyon veya ateroskleroz varlığını göstermez. Bu yollar oksidatif stres, inflamasyon, vasküler sinyalleme ve hücre-matriks etkileşimiyle ilişkili ortak genleri içermektedir. Genel sonuçlar, düz kas veya perivasküler hücrelerin lezyonlu dokudaki inflamatuvar iletişim ağına katkı sağlayabileceğini düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 32 Go Bp Global Lesional vs Nonlesional

Şekil 32’de tüm hücre tiplerinden elde edilen lesional ve non-lesional diferansiyel genler birlikte değerlendirilerek global GO Biological Process analizi yapılmıştır. Aktin filament organizasyonu en yüksek GeneRatio değerine sahip süreç olarak öne çıkmaktadır. Nükleotid metabolizması, oksidatif strese yanıt, protein lokalizasyonu, hücresel solunum ve mitokondriyal gen ekspresyonu da önemli süreçler arasında yer almaktadır.

Aktin ve sitoskeleton organizasyonuyla ilişkili terimler, hücre şekli, göç, adezyon ve doku organizasyonunun iki koşul arasında farklılaştığını göstermektedir. Oksidatif strese yanıt, lezyonlu mikroçevrede reaktif oksijen türleri ve hücresel stres mekanizmalarının rol oynayabileceğini düşündürmektedir.

Hücresel solunum, aerobik solunum, elektron taşıma zinciri, mitokondriyal translasyon ve ATP metabolizmasıyla ilişkili çok sayıda terimin birlikte zenginleşmesi, global farklılığın önemli bir bölümünün enerji metabolizması ve mitokondriyal işlevlerle ilişkili olduğunu göstermektedir.

Bununla birlikte global analiz, tüm hücre tiplerini bir arada değerlendirdiğinden hücre tipi kompozisyonundan etkilenebilir. Örneğin belirli bir hücre tipinin lesional dokuda daha fazla bulunması, o hücre tipine özgü genlerin global DEG listesinde zenginleşmesine yol açabilir. Bu nedenle global sonuçların hücre tipi bazlı analizlerle birlikte yorumlanması gerekmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 33 KEGG Global Lesional vs Nonlesional

Şekil 33’te global diferansiyel gen listesine ait KEGG sonuçları gösterilmektedir. Nörodejenerasyon, Parkinson, amyotrofik lateral skleroz, Alzheimer, Huntington ve prion hastalığı gibi yolak başlıkları üst sıralarda yer almaktadır.

Bu başlıkların birlikte bulunması doğrudan nörodejeneratif bir süreç olduğu anlamına gelmemektedir. Bu yollar mitokondriyal solunum zinciri, ribozomal proteinler, proteazom, oksidatif stres ve apoptozla ilişkili çok sayıda ortak gen içermektedir. Dolayısıyla sonuçlar öncelikle mitokondriyal ve hücresel stres genlerinin farklılaştığını göstermektedir.

Ribozom, oksidatif fosforilasyon, proteazom ve apoptoz yollarının zenginleşmesi; protein sentezi, enerji üretimi, protein yıkımı ve hücre ölümü mekanizmalarının lesional ve non-lesional dokular arasında farklılaştığını desteklemektedir. Fagositoz, antijen işlenmesi ve sunumu, bakteriyel invazyon ve viral enfeksiyon başlıkları ise inflamatuvar ve doğuştan bağışıklık mekanizmalarıyla ilişkili ortak genleri yansıtmaktadır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 34 Lesional Circle

Şekil 34’te lesional dokudaki hücre tipleri arasında CellChat tarafından tahmin edilen ligand-reseptör temelli iletişim ağı gösterilmektedir. Daire çevresindeki her düğüm bir hücre tipini, düğümler arasındaki yönlü bağlantılar ise hücre tipleri arasında tahmin edilen sinyal iletimini temsil etmektedir.

Bağlantıların kalınlığı iletişimin toplam gücüyle ilişkilidir. Kalın bağlantılar, iki hücre tipi arasında daha güçlü veya daha fazla ligand-reseptör etkileşimi bulunduğunu göstermektedir. Hücrenin kendi üzerine dönen bağlantılar ise aynı hücre tipinin otokrin sinyal iletimini temsil etmektedir.

Lesional grupta fibroblast, epitel, endotelyal, kondrosit olarak anotasyonlanan hücreler ve doku kök hücreleri arasında yoğun bağlantılar görülmektedir. Fibroblastların çok sayıda hücre tipiyle bağlantı kurması, ekstraselüler matriks ve büyüme faktörü sinyallerinin lezyonlu mikroçevrede önemli olabileceğini düşündürmektedir.

Ancak CellChat bağlantıları doğrudan deneysel olarak ölçülmüş fiziksel iletişimler değildir. Bunlar, ligand ve reseptör genlerinin ekspresyonuna dayalı hesaplamalı tahminlerdir. Dolayısıyla bağlantılar olası iletişim yollarını göstermekte, kesin hücre-hücre etkileşimi kanıtı oluşturmamaktadır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 35 Lesional Outgoing vs Incoming

Şekil 35’te lesional gruptaki hücre tiplerinin sinyal gönderme ve alma profilleri birlikte gösterilmektedir. Sol taraftaki outgoing heatmap, hücre tiplerinin hangi sinyal yollarında gönderici olarak rol aldığını; sağ taraftaki incoming heatmap ise hangi yollar üzerinden sinyal aldığını göstermektedir.

Heatmap içerisindeki koyu yeşil alanlar daha yüksek göreceli sinyal aktivitesini, açık veya beyaz alanlar ise düşük ya da belirlenemeyen aktiviteyi temsil etmektedir. Sütunlar hücre tiplerini, satırlar ise CellChat tarafından tanımlanan sinyal yollarını göstermektedir.

Kollajen, laminin, fibronektin ve diğer matriks ilişkili yolların fibroblast ve stromal hücrelerde güçlü outgoing aktivite göstermesi, bu hücrelerin lezyonlu dokuda ekstraselüler matriks sinyallerinin temel kaynaklarından biri olduğunu düşündürmektedir. Buna karşılık epitel, endotelyal ve diğer hücrelerin bu yollar için incoming aktivite göstermesi, matriks kaynaklı sinyallere yanıt verdiklerini göstermektedir.

MIF, APP, MK, PTN, galektin ve TGF-β gibi yolların farklı hücre tiplerinde değişen gönderici ve alıcı örüntüler göstermesi, lezyonlu mikroçevrede hücresel iletişimin tek bir yolak üzerinden değil, çok sayıda paralel sinyal sistemi üzerinden gerçekleştiğini göstermektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 36 Nonlesional Circle

Şekil 36’da non-lesional örneklerde tahmin edilen hücre-hücre iletişim ağı gösterilmektedir. Lesional grupta olduğu gibi düğümler hücre tiplerini, bağlantılar ise hücreler arasında öngörülen ligand-reseptör iletişimlerini temsil etmektedir.

Non-lesional dokuda da fibroblast, epitel, endotelyal ve stromal hücreler arasında çok sayıda iletişim bağlantısı bulunmaktadır. Bu durum, normal görünümlü dokuda da hücrelerin doku homeostazını korumak amacıyla yoğun bir iletişim ağı oluşturduğunu göstermektedir.

Lesional ve non-lesional ağlar görsel olarak karşılaştırıldığında bazı bağlantıların kalınlığında ve belirli hücrelerin ağ içerisindeki rolünde farklılıklar görülmektedir. Ancak daire grafiklerinin ölçekleri aynı değilse yalnızca çizgi kalınlıklarına bakılarak doğrudan nicel karşılaştırma yapılmamalıdır. Koşullar arasındaki toplam etkileşim sayısı ve gücü Şekil 38’deki karşılaştırmalı grafikle daha güvenilir biçimde değerlendirilmiştir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 37 Nonlesional Outgoing vs Incoming

Şekil 37’de non-lesional dokudaki hücre tiplerinin sinyal gönderme ve alma profilleri gösterilmektedir. Outgoing panel hücrelerin hangi sinyal yollarının kaynağı olduğunu, incoming panel ise hangi sinyalleri aldığını ortaya koymaktadır.

Fibroblastların kollajen, laminin ve diğer matriks ilişkili sinyallerde belirgin gönderici özellik gösterdiği görülmektedir. Bu durum fibroblastların normal doku yapısının korunması ve hücrelerin matrikse bağlanmasının düzenlenmesindeki rolüyle uyumludur.

Farklı hücre tiplerinin incoming profilleri, normal dokudaki hücrelerin matriks, büyüme faktörü ve immün düzenleyici sinyallere farklı düzeylerde yanıt verdiğini göstermektedir. Lesional gruptaki heatmap ile karşılaştırıldığında bazı yolların hücre tipi dağılımında değişiklikler bulunduğu görülmektedir. Bu farklılıklar, lezyon oluşumu sırasında hangi hücrelerin sinyal kaynağı veya hedefi hâline geldiğini anlamak açısından önemlidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 38 Compare Inferred Interactions and Interaction Strength

Şekil 38’de non-lesional ve lesional gruplar toplam hücre-hücre etkileşim sayısı ve toplam iletişim gücü bakımından karşılaştırılmıştır. Sol panelde tahmin edilen etkileşim sayısı, sağ panelde ise bu etkileşimlerin birleştirilmiş iletişim gücü gösterilmektedir.

Non-lesional grupta yaklaşık **421**, lesional grupta ise yaklaşık **594** hücre-hücre etkileşimi belirlenmiştir. Bu sonuç lesional dokuda tahmin edilen ligand-reseptör bağlantılarının sayısının daha fazla olduğunu göstermektedir.

Toplam iletişim gücü non-lesional grupta yaklaşık **37,4**, lesional grupta ise yaklaşık **54,9** olarak hesaplanmıştır. Dolayısıyla lesional grupta yalnızca etkileşim sayısı değil, tahmin edilen toplam sinyal gücü de daha yüksektir.

Bu bulgu, lezyonlu dokuda hücrelerin daha yoğun ve karmaşık bir iletişim ağı oluşturduğunu düşündürmektedir. Artan iletişim; inflamatuvar sitokinler, kemokinler, büyüme faktörleri ve ekstraselüler matriks sinyallerindeki değişikliklerle ilişkili olabilir. Bununla birlikte CellChat sonuçları hücre sayısından ve hücre tipi oranlarından etkilenebilir. Bu nedenle gruplardaki hücre sayılarının dengesi ve kullanılan örnek sayıları dikkate alınarak yorum yapılmalıdır.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 39 Metascape Global Lesional

Şekil 39’da lesional grupta yüksek ekspresyon gösteren global DEG listesinin Metascape analizi sonucunda elde edilen zenginleştirilmiş biyolojik terimler gösterilmektedir. Yatay çubuklar zenginleşen terimleri, x ekseni ise terimlerin istatistiksel anlamlılığını gösteren −log10(P) değerini temsil etmektedir. Çubuk uzadıkça p-değeri küçülmekte ve ilgili terimin istatistiksel desteği artmaktadır.

Translasyon ve ribozomal işlevlerle ilişkili süreçlerin üst sıralarda bulunması, lesional dokuda protein senteziyle ilişkili genlerin yüksek düzeyde temsil edildiğini göstermektedir. Mitokondriyal gen ekspresyonu, enerji üretimi ve metabolik süreçlerle ilişkili terimler, lezyonlu dokuda hücresel enerji gereksiniminin ve metabolik aktivitenin değiştiğini düşündürmektedir.

Protein katlanması, endoplazmik retikulum işlevleri ve protein transportuyla ilişkili süreçler, artan protein üretimine eşlik eden proteostaz gereksinimini yansıtabilir. Antijen işlenmesi, sitokin üretimi ve lenfosit aktivasyonuyla ilişkili terimlerin bulunması ise lezyonlu dokudaki bağışıklık aktivasyonunu desteklemektedir.

Çok sayıda terimin birbiriyle ilişkili olması nedeniyle listedeki her başlık bağımsız bir biyolojik mekanizma olarak değerlendirilmemelidir. Örneğin translasyon, ribozom ve protein senteziyle ilişkili farklı terimler büyük ölçüde aynı gen grubundan kaynaklanabilir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 40 Global Lesional upregulated DEG listesine ait Metascape enrichment analizi

Şekil 40’ta lesional grupta yüksek eksprese edilen genlere ait zenginleştirilmiş terimler, aralarındaki benzerliğe göre bir ağ içerisinde gösterilmiştir. Her düğüm bir zenginleştirilmiş terimi, düğümler arasındaki çizgiler ise terimlerin ortak genler paylaşmasını veya işlevsel olarak benzer olmasını temsil etmektedir.

Aynı renkte gösterilen düğümler aynı biyolojik kümeye atanmıştır. Büyük ve yoğun bağlantılı kümeler, çok sayıda ilişkili GO veya yolak teriminin ortak bir biyolojik tema altında toplandığını göstermektedir. Translasyon, ribozomal işlevler, mitokondriyal gen ekspresyonu, protein işlenmesi ve immün yanıtla ilişkili ayrı fakat bağlantılı kümeler görülmektedir.

Ağın merkezinde bulunan ve çok sayıda bağlantıya sahip terimler, diğer süreçlerle daha fazla gen paylaşan geniş biyolojik programları temsil eder. Buna karşılık ağın kenarında bulunan küçük kümeler daha spesifik veya daha az sayıda gen tarafından desteklenen süreçleri gösterebilir.

Bu analiz, lesional dokuda gözlenen değişimlerin tek tek bağımsız yolaklardan ziyade, birbirleriyle bağlantılı metabolik, translasyonel ve immün programlardan oluştuğunu göstermektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 41 Global Lesional colored by p-value, where terms containing more genes tend to have a more significant p-value.

Şekil 41’de Şekil 40’taki aynı zenginleştirme ağı bu kez terimlerin p-değerlerine göre renklendirilmiştir. Açık renkler daha yüksek, koyu turuncu ve kahverengi renkler ise daha düşük p-değerlerini ve daha güçlü istatistiksel anlamlılığı temsil etmektedir.

Ağın bazı merkezi bölgelerinde koyu renkli düğümlerin yoğunlaşması, bu bölgelerdeki terimlerin hem birbirleriyle bağlantılı hem de istatistiksel olarak güçlü olduğunu göstermektedir. Translasyon, mitokondriyal işlev ve protein işlenmesiyle ilişkili kümelerin yüksek anlamlılık göstermesi, lesional gruptaki temel transkripsiyonel değişimlerin bu süreçlerde yoğunlaştığını desteklemektedir.

Bu grafikte düğüm renginin biyolojik artış miktarını değil, yalnızca istatistiksel anlamlılığı gösterdiği unutulmamalıdır. Çok düşük p-değeri, etkinin mutlaka biyolojik olarak çok büyük olduğu anlamına gelmez; ilgili genlerin süreç içerisinde beklenenden daha sık bulunmasının istatistiksel olarak güçlü olduğunu gösterir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 42 Metascape Nonlesional

Şekil 42’de non-lesional grupta daha yüksek ekspresyon gösteren global genlerin Metascape zenginleştirme sonuçları sunulmuştur. Epidermis gelişimi, endotelyal gelişim, ekstraselüler matriks organizasyonu, hücre farklılaşması ve doku yapısının korunmasıyla ilişkili süreçler üst sıralarda yer almaktadır.

Epidermis gelişimi ve epidermal farklılaşmayla ilişkili terimlerin zenginleşmesi, non-lesional dokuda normal epidermal yapı ve farklılaşma programlarının daha belirgin biçimde temsil edildiğini düşündürmektedir. Ekstraselüler matriks organizasyonu ve hücre-matriks etkileşimiyle ilişkili terimler, doku mimarisinin ve stromal yapının korunmasıyla bağlantılıdır.

Dolaşım sistemi gelişimi, vasküler gelişim ve endotelyal farklılaşma gibi terimler, non-lesional dokuda damar yapısının korunmasıyla ilişkili genlerin daha yüksek temsil edildiğini gösterebilir. Yara iyileşmesi, hücre göçü ve kas yapısının gelişimiyle ilişkili süreçler ise stromal ve mezenkimal hücrelerin normal doku homeostazındaki rollerini yansıtmaktadır.

Bu sonuçlar, lesional grupta metabolik ve inflamatuvar programlar öne çıkarken non-lesional grupta doku organizasyonu, farklılaşma ve homeostazla ilişkili programların daha belirgin olduğunu düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 43 Global Nonlesional upregulated DEG listesine ait Metascape enrichment analizi

Şekil 43’te non-lesional grupta yüksek eksprese edilen genlerle ilişkili biyolojik terimler, benzerliklerine göre ağ yapısında gösterilmiştir. Aynı renkteki düğümler ortak bir biyolojik temaya ait terimleri temsil etmektedir.

Epidermal gelişim, hücre farklılaşması, ekstraselüler matriks organizasyonu, hücre hareketliliği, dolaşım sistemi gelişimi, kas yapısının gelişimi, homeostaz ve vasküler gelişimle ilişkili kümelerin birbirleriyle bağlantılı olduğu görülmektedir.

Ağın geniş ve yoğun bağlantılı olması, non-lesional dokuda doku yapısının korunmasıyla ilişkili süreçlerin birbirinden bağımsız olmadığını göstermektedir. Örneğin epidermal farklılaşma, matriks organizasyonu ve hücre adezyonu aynı doku homeostazı programının farklı yönlerini temsil edebilir.

Lesional ağdan farklı olarak non-lesional ağda yapısal gelişim ve doku organizasyonu temalarının daha belirgin olması, sağlıklı görünümlü dokuda normal farklılaşma programlarının korunabildiğini düşündürmektedir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 44 Global Nonlesional colored by p-value, where terms containing more genes tend to have a more significant p-value.

Şekil 44’te non-lesional gruptaki zenginleştirilmiş terimler p-değerlerine göre renklendirilmiştir. Daha koyu renkli düğümler daha düşük p-değerine ve daha yüksek istatistiksel anlamlılığa sahiptir.

Ağın merkezindeki yoğun bağlantılı ve koyu renkli düğümler, non-lesional grupta en güçlü biçimde temsil edilen biyolojik programları göstermektedir. Epidermal gelişim, ekstraselüler matriks organizasyonu, hücre farklılaşması ve doku yapısının korunmasıyla ilişkili bölgelerin belirgin olması, non-lesional dokuya özgü homeostatik özellikleri desteklemektedir.

Şekil 43 biyolojik terimlerin hangi temalar altında kümelendiğini gösterirken, Şekil 44 bu terimlerden hangilerinin daha güçlü istatistiksel kanıta sahip olduğunu göstermektedir. Dolayısıyla iki grafik birlikte değerlendirilmelidir.

*[Görsel: Orijinal Word dosyasında bulunmaktadır.]*

Şekil 45 Global Nonlesional Protein-protein Interaction Enrichment Analysis

Şekil 45’te non-lesional grupta yüksek eksprese edilen genlerden oluşturulan protein-protein etkileşim ağı ve bu ağ içerisindeki yoğun protein modülleri gösterilmektedir. Soldaki geniş ağda her düğüm bir proteini, çizgiler ise proteinler arasındaki bilinen veya tahmin edilen etkileşimleri temsil etmektedir.

Ağ içerisinde çok sayıda bağlantıya sahip merkezi proteinler, birden fazla biyolojik süreçle ilişkili düzenleyici veya yapısal proteinler olabilir. Kırmızıyla işaretlenen düğümler ağın yüksek bağlantılı bölgelerinde yer alan proteinleri göstermektedir. Bu proteinler, non-lesional dokuya ait biyolojik programların sürdürülmesinde merkezi rol oynayan adaylar olarak değerlendirilebilir.

Sağ tarafta gösterilen MCODE sonuçları, geniş protein ağının içerisinde birbirleriyle yoğun biçimde etkileşen alt modülleri ortaya koymaktadır. Her MCODE kümesi, ortak bir protein kompleksi, sinyal yolu veya biyolojik işlevle ilişkili olabilir. Örneğin bir modül ekstraselüler matriks proteinlerini, başka bir modül hücre iskeleti veya protein senteziyle ilişkili proteinleri içerebilir.

MCODE analizi, tek tek genlerden ziyade birlikte çalışan protein gruplarını belirlemesi açısından önemlidir. Bununla birlikte MCODE kümelerinin biyolojik anlamı, her modülde bulunan genlerin ayrıca GO ve yolak analizine alınmasıyla doğrulanmalıdır.

## Kaynakça

Afshari, M., Kolackova, M., Rosecka, M., Čelakovská, J., & Krejsek, J. (2024). Unraveling the skin; a comprehensive review of atopic dermatitis, current understanding, and approaches. *Frontiers in Immunology*, *15*, 1361005. https://doi.org/10.3389/fimmu.2024.1361005

American Academy of Dermatology Association. (2025). *Atopic dermatitis clinical guideline*. (n.d.). https://www.aad.org/member/clinical-quality/guidelines/atopic-dermatitis

Carr, S., Pratt, R., White, F., & Watson, W. (2024). Atopic dermatitis. *Allergy Asthma and Clinical Immunology*, *20*(S3), 63. https://doi.org/10.1186/s13223-024-00927-2

National Institute for Health and Care Excellence (NICE). (2025, September 22). *Atopic eczema in under 12s: diagnosis and management*. NCBI Bookshelf. https://www.ncbi.nlm.nih.gov/books/NBK571009/

Yue, C., Zhou, H., Wang, X., Yu, J., Hu, Y., Zhou, P., Zhao, F., Zeng, F., Li, G., Li, Y., Feng, Y., Sun, X., Huang, S., He, M., Wu, W., Huang, N., & Li, J. (2024). Atopic dermatitis: pathogenesis and therapeutic intervention. *MedComm*, *5*(12), e70029. https://doi.org/10.1002/mco2.70029
