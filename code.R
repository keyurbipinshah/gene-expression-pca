# Loading libraries
library(ggplot2)
library(ggrepel)
library(scales)

# Setting R options
windowsFonts(tw_cen_mt = windowsFont("Tw Cen MT"))

# Loading the data
load("data.RData")
ls() # To figure out name of the object in which data is stored

# Checking the dimensions of the data, names of samples and names of gene
# expressions
dim(leuk)
rownames(leuk)
colnames(leuk)

# Classes of each column to figure out column which holds the sample labels
col_classes <- sapply(leuk, class)
unique(col_classes)
Filter(function(x) x == "factor", col_classes)

# Extracting the labels in a vector
labs <- leuk[["V5001"]]

# Extracting all columns but the labels in a matrix
gene_exp <- as.matrix(leuk[, names(leuk) != "V5001"])

# Variance-covariance matrix of gene expressions
sigma <- cov(gene_exp)

# Eigen decomposition of the variance-covariance matrix
eigen_decomp <- eigen(sigma)

# Extracting eigenvalues and eigenvectors
eigen_vals <- eigen_decomp$values
length(eigen_vals)

eigen_vecs <- eigen_decomp$vectors
dim(eigen_vecs)

# Proportion of variance explained by PCs
prop_explained <- cumsum(eigen_vals) / sum(eigen_vals)
prop_explained[1:10]
max(which(prop_explained <= 0.9))

# Scoring the observations on the first two PCs
scores <- scale(gene_exp, center = TRUE, scale = FALSE) %*% eigen_vecs[, 1:2]
scores <- as.data.frame(scores)
scores$labels <- labs
scores$sample_num <- rownames(scores)

# Initial Bi-plot
ggplot2::ggplot(scores, ggplot2::aes(x = -V1, y = V2, colour = labels)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_minimal(base_family = "tw_cen_mt") +
  ggplot2::labs(x = "Principal Component 1", y = "Principal Component 2") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggrepel::geom_text_repel(ggplot2::aes(x = -V1, y = V2, label = sample_num),
                           inherit.aes = FALSE)

# Recreating the bi-plot with (possibly) mislabelled sample IDs
score_labs <- subset(scores, sample_num %in% c("19", "10", "2", "35"))
ggplot2::ggplot(scores, ggplot2::aes(x = -V1, y = V2, colour = labels)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_minimal(base_family = "tw_cen_mt") +
  ggplot2::labs(x = "Principal Component 1", y = "Principal Component 2",
                color = "Leukemia Type") +
  ggplot2::scale_color_manual(values = c("#01B8AA", "#374649", "#FD625E"),
                              labels = c("AML" = "Acute Myeloid",
                                         "ALL-B" = "Acute Lymphoblastic B-Cell",
                                         "ALL-T" = "Acute Lymphoblastic T-Cell")) +
  ggplot2::geom_point(data = score_labs,
                      ggplot2::aes(x = -V1, y = V2), colour = "black", size = 5, shape = 1) +
  ggrepel::geom_text_repel(ggplot2::aes(x = -V1, y = V2, label = sample_num),
                           data = score_labs,
                           inherit.aes = FALSE,
                           point.padding = 0.8) +
  ggplot2::scale_x_continuous(labels = scales::comma_format()) +
  ggplot2::scale_y_continuous(labels = scales::comma_format()) +
  ggplot2::theme(legend.position = c(0.85, 0.18),
                 legend.title.align = 0.5)

ggplot2::ggsave("Biplot.png", device = "png")