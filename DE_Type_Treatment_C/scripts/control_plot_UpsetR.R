# UpSetR - plots
# Author: Evelia Coss
# Date: 11 Dec 2023
# Make your own UpSetR
# ---

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---- Libraries ----
library(dplyr)  # Manipulacion de datos
library(ggplot2) # Graficas
library(UpSetR) # UpsetR
library(ComplexHeatmap) # Heatmap and UpsetR)


#-----load data and function ---------

# LAVIS sever
#indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
# My computer
indir <- "C:/Users/ecoss/OneDrive - CINVESTAV/Documentos/Posdoc_LIIGH/LUPUS_Proyecto_2023/bulkRNA_lupus2023/DE_Type_Treatment_C/DE_files/"
load(paste0(indir, 'DE_control_results.RData'))

# figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/figures/'
figdir <- 'C:/Users/ecoss/OneDrive - CINVESTAV/Documentos/Posdoc_LIIGH/LUPUS_Proyecto_2023/bulkRNA_lupus2023/DE_Type_Treatment_C/figures/'

#load function
functiondir <- 'C:/Users/ecoss/OneDrive - CINVESTAV/Documentos/Posdoc_LIIGH/LUPUS_Proyecto_2023/bulkRNA_lupus2023/'
source(paste(functiondir, "UpsetR_plot_function.R"))

#-----USAGE---------
# Run Function
UpsetR_data <- ConvertData_UpsetR(indir = indir)
str(UpsetR_data)

# ---- Graph with UpSetR ----------
png(file= paste0(figdir, "UpsetR_UpSetR.png"), width = 18, height = 10, units = 'in', res = 600)

upset(fromList(UpsetR_data), # Convert Dataset to matrix
      sets = names(UpsetR_data), # Dataset to show
      order.by = "freq", # ordenar por la frecuencia
      #boxplot.summary = 2,
      # colors
      sets.bar.color = c(rep("dodgerblue3",4), rep("firebrick3",4)),
      shade.color = c("blue", "blue", "red", "red"), # color de las interacciones, el color se pone de abajo hacia arriba
      set_size.show =TRUE, # mostrar los valores en el set
      set_size.scale_max = 4000, #aumentar escala del set
      text.scale = 2, # Tamano de los numeros
      line.size = 1,
      # labels
      mainbar.y.label = "Intersection Genes", # label histograma (arriba)
      sets.x.label= "Number of Genes" # label set (abajo)
)

dev.off()

# ---- Graph with Complexheatmap (upset)----------

# dataset
m1 <- ComplexHeatmap::make_comb_mat(UpsetR_data, mode = "distinct") # the default mode is `distinct`
m1 <- m1[comb_degree(m1) > 0] # eliminar datos sin interacciones

# A combination matrix with 8 sets and 48 combinations.
# ranges of combination set size: c(3, 3604).
# mode for the combination size: intersect.
# sets are on rows.
# 
# Top 8 combination sets are:
#   control_moVSmoDC_UP control_moVSmoDC_DOWN control_moVSmoDCIMQ_UP control_moVSmoDCIMQ_DOWN control_moVStolDC_UP control_moVStolDC_DOWN control_moVStolDCIMQ_UP control_moVStolDCIMQ_DOWN     code size
# x                                                                                               00010000 3604
# x 00000001 3589
# x                                                                                                                                               01000000 3584
# x                                                   00000100 3487
# x                                                                          00001000 3276
# x                                               x                                                                                               01010000 3234
# x                                                                                                                                                                     10000000 3199
# x                                                                                                                        00100000 3153
# 
# Sets are:
#   set size
# control_moVSmoDC_UP 3199
# control_moVSmoDC_DOWN 3584
# control_moVSmoDCIMQ_UP 3153
# control_moVSmoDCIMQ_DOWN 3604
# control_moVStolDC_UP 3276
# control_moVStolDC_DOWN 3487
# control_moVStolDCIMQ_UP 3099
# control_moVStolDCIMQ_DOWN 3589

set_name(m1)
set_size(m1) # Number of genes across datastes
comb_size(m1) # Number of intersect
comb_degree(m1) # Grupos y cantidad de interacciones (1 a 4 interacciones)

ss <- c("control_moVSmoDC_UP", "control_moVSmoDCIMQ_UP", "control_moVStolDC_UP", "control_moVStolDCIMQ_UP",
                 "control_moVSmoDC_DOWN", "control_moVSmoDCIMQ_DOWN", "control_moVStolDC_DOWN", "control_moVStolDCIMQ_DOWN")
cs <- comb_size(m1) # Number of intersect


# Grupos
subgroup <- c("control_moVSmoDC_UP" = "moDC",
             "control_moVSmoDCIMQ_UP" = "moDC",
             "control_moVStolDC_UP" = "tolDC",
             "control_moVStolDCIMQ_UP" = "tolDC",
             "control_moVSmoDC_DOWN" = "moDC",
             "control_moVSmoDCIMQ_DOWN" = "moDC",
             "control_moVStolDC_DOWN" = "tolDC",
             "control_moVStolDCIMQ_DOWN" = "tolDC"
)

expressed <- c("control_moVSmoDC_UP" = "UP",
              "control_moVSmoDCIMQ_UP" = "UP",
              "control_moVStolDC_UP" = "UP",
              "control_moVStolDCIMQ_UP" = "UP",
              "control_moVSmoDC_DOWN" = "DOWN",
              "control_moVSmoDCIMQ_DOWN" = "DOWN",
              "control_moVStolDC_DOWN" = "DOWN",
              "control_moVStolDCIMQ_DOWN" = "DOWN"
)


# orden de las combinaciones
combOrder <- order(comb_size(m1), decreasing = TRUE) #ordenar interacciones de mayor a menor
deletePo <- as.vector(which(comb_degree(m1) == 1))
combOrder <- combOrder[- which(combOrder %in% deletePo)] #eliminar valores con genes unicos
combOrder <- c(combOrder, which(comb_degree(m1) == 1)) # agregar las posiciones al final  
combOrder <-as.vector(combOrder)

# Grafica
png(file= paste0(figdir, "UpsetR_complexheatmap.png"), width = 12, height = 5, units = 'in', res = 600)

UpSet(m1,
      set_order = ss, # orden en los nombres de los sets
      bg_col = c("#F0F0F0", "#CCCCCC"), # relleno de filas (de abajo hacia arriba)
      bg_pt_col = "#bdbdbd", # color de los puntos
      comb_col = c("firebrick3", "black", "black", "dodgerblue3", rep("black", length(comb_degree(m1))-4)), 
      comb_order = combOrder, #ordenar interacciones de mayor a menor
      top_annotation = HeatmapAnnotation(
        degree = as.character(comb_degree(m1)), # Grupos por funciones
        "Genes \nIntersections" = anno_barplot(cs, 
                                             ylim = c(0, max(cs)*1.1),
                                             border = FALSE, 
                                             gp = gpar(fill = "black"), # color de relleno
                                             height = unit(4, "cm"), 
                                             add_numbers = TRUE # Agregar numeracion en barplot
        ), # end anno_barplot
         
        annotation_name_rot = 0, # rotar de la leyenda
        annotation_name_side = "left" # posicion de la leyenda
        ), # end HeatmapAnnotation
      # Anotaciones a la izquierda
      right_annotation = rowAnnotation(cellType = subgroup[set_name(m1)], # Grupos por tipo celular
                                       Expression = expressed[set_name(m1)], # Grupos por expresion
                                       col = list(
                                         cellType = c('moDC' = '#ae017e', 'tolDC' = '#49006a'),
                                         Expression = c('UP' = "firebrick3", 'DOWN' = "dodgerblue3")
                                         ), #end colors
                                      show_annotation_name = FALSE),
      left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
)# end Upset

dev.off()
# Manual https://rdrr.io/bioc/ComplexHeatmap/man/UpSet.html
