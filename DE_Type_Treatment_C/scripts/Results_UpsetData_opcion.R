mybiglist <- list()
for(i in 1:5){
  a <- runif(10)
  b <- rnorm(16)
  c <- rbinom(8, 5, i/10)
  name <- paste('item:',i,sep='')
  tmp <- list(uniform=a, normal=b, binomial=c)
  mybiglist[[name]] <- tmp
}

# https://stackoverflow.com/questions/12511648/building-a-list-in-a-loop-in-r-getting-item-names-correct


# Funcion 1
plot_UpsetR <- function(indir, abslogFC = 1){ 
  results <- list() # create a empty list
  # --- Load data -----
  # change directory
  setwd(indir)
  # change to place where "DE_plotName.csv" files are:
  files <- dir(getwd(), pattern = "^DE_(.+)\\.csv$") #only CSV
  
  # --- multiple files ----
  for (f in files){
    print(f)
    plot_name <- gsub("^DE_(.+)\\.csv$", "\\1", f) #name
    
    # Load data
    df <- read.csv(file = f, row.names = 'X')
    df <- na.omit(df)
    
    # --
    # 1. Add expression column
    df$genes <- row.names(df) # Gene names
    df <- df %>% 
      dplyr::mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
                                           log2FoldChange <= -(abslogFC) & padj < 0.05 ~ "Down-regulated",
                                           TRUE ~ "Unchanged")) 
    
    # Up names
    up_names <- df %>%
      filter(Expression == "Up-regulated")
    up_names <- up_names$genes # only names
    
    # Down names
    down_names <- df %>%
      filter(Expression == "Down-regulated")
    down_names <- down_names$genes # only names
    
    # Results
    # tmp <- list(UP = up_names, DOWN = down_names); results[[plot_name]] <- tmp
    # results[[paste(plot_name, "UP", sep="_")]] <- list(UP = up_names)
    # results[[paste(plot_name, "DOWN", sep="_")]] <- list(DOWN = down_names)
    name1 = paste(plot_name, "UP", sep="_")
    print(name1)
    name2 = paste(plot_name, "DOWN", sep="_")
    print(name2)
    
    # add information
    tmp <- list( name1 = up_names)
    names(tmp) <- name1 #renames
    results <- c(results, tmp)
    tmp <- list( name2 = down_names)
    names(tmp) <- name2 #renames
    results <- c(results, tmp)
    
  } #end loop
  return(results)
  
} # end function


# Opcion A
# si en resultados cambiamos
tmp <- list(UP = up_names, DOWN = down_names); 
results[[plot_name]] <- tmp

# List of 4
# $ control_moVSmoDC    :List of 2
# ..$ UP  : chr [1:3199] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# ..$ DOWN: chr [1:3584] "RBM3" "BIN2" "VNN2" "EIF2S3" ...
# $ control_moVSmoDCIMQ :List of 2
# ..$ UP  : chr [1:3153] "MRC1" "SLC7A8" "PPAP2B" "GPNMB" ...
# ..$ DOWN: chr [1:3604] "RBM3" "BIN2" "VNN2" "EIF2S3" ...
# $ control_moVStolDC   :List of 2
# ..$ UP  : chr [1:3276] "MRC1" "SLC7A8" "MAOA" "GPNMB" ...
# ..$ DOWN: chr [1:3487] "BIN2" "RBM3" "EIF2S3" "VNN2" ...
# $ control_moVStolDCIMQ:List of 2
# ..$ UP  : chr [1:3099] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# ..$ DOWN: chr [1:3589] "RBM3" "BIN2" "EIF2S3" "VNN2" ...


# Opcion B

results[[paste(plot_name, "UP", sep="_")]] <- list(UP = up_names)
results[[paste(plot_name, "DOWN", sep="_")]] <- list(DOWN = down_names)
# List of 8
# $ control_moVSmoDC_UP      :List of 1
# ..$ UP: chr [1:3199] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# $ control_moVSmoDC_DOWN    :List of 1
# ..$ DOWN: chr [1:3199] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# $ control_moVSmoDCIMQ_UP   :List of 1
# ..$ UP: chr [1:3153] "MRC1" "SLC7A8" "PPAP2B" "GPNMB" ...
# $ control_moVSmoDCIMQ_DOWN :List of 1
# ..$ DOWN: chr [1:3153] "MRC1" "SLC7A8" "PPAP2B" "GPNMB" ...
# $ control_moVStolDC_UP     :List of 1
# ..$ UP: chr [1:3276] "MRC1" "SLC7A8" "MAOA" "GPNMB" ...
# $ control_moVStolDC_DOWN   :List of 1
# ..$ DOWN: chr [1:3276] "MRC1" "SLC7A8" "MAOA" "GPNMB" ...
# $ control_moVStolDCIMQ_UP  :List of 1
# ..$ UP: chr [1:3099] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# $ control_moVStolDCIMQ_DOWN:List of 1
# ..$ DOWN: chr [1:3099] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...


# Opcion C
# Results
name1 = paste(plot_name, "UP", sep="_")
print(name1)
name2 = paste(plot_name, "DOWN", sep="_")
print(name2)

# add information
tmp <- list( name1 = up_names)
names(tmp) <- name1 #renames
results <- c(results, tmp)
tmp <- list( name2 = down_names)
names(tmp) <- name2 #renames
results <- c(results, tmp)

# List of 8
# $ control_moVSmoDC_UP      : chr [1:3199] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# $ control_moVSmoDC_DOWN    : chr [1:3584] "RBM3" "BIN2" "VNN2" "EIF2S3" ...
# $ control_moVSmoDCIMQ_UP   : chr [1:3153] "MRC1" "SLC7A8" "PPAP2B" "GPNMB" ...
# $ control_moVSmoDCIMQ_DOWN : chr [1:3604] "RBM3" "BIN2" "VNN2" "EIF2S3" ...
# $ control_moVStolDC_UP     : chr [1:3276] "MRC1" "SLC7A8" "MAOA" "GPNMB" ...
# $ control_moVStolDC_DOWN   : chr [1:3487] "BIN2" "RBM3" "EIF2S3" "VNN2" ...
# $ control_moVStolDCIMQ_UP  : chr [1:3099] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# $ control_moVStolDCIMQ_DOWN: chr [1:3589] "RBM3" "BIN2" "EIF2S3" "VNN2" ...