library(readxl)
library(openxlsx)

DoBMIQ <- function(probes_path, beta_path) {

probes <- read_excel(probes_path)
print(paste("probes_path:", probes_path))
beta_data <- read_excel(beta_path)
print(paste("beta_path:", beta_path))

# Creating data matrix
probe_ids <- sub("_.*", "", beta_data$`probe set`)
# probe_ids <- sub(".*_", "", beta_data$`probe set`)
data.m <- as.matrix(beta_data[, -1])
sample_id <- colnames(beta_data)[-1]
rownames(data.m) <- probe_ids
# save(data.m, file = "data.m.RData")

# Remove rows with any NA values
complete_rows <- complete.cases(data.m)
data.m <- data.m[complete_rows, ]

# Update probe_ids to match the filtered data.m
probe_ids <- probe_ids[complete_rows]
rownames(data.m) <- probe_ids
print(head(data.m))             
print(data.m[1:10, 1:5])

print(paste("Removed", sum(!complete_rows), "rows containing NA values"))
print(paste("Final dimensions of data.m:", nrow(data.m), "rows by", ncol(data.m), "columns"))

### DoBMIQ.R
source("./R Scripts/BMIQ_1.4.R")

print("Read Success.")


# --- Match probe_ids to probe annotation ---
index <- match(probe_ids, probes$name)
matched_probes <- probes[index, ]

# Extract targetid
probe_targetid <- matched_probes$targetid

# Extract probe type: last two characters of targetid (e.g. "11", "21")
type_ids <- substr(probe_targetid, nchar(probe_targetid) - 1, nchar(probe_targetid))
type_ids <- as.numeric(type_ids)

# Build design vector and indices
design.v <- ifelse(type_ids == 11, 1,
                   ifelse(type_ids == 21, 2, NA))

type1.idx <- which(type_ids == 11)
type2.idx <- which(type_ids == 21)

# Debug print
# cat("Probe type distribution:\n")
# print(table(type_ids, useNA = "ifany"))
# cat("Total Type 1 probes:", length(type1.idx), "\n")
# cat("Total Type 2 probes:", length(type2.idx), "\n")


pdf("Profiles.pdf", width=4, height=3)
for(s in 1:ncol(data.m)){
  check_for_fail <- grepl("FAILED", sample_id[s])
  if (check_for_fail) {
    next
  }

  cat("Sample:", sample_id[s], "\n")
  cat("Length of type1.idx:", length(type1.idx), "\n")
  cat("Non-NA values in data.m[type1.idx, s]:", sum(!is.na(data.m[type1.idx, s])), "\n")
  cat("Values:", data.m[type1.idx, s][!is.na(data.m[type1.idx, s])], "\n\n")

  # Now check if there's enough data to continue
  if (sum(!is.na(data.m[type1.idx, s])) < 2) {
    cat("⚠️ Skipping sample", sample_id[s], "- too few Type I probe values.\n\n")
    next
  }

  plot(density(data.m[type1.idx, s]), main = paste("Density Plot for Sample", sample_id[s]))
  d.o <- density(data.m[type2.idx, s])
  points(d.o$x, d.o$y, type="l", col="red")
  print(s)
}
dev.off()


for(s in 1:ncol(data.m)){
  check_for_fail <- grepl("FAILED", sample_id[s])

  if (check_for_fail) {
    next
  }

  beta.v <- data.m[,s];
  bmiq.o <- BMIQ(beta.v,design.v,sampleID=sample_id[s]);
  tmp.v <- bmiq.o$nbeta;
  save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
  print(paste("Done BMIQ for sample ",s,sep=""));
}

bmiq.m <- data.m;

for(s in 1:ncol(bmiq.m)){
  check_for_fail <- grepl("FAILED", sample_id[s])

  if (check_for_fail) {
    next
  }
  
  load(paste("bmiq",s,".Rd",sep=""));
  bmiq.m[,s] <- tmp.v;
  print(s);
}

avg_cols <- grep("AVG", sample_id, value = TRUE)
if (length(avg_cols) == 2) {
  cat("📊 Generating final Type2-BMIQ plot for:", avg_cols[1], "and", avg_cols[2], "\n")

  idx1 <- which(sample_id == avg_cols[1])  # KO AVG
  idx2 <- which(sample_id == avg_cols[2])  # Control AVG

  # Type2 BMIQ-normalized only
  type2_ko_bmiq <- bmiq.m[design.v == 2, idx1]
  type2_ctrl_bmiq <- bmiq.m[design.v == 2, idx2]

  # Densities
  d_type2_ko <- density(type2_ko_bmiq)
  d_type2_ctrl <- density(type2_ctrl_bmiq)

  ymax <- max(d_type2_ko$y, d_type2_ctrl$y)

  pdf("Compare_KO_Control_AVG_Type2Only.pdf", width=6, height=4)
  plot(d_type2_ko, type="l", col="red", lwd=2, xlim=c(0, 1), ylim=c(0, ymax),
      xlab="Beta Value", ylab="Density", main="KO vs Control AVG: BMIQ")  # blank title

  lines(d_type2_ctrl, col="green3", lwd=2)

  # # Custom multicolor title
  # mtext("Type2 BMIQ Profiles", side=3, line=1.2, font=1, col="black", cex=1)
  # mtext("KO ", side=3, line=2.2, adj=0.05, font=2, col="red", cex=1.2)
  # mtext("vs", side=3, line=2.2, adj=0.17, font=1, col="black", cex=1.2)
  # mtext(" Control", side=3, line=2.2, adj=0.27, font=2, col="green3", cex=1.2)

  # In-plot color labels (optional)
  # text(x=0.7, y=ymax*0.95, labels="KO", col="red", cex=1.2, font=2)
  # text(x=0.7, y=ymax*0.88, labels="Control", col="green3", cex=1.2, font=2)

  dev.off()
}


save(bmiq.m,file="bmiq.Rd");
bmiqdf <- as.data.frame(bmiq.m)
write.xlsx(bmiqdf, "bmiq.xlsx")

source_directory <- "./"
results_directory <- "./R Scripts/results/"

file_types <- "\\.pdf$|\\.Rd$"
list_of_files <- list.files(source_directory, pattern=file_types, full.names=TRUE)

# Move each file to the target directory
for (file in list_of_files) {
  file_name <- basename(file)  # Extract the file name
  target_path <- file.path(results_directory, file_name)  # Create target path
  file.rename(file, target_path)  # Move the file
}

# Confirmation message
cat("Moved", length(list_of_files), "files to", results_directory, "\n")

# ### select CpGs (remove CHs) - Wasn't able to do yet...
# annoindex <- which(anno450k.m[,1] %in% rownames(bmiq.m))
# annoindex <- annoindex[match(rownames(bmiq.m),anno450k.m[annoindex])]
# anno.m <- anno450k.m[annoindex,]

# cg.idx <- grep("cg",anno.m[,1]);
# anno2.m <- anno.m[cg.idx,];
# bmiq2.m <- bmiq.m[cg.idx,];
# rm(bmiq.m);
# save(bmiq2.m,file="BMIQ_RESULT.Rd")
# save(anno2.m,file="anno2.m.Rd")


}


# library(readxl)

# probes <- read_excel("./data/probesample.xlsx")
# beta_data <- read_excel("./data/beta.xlsx")

# # Creating data matrix
# probe_ids <- sub("_.*", "", beta_data$`probe set`)
# data.m <- as.matrix(beta_data[, -1])
# sample_id <- colnames(beta_data)[-1]
# rownames(data.m) <- probe_ids
# # save(data.m, file = "data.m.RData")

# # Remove rows with any NA values
# complete_rows <- complete.cases(data.m)
# data.m <- data.m[complete_rows, ]

# # Update probe_ids to match the filtered data.m
# probe_ids <- probe_ids[complete_rows]
# rownames(data.m) <- probe_ids

# print(paste("Removed", sum(!complete_rows), "rows containing NA values"))
# print(paste("Final dimensions of data.m:", nrow(data.m), "rows by", ncol(data.m), "columns"))

# ### DoBMIQ.R
# source("./BMIQ_1.4.R")

# print("Read Success.")


# index <- which(probes$name %in% rownames(data.m))
# index <- index[match(rownames(data.m),probes$name[index])]

# print("Index created successfully.")

# probe_targetid <- probes$targetid[index]
# type_ids <- substr(probe_targetid, nchar(probes$targetid) - 1, nchar(probes$targetid) - 1)

# ###
# type1.idx <- which(type_ids==1);
# type2.idx <- which(type_ids==2);
# design.v <- type_ids


# pdf("Profiles.pdf", width=4, height=3)
# for(s in 1:ncol(data.m)){
#   check_for_fail <- grepl("FAILED", sample_id[s])

#   if (check_for_fail) {
#     next
#   }

#   plot(density(data.m[type1.idx, s]), main = paste("Density Plot for Sample", sample_id[s]))
#   d.o <- density(data.m[type2.idx, s])
#   points(d.o$x, d.o$y, type="l", col="red")
#   print(s)
# }
# dev.off()


# for(s in 1:ncol(data.m)){
#   check_for_fail <- grepl("FAILED", sample_id[s])

#   if (check_for_fail) {
#     next
#   }

#   beta.v <- data.m[,s];
#   bmiq.o <- BMIQ(beta.v,design.v,sampleID=sample_id[s]);
#   tmp.v <- bmiq.o$nbeta;
#   save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
#   print(paste("Done BMIQ for sample ",s,sep=""));
# }

# bmiq.m <- data.m;
# rm(data.m);
# for(s in 1:ncol(bmiq.m)){
#   check_for_fail <- grepl("FAILED", sample_id[s])

#   if (check_for_fail) {
#     next
#   }
  
#   load(paste("bmiq",s,".Rd",sep=""));
#   bmiq.m[,s] <- tmp.v;
#   print(s);
# }


# save(bmiq.m,file="bmiq.Rd");

# source_directory <- "./"
# results_directory <- "./results/"

# file_types <- "\\.pdf$|\\.Rd$"
# list_of_files <- list.files(source_directory, pattern=file_types, full.names=TRUE)

# # Move each file to the target directory
# for (file in list_of_files) {
#   file_name <- basename(file)  # Extract the file name
#   target_path <- file.path(results_directory, file_name)  # Create target path
#   file.rename(file, target_path)  # Move the file
# }

# # Confirmation message
# cat("Moved", length(list_of_files), "files to", results_directory, "\n")


# # ### select CpGs (remove CHs) - Wasn't able to do yet...
# # annoindex <- which(anno450k.m[,1] %in% rownames(bmiq.m))
# # annoindex <- annoindex[match(rownames(bmiq.m),anno450k.m[annoindex])]
# # anno.m <- anno450k.m[annoindex,]

# # cg.idx <- grep("cg",anno.m[,1]);
# # anno2.m <- anno.m[cg.idx,];
# # bmiq2.m <- bmiq.m[cg.idx,];
# # rm(bmiq.m);
# # save(bmiq2.m,file="BMIQ_RESULT.Rd")
# # save(anno2.m,file="anno2.m.Rd")





