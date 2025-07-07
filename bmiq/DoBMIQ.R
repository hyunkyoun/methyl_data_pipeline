library(readxl)
library(openxlsx)

DoBMIQ <- function(probes_path, beta_path) {

probes <- read_excel(probes_path)
print(paste("probes_path:", probes_path))
beta_data <- read.csv(beta_path)
print(paste("beta_path:", beta_path))

# Creating data matrix - USE ALL COLUMNS EXCEPT THE FIRST ONE (probe IDs)
# Extract probe IDs from the first column (no underscore removal needed)
probe_ids <- beta_data[, 1]  # First column contains probe IDs like "cg00101675"

# Use all columns except the first one (probe ID column)
all_sample_id <- colnames(beta_data)[-1]
selected_cols <- 2:ncol(beta_data)  # All columns except the first

data.m <- as.matrix(beta_data[, selected_cols])
sample_id <- colnames(beta_data)[selected_cols]
rownames(data.m) <- probe_ids

cat("Original columns:", length(all_sample_id), "\n")
cat("Selected columns (all except TargetID):", length(sample_id), "\n")
cat("Selected column names:", paste(sample_id, collapse=", "), "\n")

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
source("./bmiq/BMIQ_1.4.R")

print("Read Success.")

# --- Match probe_ids to probe annotation ---
index <- match(probe_ids, probes$name)
matched_probes <- probes[index, ]

# Extract targetid
probe_targetid <- matched_probes$targetid

# Extract probe type: last two characters of targetid (e.g. "11", "21")
type_ids <- substr(probe_targetid, nchar(probe_targetid) - 1, nchar(probe_targetid))
type_ids <- as.numeric(type_ids)

# Build design vector and indices - FILTER OUT INVALID TYPES
valid_types <- type_ids %in% c(11, 21) & !is.na(type_ids)

# Apply the filter to all related objects
type_ids <- type_ids[valid_types]
data.m <- data.m[valid_types, ]
probe_ids <- probe_ids[valid_types]
matched_probes <- matched_probes[valid_types, ]

# Create design vector with explicit values
design.v <- rep(NA, length(type_ids))
design.v[type_ids == 11] <- 1
design.v[type_ids == 21] <- 2

# Final safety check - remove any remaining NAs
final_valid <- !is.na(design.v)
if(sum(!final_valid) > 0) {
  cat("Removing", sum(!final_valid), "additional invalid probes\n")
  design.v <- design.v[final_valid]
  data.m <- data.m[final_valid, ]
  probe_ids <- probe_ids[final_valid]
  type_ids <- type_ids[final_valid]
}

type1.idx <- which(type_ids == 11)
type2.idx <- which(type_ids == 21)

# Debug print
cat("Probe type distribution:\n")
print(table(type_ids, useNA = "ifany"))
cat("Total Type 1 probes:", length(type1.idx), "\n")
cat("Total Type 2 probes:", length(type2.idx), "\n")

pdf("Profiles.pdf", width=4, height=3)
for(s in 1:ncol(data.m)){
  cat("Sample:", sample_id[s], "\n")
  cat("Length of type1.idx:", length(type1.idx), "\n")
  cat("Non-NA values in data.m[type1.idx, s]:", sum(!is.na(data.m[type1.idx, s])), "\n")

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
  beta.v <- data.m[,s];
  
  # Add safety checks before calling BMIQ
  cat("Sample", s, ":", sample_id[s], "\n")
  cat("Design vector range:", range(design.v), "\n")
  cat("Design vector unique values:", sort(unique(design.v)), "\n")
  cat("Beta values range:", range(beta.v, na.rm=TRUE), "\n")
  cat("Type 1 probes:", sum(design.v == 1), "\n")
  cat("Type 2 probes:", sum(design.v == 2), "\n")
  
  # Ensure design vector and beta vector have same length
  if(length(design.v) != length(beta.v)) {
    cat("ERROR: Length mismatch - design.v:", length(design.v), "beta.v:", length(beta.v), "\n")
    next
  }
  
  # Check for any values other than 1 and 2 in design vector
  if(any(!design.v %in% c(1, 2))) {
    cat("ERROR: Invalid values in design vector:", unique(design.v), "\n")
    next
  }
  
  # Check if we have enough data points for each type
  type1_data <- beta.v[design.v == 1]
  type2_data <- beta.v[design.v == 2]
  
  if(sum(!is.na(type1_data)) < 50 || sum(!is.na(type2_data)) < 50) {
    cat("⚠️ Skipping sample", sample_id[s], "- insufficient data points\n")
    # Create dummy normalized data (just copy original)
    tmp.v <- beta.v
    save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
    next
  }
  
  bmiq.o <- BMIQ(beta.v,design.v,sampleID=sample_id[s]);
  tmp.v <- bmiq.o$nbeta;
  save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
  print(paste("Done BMIQ for sample ",s,sep=""));
}

bmiq.m <- data.m;

for(s in 1:ncol(bmiq.m)){
  load(paste("bmiq",s,".Rd",sep=""));
  bmiq.m[,s] <- tmp.v;
  print(s);
}

# Generate comparison plot using first two samples
if (ncol(bmiq.m) >= 2) {
  cat("📊 Generating final Type2-BMIQ plot for:", sample_id[1], "and", sample_id[2], "\n")

  # Type2 BMIQ-normalized only - comparing first two samples
  type2_sample1_bmiq <- bmiq.m[design.v == 2, 1]
  type2_sample2_bmiq <- bmiq.m[design.v == 2, 2]

  # Densities
  d_type2_sample1 <- density(type2_sample1_bmiq)
  d_type2_sample2 <- density(type2_sample2_bmiq)

  ymax <- max(d_type2_sample1$y, d_type2_sample2$y)

  pdf("Compare_First_Two_Samples_Type2Only.pdf", width=6, height=4)
  plot(d_type2_sample1, type="l", col="red", lwd=2, xlim=c(0, 1), ylim=c(0, ymax),
      xlab="Beta Value", ylab="Density", main="Sample Comparison: BMIQ (Type 2 Probes)")

  lines(d_type2_sample2, col="green3", lwd=2)
  
  # Add legend
  legend("topright", legend=c(sample_id[1], sample_id[2]), 
         col=c("red", "green3"), lwd=2, cex=0.8)

  dev.off()
} else {
  cat("📊 Less than 2 samples available - skipping comparison plot\n")
}

save(bmiq.m,file="bmiq.Rd");

# Remove X prefix from column names
colnames(bmiq.m) <- gsub("^X", "", colnames(bmiq.m))

bmiqdf <- as.data.frame(bmiq.m)

# Add CpG IDs as the first column
bmiqdf <- data.frame(CpG_ID = rownames(bmiqdf), bmiqdf, stringsAsFactors = FALSE)

write.xlsx(bmiqdf, "bmiq.xlsx")

source_directory <- "./"
results_directory <- "./bmiq/results"

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

}