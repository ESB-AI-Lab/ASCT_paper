# -------------------------------
# Workflow: All accessions per family using E-values (-log10 scale)
# -------------------------------

library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

# Load CSV
df <- read.csv("H_T_dom.csv", header = TRUE, stringsAsFactors = FALSE)

# Fix column name if needed
colnames(df) <- str_replace(colnames(df), "evolutionary_CoA._family", "evolutionary_CoA_family")

# Clean families and convert E-values
df <- df %>%
  mutate(
    traditional_CoA_family = str_trim(traditional_CoA_family),
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    evolutionary_CoA_family = case_when(
      evolutionary_CoA_family %in% c("OXCT1","Cat1") ~ evolutionary_CoA_family,
      TRUE ~ "Unknown"
    ),
    e_value = as.numeric(str_trim(e_value)),
    neglog_evalue = -log10(e_value)
  )

# -------------------------------
# Summarize mean -log10(E-value) per accession per family
# -------------------------------

# Traditional CoA family
accession_summary_traditional <- df %>%
  group_by(traditional_CoA_family, accession) %>%
  summarise(mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
            n_sequences = n(),
            .groups = "drop") %>%
  arrange(traditional_CoA_family, desc(mean_neglog_evalue))

# Evolutionary CoA family
accession_summary_evolutionary <- df %>%
  group_by(evolutionary_CoA_family, accession) %>%
  summarise(mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
            n_sequences = n(),
            .groups = "drop") %>%
  arrange(evolutionary_CoA_family, desc(mean_neglog_evalue))

# -------------------------------
# Visualization: Bar plots for all accessions
# -------------------------------

# Traditional CoA family
svglite("Traditional_CoA_family_all_accessions.svg", width = 12, height = 6)
ggplot(accession_summary_traditional, 
       aes(x = reorder(accession, mean_neglog_evalue), 
           y = mean_neglog_evalue, fill = traditional_CoA_family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~traditional_CoA_family, scales = "free_x") +
  geom_text(aes(label = round(mean_neglog_evalue,1)), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Accession") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Accessions in Traditional CoA Families (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Evolutionary CoA family
svglite("Evolutionary_CoA_family_all_accessions.svg", width = 12, height = 6)
ggplot(accession_summary_evolutionary, 
       aes(x = reorder(accession, mean_neglog_evalue), 
           y = mean_neglog_evalue, fill = evolutionary_CoA_family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~evolutionary_CoA_family, scales = "free_x") +
  geom_text(aes(label = round(mean_neglog_evalue,1)), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Accession") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Accessions in Evolutionary CoA Families (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ----------------------
# Chart values for all accessions per family (used for bar charts)lues
# ----------------------

library(dplyr)
library(stringr)

# Load CSV
df <- read.csv("H_T_dom.csv", header = TRUE, stringsAsFactors = FALSE)

# Fix column name if needed
colnames(df) <- str_replace(colnames(df), "evolutionary_CoA._family", "evolutionary_CoA_family")

# Clean families and convert E-values
df <- df %>%
  mutate(
    traditional_CoA_family = str_trim(traditional_CoA_family),
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    evolutionary_CoA_family = case_when(
      evolutionary_CoA_family %in% c("OXCT1","Cat1") ~ evolutionary_CoA_family,
      TRUE ~ "Unknown"
    ),
    e_value = as.numeric(str_trim(e_value)),
    neglog_evalue = -log10(e_value)
  )

# -------------------------------
# Mean -log10(E-value) per accession per family (used for bar charts)
# -------------------------------

# Traditional CoA family
chart_values_traditional <- df %>%
  group_by(traditional_CoA_family, accession) %>%
  summarise(
    mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
    n_sequences = n(),
    .groups = "drop"
  ) %>%
  arrange(traditional_CoA_family, desc(mean_neglog_evalue))

# Evolutionary CoA family
chart_values_evolutionary <- df %>%
  group_by(evolutionary_CoA_family, accession) %>%
  summarise(
    mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
    n_sequences = n(),
    .groups = "drop"
  ) %>%
  arrange(evolutionary_CoA_family, desc(mean_neglog_evalue))

# -------------------------------
# Save to CSV files
# -------------------------------
write.csv(chart_values_traditional, "Chart_values_Traditional_CoA.csv", row.names = FALSE)
write.csv(chart_values_evolutionary, "Chart_values_Evolutionary_CoA.csv", row.names = FALSE)

# Optional: view
head(chart_values_traditional)
head(chart_values_evolutionary)



#####




# -------------------------------
# Workflow: All domains per family using E-values (-log10 scale)
# -------------------------------

# Load libraries
install.packages(c("dplyr", "stringr", "ggplot2", "svglite"))
library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

# -------------------------------
# 1. Load CSV
# -------------------------------
df <- read.csv("H_T_dom.csv", header = TRUE, stringsAsFactors = FALSE)

# -------------------------------
# 2. Clean family columns and convert E-values
# -------------------------------
# Fix column name if needed
colnames(df) <- str_replace(colnames(df), "evolutionary_CoA._family", "evolutionary_CoA_family")

df <- df %>%
  mutate(
    # Clean family columns
    traditional_CoA_family = str_trim(traditional_CoA_family),
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    evolutionary_CoA_family = case_when(
      evolutionary_CoA_family %in% c("OXCT1","Cat1") ~ evolutionary_CoA_family,
      TRUE ~ "Unknown"
    ),
    # Convert E-value to numeric
    e_value = as.numeric(str_trim(e_value)),
    # Create -log10(E-value) for plotting
    neglog_evalue = -log10(e_value)
  )

# -------------------------------
# 3. Summarize mean -log10(E-value) per domain per family
# -------------------------------

# Traditional CoA family
domain_summary_traditional <- df %>%
  group_by(traditional_CoA_family, short_name) %>%
  summarise(mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
            n_sequences = n(),
            .groups = "drop") %>%
  arrange(traditional_CoA_family, desc(mean_neglog_evalue))

# Evolutionary CoA family
domain_summary_evolutionary <- df %>%
  group_by(evolutionary_CoA_family, short_name) %>%
  summarise(mean_neglog_evalue = mean(neglog_evalue, na.rm = TRUE),
            n_sequences = n(),
            .groups = "drop") %>%
  arrange(evolutionary_CoA_family, desc(mean_neglog_evalue))

# -------------------------------
# 4. Visualization: Bar plots for all domains
# -------------------------------

# Traditional CoA family
svglite("Traditional_CoA_family_all_domains.svg", width = 10, height = 6)
ggplot(domain_summary_traditional, aes(x = reorder(short_name, mean_neglog_evalue), 
                                       y = mean_neglog_evalue, fill = traditional_CoA_family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~traditional_CoA_family, scales = "free_x") +
  geom_text(aes(label = round(mean_neglog_evalue,1)), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Domain") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Domains in Traditional CoA Families (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Evolutionary CoA family
svglite("Evolutionary_CoA_family_all_domains.svg", width = 10, height = 6)
ggplot(domain_summary_evolutionary, aes(x = reorder(short_name, mean_neglog_evalue), 
                                        y = mean_neglog_evalue, fill = evolutionary_CoA_family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~evolutionary_CoA_family, scales = "free_x") +
  geom_text(aes(label = round(mean_neglog_evalue,1)), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Domain") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Domains in Evolutionary CoA Families (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()







########



# -------------------------------
# Top domains by E-value only
# -------------------------------

# Load libraries
install.packages(c("dplyr", "stringr", "ggplot2", "svglite"))
library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

# -------------------------------
# 1. Load CSV
# -------------------------------
df <- read.csv("H_T_dom.csv", header = TRUE, stringsAsFactors = FALSE)

# Show column names with quotes to spot extra spaces
cat(paste0('"', colnames(df), '"'), sep = "\n")

# Remove leading/trailing whitespace in column names
colnames(df) <- str_trim(colnames(df))

# Now check
colnames(df)

# Rename the column to remove the extra period
colnames(df) <- str_replace(colnames(df), "evolutionary_CoA._family", "evolutionary_CoA_family")

# Verify
colnames(df)
# -------------------------------
# 2. Clean family columns & convert e_value
# -------------------------------
###df <- df %>%
  ##mutate(
    # Trim whitespace
    traditional_CoA_family = str_trim(traditional_CoA_family),
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    
    # Standardize traditional families
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    
    # Standardize evolutionary families
    evolutionary_CoA._family = case_when(
      evolutionary_CoA._family %in% c("OXCT1","Cat1") ~ evolutionary_CoA._family,
      TRUE ~ "Unknown"
    ),
    
    # Convert e_value to numeric (handles scientific notation)
    e_value = as.numeric(str_trim(e_value))
  )
###

library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

df <- df %>%
  mutate(
    # Clean family columns
    traditional_CoA_family = str_trim(traditional_CoA_family),
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    evolutionary_CoA_family = case_when(
      evolutionary_CoA_family %in% c("OXCT1","Cat1") ~ evolutionary_CoA_family,
      TRUE ~ "Unknown"
    ),
    # Convert E-value to numeric
    e_value = as.numeric(str_trim(e_value))
  )

# Now you can safely run the E-value analysis workflow for both families
# -------------------------------
# 3. Top domain per family by MEAN E-VALUE
# -------------------------------

# Traditional CoA family
top_domain_evalue_traditional <- df %>%
  group_by(traditional_CoA_family, short_name) %>%
  summarise(mean_evalue = mean(e_value, na.rm = TRUE), .groups = "drop") %>%
  arrange(traditional_CoA_family, mean_evalue) %>%
  group_by(traditional_CoA_family) %>%
  slice_min(mean_evalue, n = 1)  # pick domain with lowest avg E-value

# Evolutionary CoA family
top_domain_evalue_evolutionary <- df %>%
  group_by(evolutionary_CoA_family, short_name) %>%
  summarise(mean_evalue = mean(e_value, na.rm = TRUE), .groups = "drop") %>%
  arrange(evolutionary_CoA_family, mean_evalue) %>%
  group_by(evolutionary_CoA_family) %>%
  slice_min(mean_evalue, n = 1)  # pick domain with lowest avg E-value

# -------------------------------
# 4. Visualization: Bar plots (Mean E-value)
# -------------------------------

# Traditional CoA family
ggplot(top_domain_evalue_traditional, 
       aes(x = reorder(traditional_CoA_family, mean_evalue), 
           y = mean_evalue, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  scale_y_log10() +   # <-- log scale for E-values
  theme_minimal() +
  xlab("Traditional CoA Family") +
  ylab("Mean E-value (log10 scale)") +
  ggtitle("Most Significant Domain per Traditional CoA Family (Mean E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

svglite("Traditional_CoA_family_evalue.svg", width = 8, height = 6)
ggplot(top_domain_evalue_traditional, aes(x = reorder(traditional_CoA_family, mean_evalue), 
                                          y = mean_evalue, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Traditional CoA Family") +
  ylab("Mean E-value") +
  ggtitle("Most Significant Domain per Traditional CoA Family (Mean E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()

# Evolutionary CoA family
svglite("Evolutionary_CoA_family_evalue.svg", width = 8, height = 6)
ggplot(top_domain_evalue_evolutionary, aes(x = reorder(evolutionary_CoA_family, mean_evalue), 
                                           y = mean_evalue, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Evolutionary CoA Family") +
  ylab("Mean E-value") +
  ggtitle("Most Significant Domain per Evolutionary CoA Family (Mean E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()

library(ggplot2)
library(svglite)

# Traditional CoA family plot
svglite("Traditional_CoA_family_evalue_neglog.svg", width = 8, height = 6)
ggplot(top_domain_evalue_traditional, 
       aes(x = reorder(traditional_CoA_family, -mean_evalue), 
           y = -log10(mean_evalue), fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Traditional CoA Family") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Most Significant Domain per Traditional CoA Family (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()

# Evolutionary CoA family plot
svglite("Evolutionary_CoA_family_evalue_neglog.svg", width = 8, height = 6)
ggplot(top_domain_evalue_evolutionary, 
       aes(x = reorder(evolutionary_CoA_family, -mean_evalue), 
           y = -log10(mean_evalue), fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Evolutionary CoA Family") +
  ylab("-log10(Mean E-value)") +
  ggtitle("Most Significant Domain per Evolutionary CoA Family (-log10 E-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()








###

# -------------------------------
# Load required libraries
# -------------------------------
install.packages(c("dplyr", "stringr", "ggplot2"))
library(dplyr)
library(stringr)
library(ggplot2)

# -------------------------------
# 1. Load your CSV
# -------------------------------
df <- read.csv("H_T_domains.csv", header = TRUE, stringsAsFactors = FALSE)

# -------------------------------
# 2. Clean family columns
# -------------------------------
df <- df %>%
  mutate(
    # Clean traditional CoA family
    traditional_CoA_family = str_trim(traditional_CoA_family),
    traditional_CoA_family = case_when(
      traditional_CoA_family %in% c("IA","IB","IC") ~ traditional_CoA_family,
      TRUE ~ "Unknown"
    ),
    # Clean evolutionary CoA family
    evolutionary_CoA_family = str_trim(evolutionary_CoA_family),
    evolutionary_CoA_family = case_when(
      evolutionary_CoA_family %in% c("OXCT1","Cat1") ~ evolutionary_CoA_family,
      TRUE ~ "Unknown"
    )
  )

# -------------------------------
# 3. Compute most common domain per family
# -------------------------------

# Traditional CoA family
most_common_traditional <- df %>%
  group_by(traditional_CoA_family, short_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(traditional_CoA_family, desc(count)) %>%
  group_by(traditional_CoA_family) %>%
  slice_max(count, n = 1)  # top domain per family

# Evolutionary CoA family
most_common_evolutionary <- df %>%
  group_by(evolutionary_CoA_family, short_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(evolutionary_CoA_family, desc(count)) %>%
  group_by(evolutionary_CoA_family) %>%
  slice_max(count, n = 1)  # top domain per family

# -------------------------------
# 4. Visualization: Bar plots
# -------------------------------

# Traditional CoA family
ggplot(most_common_traditional, aes(x = reorder(traditional_CoA_family, -count), 
                                    y = count, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Traditional CoA Family") +
  ylab("Number of Sequences") +
  ggtitle("Most Common Domain per Traditional CoA Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# Evolutionary CoA family
ggplot(most_common_evolutionary, aes(x = reorder(evolutionary_CoA_family, -count), 
                                     y = count, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Evolutionary CoA Family") +
  ylab("Number of Sequences") +
  ggtitle("Most Common Domain per Evolutionary CoA Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# svg

install.packages("svglite")
library(svglite)

# Save Traditional CoA family plot
svglite("Traditional_CoA_family_plot.svg", width = 8, height = 6)
ggplot(most_common_traditional, aes(x = reorder(traditional_CoA_family, -count), 
                                    y = count, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Traditional CoA Family") +
  ylab("Number of Sequences") +
  ggtitle("Most Common Domain per Traditional CoA Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()

# Save Evolutionary CoA family plot
svglite("Evolutionary_CoA_family_plot.svg", width = 8, height = 6)
ggplot(most_common_evolutionary, aes(x = reorder(evolutionary_CoA_family, -count), 
                                     y = count, fill = short_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = short_name), vjust = -0.5, size = 3) +
  theme_minimal() +
  xlab("Evolutionary CoA Family") +
  ylab("Number of Sequences") +
  ggtitle("Most Common Domain per Evolutionary CoA Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
dev.off()
