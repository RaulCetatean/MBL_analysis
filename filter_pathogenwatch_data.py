import polars as pl
import os
import shutil

"""
From the isolates downloaded from Pathogenwatch, I want to keep only the NDM+OXA-48 isolates
Both metadata and Kleborate data were downloaded
"""

df = pl.read_csv("/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/kleborate-pathogenwatch.csv", infer_schema_length=10000)

filtered_df = df.filter(pl.col("Bla_Carb_acquired").str.contains(r"NDM-1|NDM-5")) 
filtered_df = filtered_df.select(["Genome Name","contig_count","N50","total_size","Bla_acquired","Bla_ESBL_acquired","Bla_Carb_acquired"])

# Add metadata and create a single file
metadata_df = pl.read_csv("/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/metadata-Pathogenwatch_all.csv")
second_df = metadata_df.select(["Genome Name","Host","year","month","Country","Isolation source"])

# Combine them based on genome name
combined_df = filtered_df.join(second_df, on="Genome Name", how="left")

# Remove entries with missing Country data, and environmental samples
combined_df = combined_df.filter(pl.col("Country").is_not_null() & (pl.col("Country") != ""))
combined_df = combined_df.filter((pl.col("Host") != "Environmental"))

source_folder = "/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/Kp_Pathogenwatch"
destination_folder = "/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/KP_PW_genomes_v2"

os.makedirs(destination_folder, exist_ok=True)
genome_names = combined_df["Genome Name"].to_list()

not_found = []

for genome_name in genome_names:
    filename = f"{genome_name}.fasta"

    source_path = os.path.join(source_folder, filename)

    if os.path.exists(source_path):
        shutil.copy2(source_path, destination_folder)
    else:
        not_found.append(genome_name)

if not_found:
    print(f"Warning: {len(not_found)} genomes not found")
    print(not_found)
else:
    print("All genomes copied successfully")

combined_df.write_csv("ALL_METADATA_KP.csv")
filtered_df.write_csv("filtered_pathogenwatch_data.csv")

#---------------------------------------------------------------------
# Prepare annotation datasets for iTOL
import polars as pl
import os
import shutil

df = pl.read_csv("ALL_METADATA_KP.csv")
filtered_df = df.filter(pl.col("Bla_Carb_acquired").str.contains("OXA-48"))
unique_countries = filtered_df.select("Country").unique().to_series().to_list()

# NDM
df = filtered_df.with_columns(
    pl.when(pl.col("Bla_Carb_acquired").str.contains("NDM-1"))
    .then(pl.lit("NDM-1"))
    .when(pl.col("Bla_Carb_acquired").str.contains("NDM-5"))
    .then(pl.lit("NDM-5"))
    .otherwise(pl.lit(None))
    .alias("NDM_type")
)

df = df.with_columns(
    pl.when(pl.col("NDM_type") == "NDM-1").then(pl.lit("#806DAB")) 
    .when(pl.col("NDM_type") == "NDM-5").then(pl.lit("#749E9D")) 
    .otherwise(pl.lit(None))
    .alias("NDM_color")
)

# OXA-1
df = df.with_columns(
    pl.when(pl.col("Bla_acquired").str.contains(r"\bOXA-1\b"))
    .then(pl.lit("#C7A6BA"))
    .otherwise(pl.lit(None))
    .alias("OXA1_color")
)


# OXA-9
df = df.with_columns(
    pl.when(pl.col("Bla_acquired").str.contains(r"\bOXA-9\b"))
    .then(pl.lit("#4D2F30"))
    .otherwise(pl.lit(None))
    .alias("OXA9_color")
)

df.write_csv("KP_NDM_OXA48_TEST.csv")

# COUNTRY
country_colors = {
    "Germany": "#1f77b4",
    "Israel": "#ff7f0e",
    "UAE": "#2ca02c",
    "Saudi Arabia": "#d62728",
    "Russia": "#9467bd",
    "Belgium": "#8c564b",
    "Netherlands": "#e377c2",
    "Italy": "#7f7f7f",
    "Denmark": "#bcbd22",
    "Spain": "#17becf",
    "Iran": "#aec7e8",
    "USA": "#ffbb78",
    "United Kingdom": "#98df8a",
    "France": "#ff9896",
    "Serbia": "#c5b0d5",
    "Tunisia": "#c49c94",
    "Egypt": "#f7b6d2",
}

df = df.with_columns(
    pl.col("Country").replace(country_colors).alias("Country_color")
)

# EXTRACT DATASETS FOR ITOL ANNOTATIONS
itol_country = df.select([
    pl.col("Genome Name"),
    pl.col("Country_color")
]).drop_nulls()

itol_oxa1 = df.select([
    pl.col("Genome Name"),
    pl.col("OXA1_color")
]).drop_nulls()

itol_oxa9 = df.select([
    pl.col("Genome Name"),
    pl.col("OXA9_color")
]).drop_nulls()

itol_ndm = df.select([
    pl.col("Genome Name"),
    pl.col("NDM_color"),
    pl.col("NDM_type")
]).drop_nulls()

itol_ndm.write_csv("ANNOT_TABLE_NDM.csv", separator="\t")
itol_oxa1.write_csv("ANNOT_TABLE_OXA1.csv", separator="\t")
itol_oxa9.write_csv("ANNOT_TABLE_OXA9.csv", separator="\t")
itol_country.write_csv("ANNOT_TABLE_COUNTRY.csv", separator="\t")

source_folder = "/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/Kp_Pathogenwatch"
destination_folder = "/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/KP_NDM_OXA48"

os.makedirs(destination_folder, exist_ok=True)
genome_names = filtered_df["Genome Name"].to_list()

not_found = []

for genome_name in genome_names:
    filename = f"{genome_name}.fasta"

    source_path = os.path.join(source_folder, filename)

    if os.path.exists(source_path):
        shutil.copy2(source_path, destination_folder)
    else:
        not_found.append(genome_name)

if not_found:
    print(f"Warning: {len(not_found)} genomes not found")
    print(not_found)
else:
    print("All genomes copied successfully")