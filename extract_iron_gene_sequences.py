import os
from pathlib import Path
from Bio import SeqIO


"""
Extract gene sequences from the .gbff files from each isolates. Annotation was performed with Bakta
"""

cwd = os.getcwd()
isolate_folders = [f for f in os.listdir(cwd) if "." not in f] #bakta_dir

 gene_list = ["fhuA", "fepA", "fbpA", "efeO","exbB","fiu","baeS","two-component system sensor histidine kinase EnvZ","cirA","feoA","sitC","apbC","fepG","fepC","fetB","fetA","fecA",
               "fiu","ompR","yicI","yicJ","yicL","chrA","tonB","sugE","yicM","PmrB","pmrB","ompK37", "exbD","two-component system response regulator OmpR", "Chromate transport protein", "Outer Membrane Siderophore Receptor IroN", "pbpC"]

gene_list = []

# Deduplicate while preserving order
gene_list = list(dict.fromkeys(gene_list))

faa_files = {
    isolate: f"{cwd}/{isolate}/{isolate}.faa"
    for isolate in isolate_folders
}

output_dir = f"/run/media/raulc/Dati/Sequencing_runs_analysis/MBL_analyses_all/out/extracted_genes" #hardcoded for now
os.makedirs(output_dir, exist_ok=True)

priority_isolate = "KP1PI" # Reference genome on top

for gene in gene_list:
    output_file = f"{output_dir}/{gene}_sequences.fasta"
    gene_lower = gene.lower()

    # Reorder so priority isolate is first
    ordered_isolates = sorted(
        faa_files.items(),
        key=lambda x: (0 if x[0] == priority_isolate else 1)
    )

    with open(output_file, "w") as out_f:
        for isolate, faa_file in ordered_isolates:
            found = False

            if not os.path.exists(faa_file):
                print(f"[Error] FAA file missing for {isolate}: {faa_file}")
                continue

            for record in SeqIO.parse(faa_file, "fasta"):
                if gene_lower in record.description.lower():
                    record.id = isolate
                    record.description = ""
                    SeqIO.write(record, out_f, "fasta")
                    found = True

            if not found:
                print(f"[Warning] '{gene}' not found in {isolate}")

print("Done. Files written to: ", output_dir)