# Alexandre GERMAIN

import os
import subprocess
import csv
import re

import tools

# Objectif : 
# Associer chaque protéine du génome à une catégorie fonctionnelle COG grace à RPSBLAST
# Créer un dictionnaire: protéine->fonction biologique

def load_cog_names(cognames_file):          # Lecture du fichier cogname2003-2014.tab
    cog_info ={}

    with open(cognames_file, "r", encoding = "utf-8") as f:
        reader = csv.reader(f,delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            
            cog_id = row[0]
            category = row[1]
            description =row[2]

            cog_info[cog_id]={
                "category": category,
                "description" : description
            }
    return cog_info

def load_function_names(fun_file):            # Lecture du fichier fun2003-2014.tab
    function_names ={}

    with open(fun_file, "r", encoding = "utf-8") as f:
        reader = csv.reader(f,delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            
            category = row[0]
            name =row[1]

            function_names[category]= name
    return function_names

def run_rpsblast(protein_fasta, genome_name, cog_db_path, output_dir="cog_results", evalue ="1e-10"):       #Lancer rpsblast contre la base COG
    os.makedirs(output_dir, exist_ok = True)        #dossier résultat

    safe_name = tools.safe_filename(genome_name)
    output_file = os.path.join(output_dir, f"{safe_name}_cog.tsv")

    # Fichier vide = analyse interrompue, on relance.
    if tools.is_valid_result(output_file):  # on évite de relancer rpsblast inutilement
        print("Existing COG result reused:", output_file)
        return output_file

    # Debian et Ubuntu installent le binaire sous le nom "rpsblast+" ; partout
    # ailleurs (Windows, macOS, conda, bioconda) il s'appelle "rpsblast".
    rpsblast_bin = tools.find_blast_tool("rpsblast")

    cmd = [
        rpsblast_bin,
        "-query", protein_fasta,
        "-db", cog_db_path,
        "-num_alignments","1", # un seul meilleur hit COG par protéine
        "-evalue", str(evalue),
        "-outfmt", "6 qaccver stitle evalue",
        "-out", output_file
    ]

    print("Running COG analysis...")
    subprocess.run(cmd, check=True)             # exécution

    print("COG analysis finished:", output_file)
    return output_file

def parse_rpsblast_results(rpsblast_file, cognames_file, fun_file):     # on transforme les résultats rpsblast bruts en un dictionnaire exploitable
    cog_info = load_cog_names(cognames_file)
    function_names =load_function_names(fun_file)

    protein_to_cog ={}

    with open(rpsblast_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter ="\t")

        for row in reader:
            if len(row)< 3:
                continue
            
            protein_id = row[0]
            stitle = row[1]
            evalue = row[2]

            match = re.search(r'(COG\d+)', stitle)
            if match is None :
                continue
            cog_id = match.group(1)            

            #cog_id = stitle.split(",")[0].strip()

            if cog_id not in cog_info:
                continue
            
            category = cog_info[cog_id]["category"][0]
            description = cog_info[cog_id]["description"]

            protein_to_cog[protein_id]={
                "cog_id": cog_id,
                "category": category,
                "description": description,
                "function_name": function_names.get(category, "Fonction Inconnue"),
                "evalue": evalue
            }
    return protein_to_cog