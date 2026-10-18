# Alexandre GERMAIN

import os
import csv
import matplotlib.pyplot as plt

# Objectif
# fonctions qui permettent de générer des dotplots

def load_protein_positions(gff_file):
    protein_positions={}                                # dictionnaire qui contiendra : key=protein_id : value=gene_position_in_genome

    with open(gff_file, "r", encoding="utf-8") as f : 
        for line in f :
            if line.startswith("#"):
                continue
            
            parts = line.strip().split("\t")

            if len(parts)<9:
                continue
            
            seqid = parts[0]
            feature_type = parts[2]
            start = int(parts[3])           # début du gène
            end = int(parts[4])             # fin du gene
            attributes = parts[8]

            if feature_type != "CDS":       # on garde seulement les CDS
                continue
            
            protein_id = None

            for attribute in attributes.split(";"):
                if attribute.startswith("protein_id="):
                    protein_id=attribute.replace("protein_id=","")
                    break
            if protein_id is None : 
                continue
            
            position = (start+end)//2                   # position du gène = milieu entre start et end (pour placer le point sur le dotplot)
            protein_positions[protein_id] = position

    return protein_positions

def create_dotplot (blast_result_file, gff_gen1, gff_gen2, output_file = None):

    if output_file is None :
        output_file=os.path.join(os.path.dirname(blast_result_file), "dotplot.png")

    data = get_dotplot_data(blast_result_file, gff_gen1, gff_gen2)

    if len(data)==0:
        print("No BLAST hit found, cannot draw dotplot.")
        return None
    
    x_position = [d["x"] for d in data] #positions sur le Génome1
    y_position = [d["y"] for d in data] #positions sur le Génome2

    plt.figure(figsize=(8,8))
    plt.scatter(x_position, y_position, s=2)

    plt.xlabel("Position on genome 1")
    plt.ylabel("Position on genome 2")
    plt.title("BLAST dotplot")

    plt.savefig(output_file, dpi=300)
    plt.close()
    print("Dotplot image created:", output_file)
    
    return output_file # image.png du dotplot

def get_dotplot_data(blast_result_file, gff_gen1, gff_gen2, cog_gen1=None, cog_gen2=None):
    positions_gen1 = load_protein_positions(gff_gen1)       # chargement des positions des genes dans les génomes 
    positions_gen2 = load_protein_positions(gff_gen2)

    data =[]

    with open(blast_result_file, "r", encoding="utf-8") as f :  #lecture du Blast
        reader = csv.reader(f, delimiter="\t")

        for row in reader : 
            if len(row)< 12:
                continue
            
            qseqid = row[0]         # protéine du génome 1
            sseqid = row[1]         # protéine du génome 2
            pident = row[2]         # % identité
            length = row[3]         # longueur alignée
            evalue = row[10]        # e-value, significativité
            bitscore = row[11]      # qualité de l'alignement

            if qseqid not in positions_gen1:
                continue
            if sseqid not in positions_gen2:
                continue
            
            cog1 = cog_gen1.get(qseqid,{}) if cog_gen1 else {}  # récupération COG
            cog2 = cog_gen2.get(sseqid,{}) if cog_gen2 else {}

            category1 = cog1.get("category","NA")
            category2 = cog2.get("category","NA")

            if category1 == category2 and category1!="NA":
                cog_category = category1
            else:
                cog_category = "Mismatch"


            data.append({                           # chaque hit blast devient un dictionnaire
                "x": positions_gen1[qseqid],
                "y": positions_gen2[sseqid],
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": pident,
                "length": length,
                "evalue": evalue,
                "bitscore": bitscore,
                "cog_category": cog_category,
                "cog1": cog1.get("cog_id", "NA"),
                "cog2": cog2.get("cog_id", "NA"),
                "function1": cog1.get("function_name","NA"),
                "function2": cog2.get("function_name","NA")

            })
    return data             # contient les dictionnaires de tous les points du dotplot

def detect_diagonal_blocks(data, bin_size = 50000, min_points=30):      # détection des régions de synténie
    
    diagonals ={}

    for i,d in enumerate(data):
        x=d["x"]
        y=d["y"]

        diagonal_value = y-x            #calcul des diagonales = du décalage entre les génomes
        diagonal_bin = round(diagonal_value/bin_size)*bin_size      # comme les diagonales ne sont jamais parfaitement alignées, on les regroupe par bin

        if diagonal_bin not in diagonals:
            diagonals[diagonal_bin]=[]
        diagonals[diagonal_bin].append(i)
    detected_blocks = []

    for diagonal_bin, indices in diagonals.items():
        if len(indices)>= min_points:                       # les vraies régions de synténie contiennent beaucoup de genes donc les "petits" bins ne sont pas gardés
            x_values = [data[i]["x"] for i in indices]
            y_values = [data[i]["y"] for i in indices]

            block = {
                "diagonal_bin": diagonal_bin,
                "indices": indices,
                "nb_points": len(indices),  #taille de la région
                "x_min": min(x_values),     #x_min/x_max = région génome 1
                "x_max": max(x_values),     
                "y_min": min(y_values),     #y_min/y_max = région génome 2
                "y_max": max(y_values)}
            detected_blocks.append(block)
    detected_blocks.sort(key=lambda b: b["nb_points"], reverse = True) # tri pour avoir les diagonales les plus importantes d'abord

    return detected_blocks