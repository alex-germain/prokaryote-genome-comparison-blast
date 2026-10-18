# Alexandre GERMAIN

import csv
import os

import tools

# Objectif : Création d'une DB locale qui regroupe les métadonnées utiles dans "prokaryotes_db.tsv" de la forme : 
# assembly_accession	taxid	organism_name	    refseq_category	    assembly_level	ftp_path	domain
# GCF_000005845.2	    562	    Escherichia coli	reference genome	Complete Genome	https://...	Bacteria

# on utilise tsv parce que les colonnes sont automatiquement séparées par des tabulations

# cette DB locale est utile pour permettre àl'utilisateur de sélectionner le génome qui l'intéresse sans avoir à faire des requètes vers NCBI pour mettre à jour la listbox.
# La DB contient le ftp_path de chaque génome donc la requète pour télécharger le génome est faite uniquement une fois que le génome est sélectionné

url_bact = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
url_archaea = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'Local_DB') # on enregistre les BD dans le fichier "Local_DB"

bacteria_path = os.path.join(data_dir, 'bacteria.txt')
archaea_path = os.path.join(data_dir, 'archaea.txt')

output_file = os.path.join(data_dir, 'procaryote_db.tsv')

# Nom public du chemin de la base : l'interface le lit ici plutot que de
# reecrire le chemin de son cote, ou les deux finissaient par diverger.
DB_PATH = output_file

# Téléchargement des listes des métadonnées NCBI pour les bactéries et les archées
def download_metadata(url, output_path):
    # Ces deux fichiers font plusieurs dizaines de Mo. Un telechargement
    # coupe laissait auparavant un fichier tronque, que create_DB() voyait
    # ensuite comme deja present : l'erreur n'apparaissait qu'au moment de
    # l'analyse, sous la forme d'un en-tete introuvable.
    tools.download_file(url, output_path, timeout=120)
    print("Telecharge ici : ", output_path)

# Création de la DB "prokaryotes_db.tsv" à partir des métadonnées téléchargées
def extract_data(input_path, domain):
    records = []
    with open(input_path, "r", encoding="utf-8") as f : 
        reader = csv.reader(f, delimiter="\t")
        header = None
        
        for row in reader : 
            if not row : 
                continue
            if row[0].startswith("#assembly_accession"):
                header = [col.lstrip('#').strip() for col in row]
                break
        if header is None : 
            raise ValueError(f'En-tête introuvable dans {input_path}')
        
        col_index = {col_name : i for i, col_name in enumerate(header)}
        
        for i,row in enumerate(reader) : 
            if i%1000==0:
                print(f"{i} génomes traités...")

            if not row or len(row)<len(header):
                continue
            
            record = {
                "assembly_accession": row[col_index["assembly_accession"]],
                "taxid": row[col_index["taxid"]],
                "organism_name": row[col_index["organism_name"]],
                "refseq_category": row[col_index["refseq_category"]],
                "assembly_level": row[col_index["assembly_level"]],
                "ftp_path": row[col_index["ftp_path"]],
                "domain": domain
                }
            records.append(record)
    return records

def create_DB_file(records, output_path):
    columns= [
    "assembly_accession",
    "taxid",
    "organism_name",
    "refseq_category",
    "assembly_level",
    "ftp_path",
    "domain"]
    with open(output_path, 'w', encoding='utf-8', newline='') as f : 
        writer = csv.DictWriter(f, fieldnames= columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)

def create_DB():
    if not os.path.exists(bacteria_path):
        download_metadata(url_bact, bacteria_path)
    if not os.path.exists(archaea_path):
        download_metadata(url_archaea, archaea_path)
    
    all_records = []
    all_records.extend(extract_data(bacteria_path, "Bacteria"))
    all_records.extend(extract_data(archaea_path, "Archaea"))

    print("Nb total de records : ", len(all_records))
    print("bact path :", bacteria_path)
    print("archea path ;", archaea_path)
    create_DB_file(all_records, output_file)
