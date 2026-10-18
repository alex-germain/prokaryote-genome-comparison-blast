# Alexandre GERMAIN

import subprocess
import os

import tools

def run_blastp(query_fasta, subject_fasta,query_name, subject_name, evalue="1e-20"):

    # Les noms d'organismes NCBI contiennent des caractères interdits dans un
    # nom de fichier selon le système (parenthèses, apostrophes, deux-points).
    # safe_filename ne garde que lettres, chiffres, point, tiret et underscore.
    safe_query=tools.safe_filename(query_name)
    safe_subject=tools.safe_filename(subject_name)
    safe_evalue=str(evalue).replace("-","m").replace("+","p").replace(".","_")

    base_dir = "blast_results"
    os.makedirs(base_dir, exist_ok=True)

    output_dir= os.path.join(base_dir, f"blast_results_{safe_query}_VS_{safe_subject}_evalue_{safe_evalue}")

    os.makedirs(output_dir, exist_ok=True)

    #Nom de la base de BLAST
    db_name_dir = f"blast_db_{safe_query}_VS_{safe_subject}"
    db_name = os.path.join(output_dir, db_name_dir)

    #Fichier Résultat
    result_file_dir = f"resultats_blast_{safe_query}_VS_{safe_subject}_evalue_{safe_evalue}.tsv"
    result_file = os.path.join(output_dir, result_file_dir)

    # Un BLAST interrompu laisse un fichier vide : l'existence seule ne suffit
    # pas à conclure que le résultat est réutilisable.
    if tools.is_valid_result(result_file):
        print("Existing BLAST result reused:", result_file)
        return result_file

    # Les binaires sont résolus ici : si BLAST+ n'est pas installé, l'erreur
    # explique quoi installer au lieu d'un FileNotFoundError nu.
    makeblastdb_bin = tools.find_blast_tool("makeblastdb")
    blastp_bin = tools.find_blast_tool("blastp")

    cmd_makeblastdb = [makeblastdb_bin, "-in", subject_fasta, "-dbtype","prot","-out", db_name]
    print("Building BLAST database...")
    subprocess.run(cmd_makeblastdb, check=True)

    #Lancement du BLAST
    cmd_blastp =[blastp_bin, "-query", query_fasta, "-db", db_name, "-out", result_file, 
                 "-evalue", str(evalue), "-outfmt", 
                 "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore"]
    print("Running BLAST...")
    subprocess.run(cmd_blastp, check=True)

    print("BLAST finished:", result_file)
    return result_file