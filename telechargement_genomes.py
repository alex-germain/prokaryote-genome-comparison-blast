# Alexandre GERMAIN

import os
import gzip
import shutil

import tools

# Objectif :
# Telecharger depuis le NCBI, pour un genome selectionne, le fichier des
# proteines (.faa) et le fichier d'annotation (.gff), et les decompresser.


def download_and_unzip(url, output_path):
    """
    Telecharge `url` vers `output_path` (.gz) puis le decompresse a cote.

    Les deux etapes sont atomiques et sautees si leur resultat est deja
    present : relancer une comparaison deja faite ne retelecharge rien, mais
    un telechargement ou une decompression interrompus sont bien refaits.
    """
    # Verifie avant de telecharger : inutile de rapatrier des dizaines de Mo
    # pour echouer ensuite sur le nom du fichier.
    if not output_path.endswith(".gz"):
        raise ValueError(f"Fichier compresse attendu (.gz) : {output_path}")
    uncompressed_path = output_path[:-len(".gz")]

    if not os.path.exists(output_path):
        tools.download_file(url, output_path)

    if not os.path.exists(uncompressed_path):
        def unzip(part_path):
            with gzip.open(output_path, "rb") as f_in:
                with open(part_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        tools.atomic_write(uncompressed_path, unzip)

    return uncompressed_path


def download_genome(genome, output_dir="genomes"):
    """
    Telecharge les proteines et les annotations d'un genome du NCBI.

    `genome` est une ligne de la base locale, qui contient le ftp_path.
    Renvoie (chemin_proteines, chemin_gff), ou (None, None) si le NCBI ne
    publie pas de ftp_path pour cet assemblage.
    """
    ftp_path = genome["ftp_path"].strip()

    if ftp_path in ("na", ""):
        print("No ftp_path available for this genome")
        return None, None

    os.makedirs(output_dir, exist_ok=True)

    # Le NCBI nomme ses fichiers d'apres le dernier segment du ftp_path :
    # .../GCF_000005845.2_ASM584v2 donne GCF_000005845.2_ASM584v2_protein.faa.gz
    ftp_path = ftp_path.rstrip("/")
    assembly_name = ftp_path.split("/")[-1]

    protein_filename = assembly_name + "_protein.faa.gz"
    protein_url = ftp_path + "/" + protein_filename
    protein_output_path = os.path.join(output_dir, protein_filename)

    # Le GFF donne la position des genes sur le genome : c'est lui qui permet
    # de placer les points du dotplot selon leurs coordonnees reelles.
    gff_filename = assembly_name + "_genomic.gff.gz"
    gff_url = ftp_path + "/" + gff_filename
    gff_output_path = os.path.join(output_dir, gff_filename)

    protein_path = download_and_unzip(protein_url, protein_output_path)
    gff_path = download_and_unzip(gff_url, gff_output_path)

    print("Proteins:", protein_path)
    print("GFF annotations:", gff_path)

    return protein_path, gff_path
