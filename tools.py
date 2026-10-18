# Alexandre GERMAIN

"""
Fonctions utilitaires partagees par les modules qui lancent des programmes
externes ou qui construisent des noms de fichiers.

Deux problemes sont traites ici :

1. Le nom des binaires BLAST n'est pas le meme partout. Sur Debian et Ubuntu,
   le paquet ncbi-blast+ installe rpsblast sous le nom "rpsblast+", parce que
   "rpsblast" appartenait deja a l'ancien paquet blast. Sur Windows, macOS,
   conda et bioconda, le binaire s'appelle "rpsblast". Un nom code en dur rend
   donc le projet dependant d'une seule famille de distributions.

2. Les caracteres autorises dans un nom de fichier dependent du systeme. Sous
   Linux tout passe sauf "/" et l'octet nul ; sous Windows, les caracteres
   < > : " / \\ | ? * sont interdits. Les noms d'organismes NCBI contiennent
   regulierement des parentheses, des apostrophes, des crochets et parfois des
   deux-points.
"""

import os
import re
import shutil
import time


class BlastToolNotFound(RuntimeError):
    """Un programme BLAST requis est introuvable dans le PATH."""


# Pour chaque outil, les noms possibles dans l'ordre de preference.
BLAST_TOOL_NAMES = {
    "blastp": ["blastp"],
    "makeblastdb": ["makeblastdb"],
    "rpsblast": ["rpsblast", "rpsblast+"],
}

INSTALL_HINT = (
    "Installer BLAST+ :\n"
    "  conda install -c bioconda blast\n"
    "  Debian/Ubuntu : sudo apt install ncbi-blast+\n"
    "  Windows/macOS : https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
)


def find_blast_tool(tool):
    """
    Renvoie le chemin du binaire correspondant a `tool`.

    `tool` est une cle de BLAST_TOOL_NAMES ("blastp", "makeblastdb",
    "rpsblast"). Les noms possibles sont essayes dans l'ordre et le premier
    present dans le PATH est renvoye.

    Leve BlastToolNotFound si aucun n'est trouve, avec un message qui indique
    quoi installer plutot qu'un FileNotFoundError nu.
    """
    candidates = BLAST_TOOL_NAMES.get(tool, [tool])

    for name in candidates:
        path = shutil.which(name)
        if path is not None:
            return path

    noms = " ou ".join(candidates)
    raise BlastToolNotFound(
        f"{tool} est introuvable (cherche sous : {noms}).\n{INSTALL_HINT}"
    )


def check_blast_tools(tools=("blastp", "makeblastdb", "rpsblast")):
    """
    Verifie la presence de plusieurs outils d'un coup.

    Renvoie un dict {outil: chemin} si tout est present. Leve
    BlastToolNotFound en listant *tous* les manquants, pour que l'utilisateur
    n'ait pas a relancer le programme trois fois de suite.
    """
    found, missing = {}, []

    for tool in tools:
        try:
            found[tool] = find_blast_tool(tool)
        except BlastToolNotFound:
            missing.append(tool)

    if missing:
        raise BlastToolNotFound(
            "Programmes BLAST introuvables : " + ", ".join(missing) + "\n" + INSTALL_HINT
        )
    return found


# Tout ce qui n'est pas lettre, chiffre, point, tiret ou underscore est remplace.
_ILLEGAL = re.compile(r"[^A-Za-z0-9._-]+")


def safe_filename(name, max_length=100):
    """
    Transforme un nom quelconque en nom de fichier valide sur tous les systemes.

    Exemple : "Escherichia coli O157:H7 (strain Sakai)" devient
    "Escherichia_coli_O157_H7_strain_Sakai".

    La troncature a max_length evite de depasser la limite de longueur de
    chemin de Windows, atteinte facilement avec des noms du type
    "blast_results_<organisme1>_VS_<organisme2>_evalue_<valeur>".
    """
    if name is None:
        return "sans_nom"

    cleaned = _ILLEGAL.sub("_", str(name))
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")

    if not cleaned:
        return "sans_nom"
    return cleaned[:max_length].rstrip("._") or "sans_nom"


def is_valid_result(path, min_size=1):
    """
    Un fichier de resultat est reutilisable s'il existe et n'est pas vide.

    Sert a ne pas reservir le fichier laisse par un BLAST interrompu : son
    existence seule ne prouve pas qu'il est complet.
    """
    return os.path.exists(path) and os.path.getsize(path) >= min_size


# --------------------------------------------------------------------------
# Ecriture atomique
# --------------------------------------------------------------------------
#
# Tout le projet decide de refaire ou non un travail couteux en testant
# os.path.exists(). Ce test n'est fiable que si un fichier n'apparait a son
# emplacement definitif qu'une fois complet. Sinon un telechargement coupe
# au milieu laisse un fichier tronque que toutes les executions suivantes
# considereront comme valide, et l'erreur se manifestera bien plus loin,
# sous la forme d'un gzip illisible ou d'un FASTA ampute.
#
# La parade est d'ecrire dans un fichier temporaire ".part" puis de le
# renommer. os.replace() est atomique sur un meme systeme de fichiers : le
# nom final designe l'ancien fichier ou le nouveau, jamais un fichier a
# moitie ecrit.

PART_SUFFIX = ".part"


def atomic_write(output_path, write_callback):
    """
    Ecrit `output_path` via un fichier temporaire, puis le renomme.

    `write_callback` recoit le chemin temporaire et doit y ecrire le contenu.
    Si elle leve une exception, le fichier temporaire est supprime et
    `output_path` reste tel qu'il etait : jamais de fichier a moitie ecrit.
    """
    directory = os.path.dirname(output_path)
    if directory:
        os.makedirs(directory, exist_ok=True)

    part_path = output_path + PART_SUFFIX
    try:
        write_callback(part_path)
        os.replace(part_path, output_path)
    except BaseException:
        # BaseException et non Exception : un Ctrl-C pendant un telechargement
        # est precisement le cas que l'on veut nettoyer.
        if os.path.exists(part_path):
            os.remove(part_path)
        raise
    return output_path


def download_file(url, output_path, timeout=60, retries=3, chunk_size=8192):
    """
    Telecharge `url` vers `output_path`, de facon atomique et avec reessais.

    Trois protections, absentes d'un requests.get() nu :

    1. timeout : sans lui, requests attend indefiniment. Dans une interface
       graphique, le programme parait simplement fige.
    2. ecriture atomique : une coupure ne laisse pas de fichier tronque qui
       serait pris pour un telechargement reussi a l'execution suivante.
    3. reessais avec attente croissante : les serveurs FTP/HTTP du NCBI
       refusent regulierement une connexion quand ils sont charges, et le
       meme appel reussit quelques secondes plus tard.

    `timeout` s'applique a chaque etape du transfert, pas a sa duree totale :
    un gros genome n'est donc pas interrompu parce qu'il est long, seulement
    parce que le serveur ne repond plus.
    """
    import requests  # importe ici pour que tools reste utilisable sans requests

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            def write(part_path):
                with requests.get(url, stream=True, timeout=timeout) as response:
                    response.raise_for_status()
                    with open(part_path, "wb") as f:
                        for chunk in response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                f.write(chunk)

            return atomic_write(output_path, write)

        except requests.RequestException as error:
            last_error = error
            if attempt < retries:
                attente = 2 ** attempt  # 2 s, puis 4 s
                print(f"Echec du telechargement ({error}). "
                      f"Nouvel essai dans {attente} s "
                      f"[{attempt}/{retries - 1}]")
                time.sleep(attente)

    raise RuntimeError(
        f"Telechargement impossible apres {retries} tentatives :\n"
        f"  {url}\n"
        f"  derniere erreur : {last_error}"
    )
