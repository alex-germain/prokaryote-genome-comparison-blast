# Alexandre GERMAIN

"""
Source unique de vérité pour l'affichage des catégories COG.

Avant ce module, le code couleur COG était défini deux fois : une fois dans
app.get_cog_colors() et une fois dans
export_results.get_cog_colors(). Les deux copies avaient déjà divergé
(la categorie du bas valait "lightgray" dans l'une et
"#D3D3D3" dans l'autre). Tout ce qui concerne le rendu d'une catégorie COG
(ordre, couleur, libellé) est donc regroupé ici et importé par les modules
d'affichage.

Le choix des couleurs n'est pas arbitraire : chaque grand groupe COG a sa
famille de teintes, ce qui permet de reconnaître le groupe d'un point du
dotplot même sans lire la légende.

    Information storage and processing  -> bleus
    Cellular processes and signaling    -> verts et turquoises
    Metabolism                          -> oranges, rouges et bruns
    Poorly characterized                -> violets
    Autre                               -> gris

Les 26 couleurs sont données en hexadécimal et non par leur nom
(matplotlib et Tk acceptent les deux, mais les noms rapprochaient des
teintes indistinguables : "blue"/"darkblue", "turquoise"/"darkturquoise",
"coral"/"tomato").
"""

# Catégories groupées, dans l'ordre officiel des groupes fonctionnels COG.
# C'est cet ordre qui fixe l'ordre de la légende et des camemberts.
COG_GROUPS = [
    ("Information storage and processing", [
        ("J", "Translation, ribosomal structure and biogenesis", "#1f4e9c"),
        ("A", "RNA processing and modification",                 "#3b7dd8"),
        ("K", "Transcription",                                   "#6ba3e8"),
        ("L", "Replication, recombination and repair",           "#0d2d5e"),
        ("B", "Chromatin structure and dynamics",                "#9cc5f0"),
    ]),
    ("Cellular processes and signaling", [
        ("D", "Cell cycle control, cell division, chromosome partitioning", "#0b6b4f"),
        ("Y", "Nuclear structure",                               "#17a673"),
        ("V", "Defense mechanisms",                              "#5ed0a4"),
        ("T", "Signal transduction mechanisms",                  "#0f8a8a"),
        ("M", "Cell wall/membrane/envelope biogenesis",          "#45b8c4"),
        ("N", "Cell motility",                                   "#8fd4b0"),
        ("Z", "Cytoskeleton",                                    "#006d5b"),
        ("W", "Extracellular structures",                        "#b7e3cf"),
        ("U", "Intracellular trafficking, secretion, and vesicular transport", "#2f8f5b"),
        ("O", "Posttranslational modification, protein turnover, chaperones",  "#7fc241"),
        ("X", "Mobilome: prophages, transposons",                "#c3e08a"),
    ]),
    ("Metabolism", [
        ("C", "Energy production and conversion",                "#d94f1e"),
        ("E", "Amino acid transport and metabolism",             "#f07d2a"),
        ("F", "Nucleotide transport and metabolism",             "#f5a94e"),
        ("G", "Carbohydrate transport and metabolism",           "#a8321e"),
        ("H", "Coenzyme transport and metabolism",               "#c96a3c"),
        ("I", "Lipid transport and metabolism",                  "#e8c547"),
        ("P", "Inorganic ion transport and metabolism",          "#8c4a1f"),
        ("Q", "Secondary metabolites biosynthesis, transport and catabolism", "#f2d08a"),
    ]),
    ("Poorly characterized", [
        ("R", "General function prediction only",                "#7a5aa8"),
        ("S", "Function unknown",                                "#b09cd0"),
    ]),
    ("Autre", [
        ("Mismatch", "Different or unknown category", "#c8c8c8"),
    ]),
]

# Couleur utilisée pour une catégorie absente des tables COG.
DEFAULT_COLOR = "#000000"

# Dictionnaires plats construits une seule fois a partir de COG_GROUPS,
# pour que l'ordre et le contenu ne puissent pas se desynchroniser.
COG_COLORS = {cat: color for _, cats in COG_GROUPS for cat, _, color in cats}
COG_DESCRIPTIONS = {cat: desc for _, cats in COG_GROUPS for cat, desc, _ in cats}

# Ordre d'affichage canonique des catégories.
COG_ORDER = [cat for _, cats in COG_GROUPS for cat, _, _ in cats]


def color_for(category):
    """Couleur d'une catégorie COG, noir si la catégorie est inconnue."""
    return COG_COLORS.get(category, DEFAULT_COLOR)


def colors_for(categories):
    """Liste de couleurs alignée sur une liste de catégories."""
    return [color_for(category) for category in categories]


def description_for(category):
    """Description longue d'une catégorie COG."""
    return COG_DESCRIPTIONS.get(category, "Catégorie inconnue")


def label_for(category):
    """Libellé de légende : "J : Translation, ribosomal structure..."."""
    return f"{category} : {description_for(category)}"


def sort_categories(categories):
    """
    Trie des catégories dans l'ordre COG officiel.

    Les catégories inconnues sont placées a la fin, par ordre alphabétique,
    afin qu'une catégorie inattendue reste visible au lieu d'etre perdue au
    milieu de la légende.
    """
    known = [c for c in COG_ORDER if c in categories]
    unknown = sorted(c for c in categories if c not in COG_COLORS)
    return known + unknown


def draw_pie(ax, distribution, title, fontsize=6, short_labels=False):
    """
    Dessine un camembert COG sur l'axe matplotlib `ax`.

    `distribution` est un dictionnaire {lettre_de_categorie: nombre}, tel que
    le renvoient genome_stats.compute_cog_distribution() et
    genome_stats.compute_common_cog().

    Trois choses sont garanties ici, et c'est la raison d'etre de la fonction :

      - les parts sont triées dans l'ordre COG officiel, donc les camemberts
        de deux génomes se lisent dans le meme sens ;
      - chaque part recoit la couleur de sa catégorie, la meme que dans le
        dotplot COG et dans la légende ;
      - la légende donne la description complete de la catégorie et son
        pourcentage.

    Renvoie l'axe, ou None si la distribution est vide (cas d'un génome sans
    aucune protéine annotée, qui ferait planter ax.pie()).
    """
    if not distribution:
        ax.set_visible(False)
        return None

    categories = sort_categories(distribution.keys())
    values = [distribution[c] for c in categories]
    total = sum(values)

    wedges, _ = ax.pie(
        values,
        colors=colors_for(categories),
        startangle=90,
        counterclock=False,
        wedgeprops={"linewidth": 0.3, "edgecolor": "white"},
    )
    ax.set_title(title, fontsize=fontsize + 3)

    # Les descriptions completes sont longues : cote a cote, deux legendes
    # completes sortent de la figure. short_labels=True n'affiche que la
    # lettre, la description restant disponible dans la legende COG.
    if short_labels:
        labels = [f"{c} ({round(v / total * 100, 1)}%)"
                  for c, v in zip(categories, values)]
    else:
        labels = [f"{label_for(c)} ({round(v / total * 100, 1)}%)"
                  for c, v in zip(categories, values)]
    ax.legend(
        wedges, labels,
        loc="center left", bbox_to_anchor=(1.0, 0.5),
        fontsize=fontsize, frameon=True,
    )
    return ax
