# Alexandre GERMAIN

import csv
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak, KeepTogether)
from reportlab.lib.enums import TA_CENTER, TA_LEFT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tempfile
import genome_stats as gs
import dotplot as dp
import cog_style

# Objectif : 
# Fonctions qui sont capables de créer des fichiers csv et python qui contiennent les informations 
# importantes/ qui résument la comparaison des génomes 

# Le code couleur et les descriptions sont definis une seule fois, dans
# cog_style. Ces deux fonctions sont conservees pour ne pas casser le code
# appelant, mais elles ne font plus que relayer.
# Cote des dotplots dans le PDF. Chaque figure forme un bloc insecable avec
# son titre : titre (~0,8 cm) + figure + espacement (0,4 cm). La hauteur utile
# d'une page A4 avec les marges de ce rapport est de 25,7 cm, donc deux blocs
# tiennent tant que la figure ne depasse pas 12,2 cm. On prend 10,5 cm pour
# garder une marge confortable et ne pas dependre d'un calcul au millimetre.
DOTPLOT_SIZE = 10.5*cm


def get_cog_colors():
    return cog_style.COG_COLORS

def get_cog_description():
    return cog_style.COG_DESCRIPTIONS

def export_csv(blast_result_file, output_path, cog_gen1=None, cog_gen2=None): #Exporte les résultats BLAST (+ COG si disponible) dans un fichier CSV
    
    # On relit le fichier BLAST directement
    rows = []
    with open(blast_result_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            qseqid  = row[0]
            sseqid  = row[1]
            pident  = row[2]
            length  = row[3]
            evalue  = row[10]
            bitscore = row[11]

            cog1_id  = cog_gen1.get(qseqid, {}).get("cog_id",  "NA") if cog_gen1 else "NA"
            cog2_id  = cog_gen2.get(sseqid, {}).get("cog_id",  "NA") if cog_gen2 else "NA"
            cat1     = cog_gen1.get(qseqid, {}).get("category","NA") if cog_gen1 else "NA"
            cat2     = cog_gen2.get(sseqid, {}).get("category","NA") if cog_gen2 else "NA"
            func1    = cog_gen1.get(qseqid, {}).get("function_name","NA") if cog_gen1 else "NA"
            func2    = cog_gen2.get(sseqid, {}).get("function_name","NA") if cog_gen2 else "NA"

            rows.append([qseqid, sseqid, pident, length, evalue, bitscore,
                         cog1_id, cog2_id, cat1, cat2, func1, func2])               # contenu de chaque ligne du fichier csv

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                         "COG1", "COG2", "category1", "category2", "function1", "function2"])
        writer.writerows(rows)

    print(f"CSV exported: {output_path}")
    return output_path


def make_dotplot_image(blast_result_file, gff_gen1, gff_gen2,
                       cog_gen1=None, cog_gen2=None, mode="simple",
                       pident_threshold=0, name1=None, name2=None):
    #Génère une image PNG du dotplot (simple, synteny ou cog)
    
    data = dp.get_dotplot_data(blast_result_file, gff_gen1, gff_gen2, cog_gen1, cog_gen2)
    data = [d for d in data if float(d["pident"]) >= pident_threshold]

    if not data:
        return None

    cog_colors = get_cog_colors()

    x = [d["x"] for d in data]
    y = [d["y"] for d in data]

    if mode == "cog" and cog_gen1 is not None:
        c = [cog_colors.get(d["cog_category"], "black") for d in data]
    else:
        c = "black"

    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.scatter(x, y, s=1, c=c, alpha=0.7, linewidths=0)

    if mode == "synteny":
        diagonal_blocks = dp.detect_diagonal_blocks(data, bin_size=50000, min_points=30)
        for block in diagonal_blocks[:10]:
            x_block = [data[i]["x"] for i in block["indices"]]
            y_block = [data[i]["y"] for i in block["indices"]]
            ax.scatter(x_block, y_block, s=3, c="red", alpha=0.7, linewidths=0)
    
    if mode =="simple":
        ax.set_title("Simple dotplot", fontsize=9)
    elif mode =="synteny":
        ax.set_title("Dotplot with synteny regions", fontsize=9)
    elif mode=="cog":
        ax.set_title("Dotplot coloured by COG category", fontsize=9)

    ax.set_xlabel(f"Position on {name1}" if name1 else "Position on genome 1", fontsize=8)
    ax.set_ylabel(f"Position on {name2}" if name2 else "Position on genome 2", fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()

    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(tmp.name, dpi=150)
    plt.close(fig)
    return tmp.name


def compute_cog_distributions(cog_gen1, cog_gen2, blast_result_file, gff_gen1, gff_gen2):
    """Les trois distributions COG, calculees une seule fois par rapport."""
    dist1 = gs.compute_cog_distribution(cog_gen1)
    dist2 = gs.compute_cog_distribution(cog_gen2)
    blast_data = dp.get_dotplot_data(blast_result_file, gff_gen1, gff_gen2)
    dist_common = gs.compute_common_cog(cog_gen1, cog_gen2, blast_data)
    return dist1, dist2, dist_common


def make_cog_pie_images(dist1, dist2, dist_common, name1, name2):
    #Génère une image PNG avec les 3 camemberts COG.

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=100)

    # Meme fonction de dessin que dans l'interface, donc memes couleurs et
    # meme ordre de categories entre le PDF et l'ecran.
    cog_style.draw_pie(axes[0][0], dist1, name1, short_labels=True)
    cog_style.draw_pie(axes[0][1], dist2, name2, short_labels=True)

    # Camembert commun centré en bas
    axes[1][0].set_visible(False)
    axes[1][1].set_visible(False)
    cog_style.draw_pie(fig.add_subplot(2, 1, 2), dist_common, "Shared proteins")

    fig.tight_layout()

    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(tmp.name, dpi=100)
    plt.close(fig)
    return tmp.name

def make_cog_legend_table():
    cog_colors = get_cog_colors()
    descriptions = get_cog_description()

    table_data = [["Category", "Colour", "Description"]]

    for category, description in descriptions.items():
        table_data.append([category, "", description])
    table = Table(table_data, colWidths=[2*cm, 2*cm, 11*cm])
    style = [
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#4472C4")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("FONTSIZE", (0,0), (-1,-1), 7),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]
    for i, category in enumerate(descriptions.keys(), start=1):
        style.append(("BACKGROUND", (1,i), (1,i), colors.toColor(cog_colors.get(category, "black"))))
    
    table.setStyle(TableStyle(style))
    return table

def export_pdf(output_path, name1, name2, stats, interpretation, blast_result_file, gff_gen1, gff_gen2, cog_gen1=None, cog_gen2=None):
    #Génère un PDF complet avec stats, dotplot et camemberts COG.
    
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle("title",
                                  fontSize=16, alignment=TA_CENTER,
                                  spaceAfter=15, fontName="Helvetica-Bold")
    section_style = ParagraphStyle("section",
                                    fontSize=12, spaceAfter=8,
                                    fontName="Helvetica-Bold", textColor=colors.HexColor("#333333"))
    body_style = ParagraphStyle("body",
                                 fontSize=10, spaceAfter=7,
                                 fontName="Helvetica")

    doc = SimpleDocTemplate(output_path, pagesize=A4,
                            rightMargin=2*cm, leftMargin=2*cm,
                            topMargin=2*cm, bottomMargin=2*cm)

    content = []

    # ----- Titre -----
    content.append(Paragraph("Genome comparison", title_style))
    content.append(Paragraph(f"{name1} vs {name2}", title_style))
    content.append(Spacer(1, 0.5*cm))

    # ----- Stats BLAST -----
    content.append(Paragraph("BLAST statistics", section_style))

    table_data = [
        ["Metric", "Value"],
        ["Mean identity",                  f"{stats['pident_mean']} %"],
        ["Minimum identity",               f"{stats['pident_min']} %"],
        ["Maximum identity",               f"{stats['pident_max']} %"],
        ["Total number of alignments",     str(stats['nb_hits'])],
        [f"Conserved proteins ({name1})",  str(stats['nb_proteins_gen1'])],
        [f"Conserved proteins ({name2})",  str(stats['nb_proteins_gen2'])],
    ]

    table = Table(table_data, colWidths=[10*cm, 6*cm])
    table.setStyle(TableStyle([
        ("BACKGROUND",  (0, 0), (-1, 0), colors.HexColor("#4472C4")),
        ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
        ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",    (0, 0), (-1, -1), 10),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#EBF1FA")]),
        ("GRID",        (0, 0), (-1, -1), 0.5, colors.grey),
        ("ALIGN",       (0, 0), (-1, -1), "LEFT"),
        ("PADDING",     (0, 0), (-1, -1), 6),
    ]))
    content.append(table)
    content.append(Spacer(1, 1*cm))

    # ----- Interprétation -----
    content.append(Paragraph("Interpretation", section_style))
    content.append(Paragraph(interpretation, body_style))
    content.append(Spacer(1, 1*cm))

    # ----- Dotplot -----
    temp_files=[]
    dotplot_simple = make_dotplot_image(blast_result_file, gff_gen1, gff_gen2,
                                        cog_gen1, cog_gen2, mode="simple",
                                        name1=name1, name2=name2)
    dotplot_synteny = make_dotplot_image(blast_result_file, gff_gen1, gff_gen2,
                                         cog_gen1, cog_gen2, mode="synteny",
                                         name1=name1, name2=name2)

    dotplot_cog = None
    if cog_gen1 is not None and cog_gen2 is not None:
        dotplot_cog = make_dotplot_image(blast_result_file, gff_gen1, gff_gen2,
                                         cog_gen1, cog_gen2, mode="cog",
                                         name1=name1, name2=name2)
    
    for title, img_path in [
        ("Simple dotplot", dotplot_simple),
        ("Dotplot with synteny regions", dotplot_synteny),
        ("Dotplot coloured by COG category", dotplot_cog)]:
        if img_path:
            temp_files.append(img_path)
            # KeepTogether empeche reportlab de laisser le titre en bas d'une
            # page et la figure en haut de la suivante : le bloc entier passe
            # a la page suivante s'il ne tient pas.
            content.append(KeepTogether([
                Paragraph(title, section_style),
                Image(img_path, width=DOTPLOT_SIZE, height=DOTPLOT_SIZE),
                Spacer(1, 0.4*cm),
            ]))

    # ----- Camemberts COG -----
    pie_img_path = None

    if cog_gen1 is None or cog_gen2 is None:
        content.append(Paragraph("COG category distribution", section_style))
        content.append(Paragraph(
            "Functional analysis was not run. Use the 'COG colour code' button "
            "to include this section.",
            body_style))
    else:
        content.append(PageBreak())
        content.append(Paragraph("COG category distribution", section_style))

        dist1, dist2, dist_common = compute_cog_distributions(
            cog_gen1, cog_gen2, blast_result_file, gff_gen1, gff_gen2)

        pie_img_path = make_cog_pie_images(dist1, dist2, dist_common, name1, name2)
        if pie_img_path:
            content.append(Image(pie_img_path, width=16*cm, height=13*cm))
            content.append(Spacer(1, 0.4*cm))

        # Lecture automatique des camemberts, calculee sur les distributions
        # ci-dessus : fonctions les mieux conservees, categorie la plus
        # desequilibree entre les deux proteomes, couverture de l'annotation.
        for paragraph in gs.interpret_cog(dist1, dist2, dist_common, name1, name2):
            content.append(Paragraph(paragraph, body_style))

        # La legende des couleurs n'apparait qu'une fois, en fin de rapport :
        # elle servait a la fois aux dotplots COG et aux camemberts, et etait
        # auparavant inseree deux fois.
        # Le tableau mesure environ 18 cm pour 26 cm de hauteur utile : il
        # tient donc sur une page, et KeepTogether garantit qu'il n'est ni
        # coupe en deux ni separe de son titre.
        content.append(Spacer(1, 0.4*cm))
        content.append(KeepTogether([
            Paragraph("COG colour code", section_style),
            make_cog_legend_table(),
        ]))

    # Construction du PDF
    doc.build(content)

    # Nettoyage des fichiers temporaires
    for tmp_path in temp_files + [pie_img_path]:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

    print(f"PDF exported: {output_path}")
    return output_path