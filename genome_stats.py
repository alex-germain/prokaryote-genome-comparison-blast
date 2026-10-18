# Alexandre GERMAIN

import csv


def compute_blast_stats(blast_result_file):
    """
    Statistiques descriptives du fichier BLAST tabulaire.

    Trois quantités differentes sont comptées, et il est important de ne pas
    les confondre :

      - nb_hits          : nombre total d'alignements retenus par BLAST.
                           Une meme protéine peut en produire plusieurs.
      - nb_proteins_gen1 : nombre de protéines distinctes du génome 1 qui ont
                           au moins un alignement, donc conservées.
      - nb_proteins_gen2 : idem pour le génome 2.

    Les pourcentages d'identité sont calculés sur le meilleur alignement de
    chaque protéine du génome 1 (le plus fort bitscore), pour qu'une protéine
    a nombreux paralogues ne pese pas plus qu'une autre dans la moyenne.
    """
    best_hits = {}
    proteins_gen1 = set()
    proteins_gen2 = set()
    nb_alignments = 0

    with open(blast_result_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            qseqid = row[0]
            sseqid = row[1]
            pident = float(row[2])
            bitscore = float(row[11])

            nb_alignments += 1
            proteins_gen1.add(qseqid)
            proteins_gen2.add(sseqid)

            if qseqid not in best_hits or bitscore > best_hits[qseqid]["bitscore"]:
                best_hits[qseqid] = {
                    "sseqid": sseqid,
                    "pident": pident,
                    "bitscore": bitscore,
                }

    if not best_hits:
        return None

    pident_values = [hit["pident"] for hit in best_hits.values()]

    return {
        "pident_mean": round(sum(pident_values) / len(pident_values), 2),
        "pident_min": round(min(pident_values), 2),
        "pident_max": round(max(pident_values), 2),
        "nb_hits": nb_alignments,
        "nb_proteins_gen1": len(proteins_gen1),
        "nb_proteins_gen2": len(proteins_gen2),
    }


def get_interpretation(pident_mean):
    """Lecture en langage courant du pourcentage d'identité moyen."""
    if pident_mean >= 90:
        return (f"Mean identity {pident_mean}% - the two genomes are very close. "
                "They are most likely strains of the same species.")
    elif pident_mean >= 70:
        return (f"Mean identity {pident_mean}% - the two genomes are moderately close. "
                "They may be species of the same genus.")
    elif pident_mean >= 50:
        return (f"Mean identity {pident_mean}% - the two genomes are relatively distant. "
                "They are most likely species of different genera.")
    else:
        return (f"Mean identity {pident_mean}% - the two genomes are very distant. "
                "They are most likely organisms of different families.")


def compute_cog_distribution(cog_gen):
    """
    Nombre de protéines par catégorie COG pour un génome.

    Renvoie {lettre_de_categorie: nombre}. La clé est la lettre seule, pas un
    libellé : c'est elle qui permet au module d'affichage d'associer la même
    couleur a une catégorie dans le dotplot et dans les camemberts. Le libellé
    lisible est reconstruit a l'affichage par cog_style.label_for().
    """
    distribution = {}

    for info in cog_gen.values():
        category = info["category"]
        distribution[category] = distribution.get(category, 0) + 1
    return distribution


def compute_common_cog(cog_gen1, cog_gen2, blast_data):
    """
    Nombre de paires de protéines alignées qui partagent la même catégorie COG.

    Renvoie {lettre_de_categorie: nombre}, dans le même format que
    compute_cog_distribution(), pour que les trois camemberts soient
    directement comparables.

    Une paire (query, subject) n'est comptée qu'une fois, même si BLAST l'a
    rapportée sur plusieurs lignes.
    """
    distribution = {}
    seen_pair = set()

    for d in blast_data:
        pair = (d["qseqid"], d["sseqid"])
        if pair in seen_pair:
            continue
        seen_pair.add(pair)

        if pair[0] not in cog_gen1 or pair[1] not in cog_gen2:
            continue

        category1 = cog_gen1[pair[0]]["category"]
        category2 = cog_gen2[pair[1]]["category"]

        # On ne garde que les paires dont les deux protéines sont rangées dans
        # la même catégorie fonctionnelle : ce sont les fonctions conservées.
        if category1 == category2:
            distribution[category1] = distribution.get(category1, 0) + 1
    return distribution


def interpret_cog(dist1, dist2, dist_common, name1, name2, top_n=3):
    """
    Lecture automatique des trois distributions COG, en anglais.

    Trois choses sont dites, et seulement celles-la, parce que ce sont les
    seules que les donnees permettent d'affirmer :

      1. quelles categories fonctionnelles dominent parmi les proteines
         conservees : ce sont les fonctions les mieux partagees par les deux
         organismes ;
      2. quelle categorie est la plus sur-representee chez l'un par rapport
         a l'autre, en points de pourcentage : c'est la piste d'une
         specialisation propre a une lignee ;
      3. la part du proteome de chaque genome qui a effectivement recu une
         annotation COG, car tout le reste est calcule sur cette part-la.

    Renvoie une liste de paragraphes. Les descriptions de categories sont
    importees de cog_style pour rester identiques partout.
    """
    import cog_style

    paragraphs = []

    total_common = sum(dist_common.values())
    total1 = sum(dist1.values())
    total2 = sum(dist2.values())

    if not total_common:
        return ["No protein pair received the same COG category in both "
                "genomes, so no functional interpretation can be given. This "
                "usually means the COG analysis was run on only one of the "
                "two genomes."]

    # 1. Fonctions les mieux conservees
    ranked = sorted(dist_common.items(), key=lambda kv: kv[1], reverse=True)
    pieces = [
        f"{cog_style.description_for(cat)} ({cat}, {count / total_common * 100:.1f}%)"
        for cat, count in ranked[:top_n]
    ]
    # Categories dites "de menage" : les fonctions essentielles, attendues en
    # tete d'un classement de proteines conservees. On ne commente leur
    # presence que si elles y sont reellement, plutot que d'affirmer une
    # attente que les donnees peuvent contredire.
    HOUSEKEEPING = {"J", "L", "C", "O", "F", "K"}
    top_categories = [cat for cat, _ in ranked[:top_n]]

    commentaire = (
        " Essential functions are expected to lead this ranking, since they "
        "are under the strongest selective pressure."
        if HOUSEKEEPING & set(top_categories)
        else " None of these are core housekeeping categories, which is worth "
             "noting: it may reflect a genuine shared specialisation, or simply "
             "that the COG annotation covers a small fraction of the proteome."
    )

    paragraphs.append(
        f"Among the {total_common} aligned protein pairs that share a COG "
        f"category, the most represented functions are: " + "; ".join(pieces) + ". "
        "These are the functions most conserved between the two organisms." +
        commentaire
    )

    # 2. Categorie la plus desequilibree entre les deux genomes
    if total1 and total2:
        ecarts = []
        for cat in set(dist1) | set(dist2):
            part1 = dist1.get(cat, 0) / total1 * 100
            part2 = dist2.get(cat, 0) / total2 * 100
            ecarts.append((part1 - part2, cat, part1, part2))

        ecart, cat, part1, part2 = max(ecarts, key=lambda e: abs(e[0]))
        if abs(ecart) >= 1.0:
            riche, pauvre = (name1, name2) if ecart > 0 else (name2, name1)
            paragraphs.append(
                f"The largest difference between the two proteomes is in "
                f"{cog_style.description_for(cat)} ({cat}): "
                f"{part1:.1f}% of annotated proteins in {name1} against "
                f"{part2:.1f}% in {name2}, a gap of {abs(ecart):.1f} points in "
                f"favour of {riche}. Such an imbalance often reflects an "
                f"expansion of this functional family in {riche}, or its loss "
                f"in {pauvre}."
            )
        else:
            paragraphs.append(
                "The two proteomes have very similar COG profiles: no "
                "category differs by more than one percentage point between "
                "them."
            )

    # 3. Couverture de l'annotation
    paragraphs.append(
        f"These proportions are computed on annotated proteins only: "
        f"{total1} proteins received a COG category in {name1} and {total2} in "
        f"{name2}. Proteins with no match in the COG database, and those whose "
        f"two members fall into different categories, are excluded, so the "
        f"figures above describe the annotated fraction rather than the whole "
        f"proteome."
    )

    return paragraphs
