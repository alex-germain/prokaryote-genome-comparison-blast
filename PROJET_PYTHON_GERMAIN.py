# Interface graphique de comparaison de genomes par BLAST
# Alexandre Germain

import create_local_DB
import telechargement_genomes
import lancer_blast
import dotplot
import cog_analysis
import genome_stats
import export_results
import cog_style

import subprocess 
import os
import csv
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter import messagebox
import sys 
import threading # pour les "..." de chargement
import mplcursors

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

class Genomes():
    def __init__(self):

        # Le projet utilise partout des chemins relatifs ("genomes/",
        # "blast_results/", "./COG/Cog"...), qui sont resolus par rapport au
        # repertoire courant, pas au script. Lance depuis un autre dossier,
        # le programme creait donc ses fichiers ailleurs et ne retrouvait plus
        # la base COG. On se place une fois pour toutes dans le dossier du
        # projet.
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # Le chemin vient de create_local_DB : auparavant l'interface lisait
        # "data/procaryote_db.tsv" alors que create_DB() ecrit dans
        # "Local_DB/". Sur une installation neuve, la base etait donc
        # retelechargee a chaque lancement, puis introuvable.
        self.DB_Path = create_local_DB.DB_PATH
        if not os.path.exists(self.DB_Path):
            create_local_DB.create_DB()
        else :
            print("Local database already present:", self.DB_Path)

        self.genomes=[]
        self.gen1 = None
        self.gen2 = None
        self.blast_result_file=None
        self.dotplot_file=None
        self.current_plot_mode = "simple"  # "simple", "synteny", "cog"
        self.path_gen1 = None
        self.path_gen2 = None
        self.gff_gen1 = None
        self.gff_gen2 = None

        self.load_database()
        self.creeInterface()
        self.fenetre.mainloop()

    # On charge la base de donnée "procaryote_db.tsv" dans la liste self.genomes
    def load_database(self):
          with open(self.DB_Path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter = "\t")
            for row in reader :
                self.genomes.append(row)
                 
    def creeInterface(self):
        self.fenetre = tk.Tk()
        self.fenetre.title('Prokaryotic Genome Comparison')
        self.fenetre.geometry("1500x850")

        # ====================
        # STRUCTURE PRINCIPALE 
        # ====================
        self.left_panel = tk.Frame(self.fenetre, width = 420, bg="#f2f2f2")
        self.left_panel.pack(side=tk.LEFT, fill=tk.Y)
        self.left_panel.pack_propagate(False)

        self.right_panel = tk.Frame(self.fenetre, bg="white")
        self.right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ====================
        # PANNEAU GAUCHE
        # ====================

        title = tk.Label(self.left_panel, text="Genome Comparison", font=("Arial", 13,"bold"), bg="#f2f2f2")
        title.pack(pady=5)

        # ----- Génome 1 -----
        frame1 = tk.LabelFrame(self.left_panel, text="Genome 1", bg="#f2f2f2", padx=5, pady=1)
        frame1.pack(fill=tk.X, padx=8, pady=1)
        
        # Domaine du Génome 1
        tk.Label(frame1, text="Domain", bg="#f2f2f2").pack(anchor='w')
        self.taxonomy1 = ttk.Combobox(frame1, values=["", "Bacteria", "Archaea"], state="readonly")
        self.taxonomy1.pack(fill=tk.X)
        self.taxonomy1.bind("<<ComboboxSelected>>", self.update_list1)  # si une taxonomie est choisie, la listbox est mise à jour avec seulement les espèces de la taxonomie sélectionnée
        # Entrée du Génome 1
        self.search1 = tk.Entry(frame1)
        self.search1.insert(0, "Search species 1")
        self.search1.config(fg="gray")
        self.search1.pack(fill=tk.X, pady=2)
        self.search1.bind("<FocusIn>", self.clear_placeholder1)
        self.search1.bind("<FocusOut>", self.restore_placeholder1)
        self.search1.bind("<KeyRelease>", self.update_list1) # à chaque lettre du mot écrit, le search1 est mis à jour
        # Listbox des génomes proposés en fonction de ce qui est écrit dans search1
        self.listbox1 = tk.Listbox(frame1, height =5, width=80)
        self.listbox1.bind("<<ListboxSelect>>", self.select_genome1)

        # ----- Génome 2 -----
        frame2 = tk.LabelFrame(self.left_panel, text="Genome 2", bg="#f2f2f2", padx=5, pady=1)
        frame2.pack(fill=tk.X, padx=8, pady=1)

        # Domaine du Génome 2
        tk.Label(frame2, text="Domain", bg="#f2f2f2").pack(anchor="w")
        self.taxonomy2 = ttk.Combobox(frame2, values=["", "Bacteria", "Archaea"],state="readonly")
        self.taxonomy2.pack(fill=tk.X)
        self.taxonomy2.bind("<<ComboboxSelected>>", self.update_list2)
        # Entrée du Génome 2
        self.search2 = tk.Entry(frame2)
        self.search2.insert(0, "Search species 2")
        self.search2.config(fg="gray")
        self.search2.pack(fill=tk.X, pady=3)
        self.search2.bind("<FocusIn>", self.clear_placeholder2)
        self.search2.bind("<FocusOut>", self.restore_placeholder2)
        self.search2.bind("<KeyRelease>", self.update_list2)
        # Listbox du génome 2
        self.listbox2 = tk.Listbox(frame2, height =5, width = 80)
        self.listbox2.bind("<<ListboxSelect>>", self.select_genome2)
        
        # ----- Actions Principales -----

        action_frame=tk.LabelFrame(self.left_panel, text="Actions", bg="#f2f2f2", padx=5, pady=1)
        action_frame.pack(fill=tk.X, padx=8, pady=1)

        # Download Button

        self.download_button = tk.Button(action_frame, text = "Download genomes", command= self.download_selected_genomes)
        self.download_button.pack(fill=tk.X, pady=1) 

        self.download_status = tk.Label(action_frame, text="", bg="#f2f2f2",font=("Arial", 5), fg="#358642")
        self.download_status.pack() 

        # Evalue du BLAST
        self.evalue_label = tk.Label(action_frame,text="BLAST e-value", bg="#d9d9d9",anchor="w",padx=150, font=("Arial", 9))
        self.evalue_label.pack(fill=tk.X, pady=(3,2))

        evalue_frame = tk.Frame(action_frame, bg="#f2f2f2")
        evalue_frame.pack(pady=1)

        self.evalue_entry = tk.Entry(evalue_frame, width=15, font=("Arial", 9))
        self.evalue_entry.pack(side=tk.LEFT)
        self.evalue_entry.insert(0,"1e-20")                             # la evalue est fixée par défaut à 1e-20
        self.evalue_entry.bind("<FocusOut>", self.restore_default_evalue)

        self.default_label = tk.Label(evalue_frame, text="(default)", bg="#f2f2f2", fg="gray", font=("Arial",9))
        self.default_label.pack(side=tk.LEFT, padx=1)

        # Blast Button
        self.blast_button = tk.Button(action_frame, text = "Run BLAST", command =self.launch_blast ) # Bouton pour lancer BLAST
        self.blast_button.pack(fill=tk.X, pady=1)

        self.blast_info_frame = tk.Frame(action_frame, bg="#f2f2f2")                                   # Indique si Blast est en cours ou terminé
        self.blast_info_frame.pack(pady=2)

        self.loading_label = tk.Label(self.blast_info_frame, text="", bg="#f2f2f2", font=("Arial", 5)) # les "..." de chargement pendant que BLAST tourne
        self.loading_label.pack(side=tk.LEFT, padx=3)

        self.blast_status = tk.Label(self.blast_info_frame, text="", bg="#f2f2f2", font=("Arial", 5))
        self.blast_status.pack(side=tk.LEFT, padx=3)

        # ----- Analyse/Affichage ------
        analysis_frame =tk.LabelFrame(self.left_panel, text='Visualisation', bg="#f2f2f2", padx=5, pady=1)
        analysis_frame.pack(fill=tk.X, padx=8, pady=1)

        # Dotplot Button (simple, synteny and cog)
        self.dotplot_button = tk.Button(analysis_frame, text="Show dotplot", command=self.show_basic_dotplot)
        self.dotplot_button.pack(fill=tk.X, pady=1)

        self.synteny_button = tk.Button(analysis_frame, text="Show synteny regions", command=self.show_synteny_dotplot)
        self.synteny_button.pack(fill=tk.X, pady=1)

        self.cog_button = tk.Button(analysis_frame, text="COG colour code", command=self.show_cog_dotplot)
        self.cog_button.pack(fill=tk.X, pady=1)

        # Filtre Seuil Identité
        filter_frame = tk.Frame(analysis_frame, bg='#f2f2f2')
        filter_frame.pack(fill=tk.X, pady=1)

        tk.Label(filter_frame, text="Min. identity threshold:", bg="#f2f2f2", font=("Arial", 10)).pack(side=tk.LEFT, padx=5)

        self.pident_filter_entry = tk.Entry(filter_frame, width=6)
        self.pident_filter_entry.insert(0, "0")
        self.pident_filter_entry.pack(side=tk.LEFT)

        tk.Label(filter_frame, text="%", bg="#f2f2f2", font=("Arial",10)).pack(side=tk.LEFT)

        self.apply_filter_button = tk.Button(filter_frame, text="Apply", command=self.apply_pident_filter)
        self.apply_filter_button.pack(side=tk.LEFT, padx=5)

        # Plus d'Infos Button
        self.info_button = tk.Button(analysis_frame, text="More info", command=self.show_info_popup)
        self.info_button.pack(fill=tk.X, pady=3)

        # Export (CSV or PDF) Button
        self.export_button = tk.Button(analysis_frame, 
                                       text ="Save results", 
                                       command= self.save_results,
                                       font=("Arial", 10), pady=1)
        self.export_button.pack(fill=tk.X, pady=1)

        # ----- Zone de Statut -----
        self.status_label = tk.Label(self.left_panel, text="Ready.", bg="#f2f2f2", fg="gray", font=("Arial", 9), wraplength=380, justify="left")
        self.status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=8, pady=1)

        # =============
        # PANNEAU DROIT  (dotplot)
        # =============

        self.plot_title = tk.Label(self.right_panel, text = "No dotplot displayed", bg="white", font=("Arial",16,"bold"))
        self.plot_title.pack(pady=5)

        self.plot_container = tk.Frame(self.right_panel, bg="white")
        self.plot_container.pack(fill=tk.BOTH, expand=True)

    #  -----------------------
    #   Sélection des génomes 
    #  -----------------------
    def clear_placeholder1(self, event):                    # sert à éffacer "Search species 1" quand l'utilisateur écrit
        if self.search1.get()=="Search species 1":
            self.search1.delete(0, tk.END)
            self.search1.config(fg='black')
    def restore_placeholder1(self, event):                  # sert à ré-afficher "Search species 1" si rien n'est écrit par l'utilisateur
        if self.search1.get()=="":
            self.search1.insert(0, "Search species 1")
            self.search1.config(fg='gray')
    def clear_placeholder2(self, event):                    # sert à éffacer "Search species 2"
        if self.search2.get()=="Search species 2":
            self.search2.delete(0, tk.END)
            self.search2.config(fg='black')
    def restore_placeholder2(self, event):                  # sert à ré-afficher "Search species 2"
        if self.search2.get()=="":
            self.search2.insert(0, "Search species 2")
            self.search2.config(fg='gray')


    def search_genomes(self, query, taxonomy_name=""):      # Fonction qui permet de mettre à jour la listbox (en fonction de ce qu'écrit l'utilisateur)
        query = query.lower().strip()                       # le query correspond à ce qu'écrit l'utilisateur (lettre par lettre)

        if len(query)<1 and taxonomy_name == "":
            return[]
        
        results = []
        for genome in self.genomes :
            organism_name = genome["organism_name"].lower()
            
            if taxonomy_name!= "":
                if genome.get("domain","") != taxonomy_name:
                    continue
            
            if query != "" :
                if query not in organism_name:
                    continue 
            results.append(genome)
        return results[:200]                                # renvoie tous les noms de génomes qui contiennent le mot qu'écrit l'utilisateur. 
    
    def update_list1(self, event=None):                     # Mise à jour des listbox en affichant "organism_name|assembly_accession|assembly_level"
        query = self.search1.get().strip()
        taxonomy_name = self.taxonomy1.get().strip()

        results = self.search_genomes(query, taxonomy_name)

        self.search1.configure(bg="white")
        self.listbox1.delete(0, tk.END)

        if len(query)==0 and taxonomy_name=="":
            self.listbox1.forget()                          # on cache la listbox s'il n'y a rien d'écrit
            return
        self.listbox1.pack()

        for genome in results : 
            text = (genome["organism_name"]+"|"+genome["assembly_accession"]+"|"+genome["assembly_level"])
            self.listbox1.insert(tk.END, text)

        if len(results)==0 : 
            self.listbox1.pack_forget()                     # si on a aucun résultat à afficher, on cache la listbox

    def update_list2(self, event=None):                     # Idem pour la liste du Génome 2
        query = self.search2.get()
        taxonomy_name = self.taxonomy2.get()

        results = self.search_genomes(query, taxonomy_name)
        self.search2.configure(bg="white")
        self.listbox2.delete(0, tk.END)

        if len(query)==0 and taxonomy_name=="":
            self.listbox2.forget() # on cache la listbox s'il n'y a rien d'écrit
            return
        self.listbox2.pack()

        for genome in results : 
            text = (genome["organism_name"]+"|"+genome["assembly_accession"]+"|"+genome["assembly_level"])
            self.listbox2.insert(tk.END, text)

        if len(results)==0 : 
            self.listbox2.pack_forget() 

    def find_genome_by_accession(self, accession):          # Permet de retenir le génome choisi via son assembly_accession (unique) plutot
        for genome in self.genomes:                         # que par son nom (il peut il y avoir plusieurs séquences pour une même espèce)
            if genome["assembly_accession"]== accession:
                return genome
        return None
    
    def select_genome1(self, event=None):                   # Quand un génome est sélectionné dans la listbox
        selection = self.listbox1.curselection()
        if not selection: 
            return
        selected_text =self.listbox1.get(selection[0])
        accession = selected_text.split("|")[1].strip()
        self.gen1 = self.find_genome_by_accession(accession)  # on sauvegarde l'assembly_accession du génome sélectionné

        if self.gen1 :                                      # la listbox disparait et le nom est affiché dans la barre de recherche
            self.search1.delete(0, tk.END)
            self.search1.insert(0, self.gen1["organism_name"])
            self.search1.configure(bg="#D8F3DC")
            self.listbox1.pack_forget()
            print("Genome 1:", self.gen1["organism_name"])            
    
    def select_genome2(self, event=None):                   # Idem pour le génome 2
        selection = self.listbox2.curselection()
        if not selection: 
            return
        selected_text =self.listbox2.get(selection[0])
        accession = selected_text.split("|")[1].strip()
        self.gen2 = self.find_genome_by_accession(accession)

        if self.gen2 : 
            self.search2.delete(0, tk.END)
            self.search2.insert(0, self.gen2["organism_name"])
            self.search2.configure(bg="#D8F3DC")
            self.listbox2.pack_forget()
            print("Genome 2:", self.gen2["organism_name"])
        
    def download_selected_genomes(self):        # téléchargement des génomes sélectionnés à partir de leur ftp_path
        if self.gen1 is None or self.gen2 is None : 
            messagebox.showwarning("Incomplete selection", "Select two genomes before downloading.")
            return
        path_gen1, gff_gen1 = telechargement_genomes.download_genome(self.gen1)
        path_gen2, gff_gen2 = telechargement_genomes.download_genome(self.gen2)
        self.path_gen1 = path_gen1
        self.path_gen2 = path_gen2
        self.gff_gen1 = gff_gen1
        self.gff_gen2 = gff_gen2

        self.download_status.config(text='(Genomes downloaded)')

        print("Genome 1 downloaded:", path_gen1)
        print("Genome 2 downloaded:", path_gen2)

    #  -----------------------
    #          E-value 
    #  -----------------------
    def restore_default_evalue(self, event=None):           # Permettre à l'utilisateur de fixer le seuil de e-value qu'il veut pour le blast
        if self.evalue_entry.get().strip() =="":
            self.evalue_entry.insert(0,"1e-20")
            self.default_label.config(text="(default)")
        elif self.evalue_entry.get().strip()=="1e-20":
            self.default_label.config(text="(default)")
        else :
            self.default_label.config(text="")
    
    #  -----------------------
    #           BLAST 
    #  -----------------------
    def animate_loading(self, step=0):      # fonction pour animer les "..." pendant que BLAST tourne
        if not self.blast_running:
            self.loading_label.config(text="")
            return
        dots = "."*((step%3)+1)
        self.loading_label.config(text=f"BLAST running{dots}")
        self.fenetre.after(500, self.animate_loading, step+1)

    def launch_blast(self):                 # Lancer BLAST
        if not self.gen1 or not self.gen2 :
            messagebox.showerror("Error", "Please select two genomes and download them.")
            return
        self.blast_running = True
        self.blast_status.config(text="")
        self.animate_loading()
        thread = threading.Thread(target=self.run_blast_thread)
        thread.start()
    
    def run_blast_thread(self):             
        evalue = self.evalue_entry.get().strip()
        query_name= self.gen1["assembly_accession"]
        subject_name= self.gen2["assembly_accession"]

        try : 
            path_gen1, gff_gen1 = telechargement_genomes.download_genome(self.gen1)
            path_gen2, gff_gen2 = telechargement_genomes.download_genome(self.gen2)

            self.path_gen1 = path_gen1
            self.path_gen2 = path_gen2
            self.gff_gen1 = gff_gen1
            self.gff_gen2 = gff_gen2

            resultat = lancer_blast.run_blastp(path_gen1, path_gen2,query_name,subject_name, evalue=evalue)
            self.fenetre.after(0, self.blast_finished, resultat)

        except Exception as e:
            self.fenetre.after(0, self.blast_error, str(e))
    
    def blast_finished(self, resultat):             # Permet d'afficher quand BLAST est terminé
        self.blast_running=False
        self.loading_label.config(text="")
        self.blast_status.config(text="(BLAST finished)", fg="#358642")
        self.blast_result_file = resultat
        print("Result:", resultat)
    
    def blast_error(self, error_message):           # Permet d'afficher q'il y a une erreur pour laquelle BLAST ne marche pas
        self.blast_running=False
        self.loading_label.config(text="")
        messagebox.showerror("BLAST error", error_message)
    
    #  -----------------------
    #          DOTPLOT 
    #  -----------------------
    def clear_plot_area(self):                      # permet de vider l'espace pour chaque nouveau plot
        for widget in self.plot_container.winfo_children():
            widget.destroy()
    
    def display_dotplot_in_interface(self, data, title, colors=None, diagonal_blocks=None, show_cog_info=False):
        self.clear_plot_area()
        self.plot_title.config(text=title)

        fig = Figure(figsize=(8,8), dpi=100)
        ax = fig.add_subplot(111)

        x = [d["x"] for d in data]
        y = [d["y"] for d in data]

        if colors is None:
            colors="black"

        scatter = ax.scatter(x,y,s=3, c=colors, alpha =0.7, linewidths=0, picker = True)

        if diagonal_blocks is not None :
            for block in diagonal_blocks[:10]:
                x_block = [data[i]["x"] for i in block["indices"]]
                y_block = [data[i]["y"] for i in block["indices"]]
                ax.scatter(x_block, y_block, s=3, c='red', linewidths=0)

        ax.set_xlabel(f"Position on {self.gen1['organism_name']}")
        ax.set_ylabel(f"Position on {self.gen2['organism_name']}")
        ax.set_aspect('equal', adjustable='box')

        plot_and_toolbar = tk.Frame(self.plot_container, bg="white")
        plot_and_toolbar.pack(fill=tk.BOTH, expand=True)

        toolbar_frame = tk.Frame(self.plot_container, bg='white')
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        #canvas Matplotlib au dessus
        canvas = FigureCanvasTkAgg(fig, master = self.plot_container)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH, expand=True)
        #toolbar en bas
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame, pack_toolbar = False)
        toolbar.update()
        toolbar.pack(side=tk.LEFT, fill=tk.X)
        
        #Curseur interactif
        cursor = mplcursors.cursor(scatter)
        @cursor.connect("add")
        def on_add(sel):
            index=sel.index
            d=data[index]
            # Le bloc COG etait ajoute deux fois : une fois inconditionnellement
            # et une fois sous show_cog_info. Il n'apparait plus que lorsque
            # l'analyse COG a ete lancee, ce qui est le role du parametre.
            texte = (
                f"Protein 1: {d['qseqid']}\n"
                f"Protein 2: {d['sseqid']}\n"
                f"Identity: {d['pident']}%\n"
                f"Length: {d['length']}\n"
                f"E-value: {d['evalue']}\n"
                f"Bitscore: {d['bitscore']}")

            if show_cog_info:
                category = d['cog_category']
                texte += (
                    f"\n\nCOG 1: {d['cog1']}"
                    f"\nCOG 2: {d['cog2']}"
                    f"\nCOG category: {category}"
                    f"\n{cog_style.description_for(category)}"
                )
            sel.annotation.set_text(texte)

    def show_basic_dotplot(self):                   #Afficher le dotplot Simple
        if self.blast_result_file is None :
            messagebox.showwarning("Dotplot unavailable", "Please run BLAST first.")
            return
        
        self.current_plot_mode = "simple"
        self.status_label.config(text="Drawing simple dotplot...")
        self.fenetre.update_idletasks()
        
        self.dotplot_file = dotplot.create_dotplot(self.blast_result_file, self.gff_gen1, self.gff_gen2)

        if self.dotplot_file is None : 
            messagebox.showwarning("Dotplot unavailable", "No BLAST hit found.")
            return
        
        data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2)
        self.display_dotplot_in_interface(data=data,
                                          title = f"Dotplot: {self.gen1['organism_name']} vs {self.gen2['organism_name']}",
                                          colors = 'black',
                                          diagonal_blocks=None,
                                          show_cog_info=False)
        self.status_label.config(text="Simple dotplot displayed.")
    
    def show_synteny_dotplot(self):                 # Afficher le dotplot avec les régions de synténie
        if self.blast_result_file is None :
            messagebox.showwarning("Dotplot unavailable", "Please run BLAST first.")
            return

        self.current_plot_mode = "synteny"
        self.status_label.config(text="Detecting synteny regions...")
        self.fenetre.update_idletasks()

        data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2)
        diagonal_blocks = dotplot.detect_diagonal_blocks(data, bin_size=50000, min_points=30)

        self.display_dotplot_in_interface(data=data,
                                          title = f"Synteny: {self.gen1['organism_name']} vs {self.gen2['organism_name']}",
                                          colors = 'black',
                                          diagonal_blocks=diagonal_blocks,
                                          show_cog_info=False)
        self.status_label.config(text=f"Synteny regions displayed. {len(diagonal_blocks)} diagonal(s) detected.")
    
    def get_cog_colors(self):
        # Le code couleur est defini une seule fois, dans cog_style, et partage
        # par le dotplot, la legende, les camemberts et l'export PDF.
        return cog_style.COG_COLORS

    def show_cog_dotplot(self):                 # Afficher le dotplot COG
        if self.blast_result_file is None :
            messagebox.showwarning("Dotplot unavailable", "Please run BLAST first.")
            return

        self.current_plot_mode = "cog"
        self.status_label.config(text="Running COG analysis...")
        self.fenetre.update_idletasks()

        cog_db_path = "./COG/Cog"
        cognames_file = "cognames2003-2014.tab"
        fun_file = "fun2003-2014.tab"

        cog_file_gen1 = cog_analysis.run_rpsblast(self.path_gen1, self.gen1["assembly_accession"], cog_db_path)
        cog_file_gen2 = cog_analysis.run_rpsblast(self.path_gen2, self.gen2["assembly_accession"], cog_db_path)

        cog_gen1 = cog_analysis.parse_rpsblast_results(cog_file_gen1, cognames_file, fun_file)
        cog_gen2 = cog_analysis.parse_rpsblast_results(cog_file_gen2, cognames_file, fun_file)

        self.cog_gen1 = cog_gen1
        self.cog_gen2 = cog_gen2

        data_test = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2)
        print("exemples prot blast gen1 :",[d["qseqid"] for d in data_test[:5]])
        print("exemples prot blast gen2 :",[d["sseqid"] for d in data_test[:5]])


        data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2, cog_gen1, cog_gen2)

        cog_colors = self.get_cog_colors()

        print(set(d["cog_category"] for d in data))

        colors = [cog_colors.get(d["cog_category"],"black") for d in data]

        self.display_dotplot_in_interface(data=data,
                                          title = f"COG dotplot: {self.gen1['organism_name']} vs {self.gen2['organism_name']}",
                                          colors = colors,
                                          diagonal_blocks=None,
                                          show_cog_info=True)
        self.show_cog_legend(cog_colors)

        self.status_label.config(text="Dotplot coloured by COG category displayed.")

    def show_cog_legend(self, cog_colors):              # Afficher la légende des couleurs dans un popup
        legend_window = tk.Toplevel(self.fenetre)
        legend_window.title("COG colour code")
        legend_window.geometry("520x700")

        title = tk.Label(legend_window, 
                         text="COG category reference",
                         font=("Arial", 12))
        title.pack(pady=10)

        # Les groupes, les descriptions et les couleurs viennent de cog_style :
        # la legende ne peut donc plus diverger du dotplot ni des camemberts.
        groups = cog_style.COG_GROUPS

        canvas = tk.Canvas(legend_window)
        scrollbar = tk.Scrollbar(legend_window, orient = "vertical", command = canvas.yview)
        scroll_frame = tk.Frame(canvas)

        scroll_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0,0), window=scroll_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side = tk.RIGHT, fill = tk.Y)

        for group_name, categories in groups:
            #titre du groupe
            group_label = tk.Label(scroll_frame, text=group_name, font=("Arial", 10), fg="white", bg="#444444", anchor='w',padx=5)
            group_label.pack(fill=tk.X, padx=10, pady=(8,2))

            for category, description, _ in categories :

                row= tk.Frame(scroll_frame)
                row.pack(fill=tk.X, padx=15, pady=2)

                color_box = tk.Label(row, text=" ", bg=cog_colors.get(category, cog_style.DEFAULT_COLOR), width=3)
                color_box.pack(side=tk.LEFT, padx=5)

                label=tk.Label(row, text=f"{category} : {description}", anchor="w", justify="left")
                label.pack(side=tk.LEFT)
    #  -----------------------
    #       Plus d'Infos 
    #  -----------------------
    def show_info_popup(self):
        if self.blast_result_file is None:
            messagebox.showwarning("Info unavailable", "Please run BLAST first.")
            return
        
        #calcul des stats BLAST
        stats = genome_stats.compute_blast_stats(self.blast_result_file)
        if stats is None :
            messagebox.showwarning("Info unavailable", "No BLAST result found.")
            return
        
        interpretation = genome_stats.get_interpretation(stats["pident_mean"])

        #Création de la fenetre popup
        popup = tk.Toplevel(self.fenetre)
        popup.title("More info")
        popup.geometry("1100x700")

        #Scrollbar
        canvas_scroll = tk.Canvas(popup)
        scrollbar = tk.Scrollbar(popup, orient="vertical", command=canvas_scroll.yview)
        scroll_frame = tk.Frame(canvas_scroll)
        scroll_frame.bind("<Configure>", lambda e: canvas_scroll.configure(scrollregion=canvas_scroll.bbox("all")))
        canvas_scroll.create_window((0,0), window=scroll_frame, anchor='nw')
        canvas_scroll.configure(yscrollcommand=scrollbar.set)
        canvas_scroll.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        #titre
        tk.Label(scroll_frame,
                 text=f"Comparison: {self.gen1['organism_name']} vs {self.gen2['organism_name']}",
                 font = ("Arial", 13)).pack(pady=10)
        #stats Blast
        stats_frame = tk.LabelFrame(scroll_frame, text="BLAST statistics", font=("Arial",10), padx=10, pady=10)
        stats_frame.pack(fill=tk.X, padx=20, pady=5)

        stats_data= [
            ("Mean identity", f"{stats['pident_mean']}%"),
            ("Minimum identity", f"{stats['pident_min']}%"),
            ("Maximum identity", f"{stats['pident_max']}%"),
            ("Total number of alignments", str(stats['nb_hits'])),
            (f"Conserved proteins ({self.gen1['organism_name']})", str(stats['nb_proteins_gen1'])),
            (f"Conserved proteins ({self.gen2['organism_name']})", str(stats['nb_proteins_gen2'])),
        ]

        for label, value in stats_data:
            row=tk.Frame(stats_frame)
            row.pack(fill=tk.X, pady=2)
            tk.Label(row, text=label+" :", width=35, anchor="w", font=("Arial", 10)).pack(side=tk.LEFT)
            tk.Label(row, text = value, anchor='w', font=("Arial", 10)).pack(side=tk.LEFT)

        #interpretations
        interpr_frame= tk.LabelFrame(scroll_frame, text='Interpretation', font=("Arial",11), padx=10, pady=10)
        interpr_frame.pack(fill=tk.X, padx=20, pady=5)
        tk.Label(interpr_frame, text=interpretation, wraplength=820, justify="left", font=("Arial", 10)).pack(anchor='w')

        # camemberts COG
        cog_frame = tk.LabelFrame(scroll_frame, text="COG category distribution", font=("Arial", 11), padx=10, pady=10)
        cog_frame.pack(fill=tk.X, padx=20, pady=20)

        #vérif si l'analyse COG a été faite
        if not hasattr(self, 'cog_gen1') or self.cog_gen1 is None :
            tk.Label(cog_frame,
                     text ="Run the COG colour code first to see the pie charts.",
                     fg="gray", font=("Arial",10, "italic")).pack(pady=20)
        else :
            dist1 = genome_stats.compute_cog_distribution(self.cog_gen1)
            dist2 = genome_stats.compute_cog_distribution(self.cog_gen2)
            blast_data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2)
            dist_common = genome_stats.compute_common_cog(self.cog_gen1, self.cog_gen2, blast_data)

            fig = Figure(figsize=(10, 10), dpi=100)

            # Les trois camemberts sont dessines par la meme fonction, avec
            # le code couleur de cog_style : une categorie a donc la meme
            # couleur dans les trois camemberts, dans le dotplot COG et dans
            # la legende.
            cog_style.draw_pie(fig.add_subplot(2, 2, 1), dist1,
                               self.gen1['organism_name'], short_labels=True)
            cog_style.draw_pie(fig.add_subplot(2, 2, 2), dist2,
                               self.gen2['organism_name'], short_labels=True)
            cog_style.draw_pie(fig.add_subplot(2, 1, 2), dist_common, "Shared proteins")

            fig.tight_layout()

            canvas_fig = FigureCanvasTkAgg(fig, master=cog_frame)
            canvas_fig.draw()
            canvas_fig.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    #  -----------------------
    #     Seuil d'Identité
    #  -----------------------
    def apply_pident_filter(self):
        if self.blast_result_file is None :
            messagebox.showwarning("Filter", "Please run BLAST first.")
            return

        # récupérer le seuil
        try :
            threshold = float(self.pident_filter_entry.get())
            if threshold < 0 or threshold > 100 :
                raise ValueError
        except ValueError :
            messagebox.showwarning("Filter", "Please enter a number between 0 and 100.")
            return

        mode = self.current_plot_mode
        cog_gen1 = getattr(self, 'cog_gen1', None)
        cog_gen2 = getattr(self, 'cog_gen2', None)

        # Charger les données avec COG uniquement si on est en mode COG
        if mode == "cog" and cog_gen1 is not None:
            data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2, cog_gen1, cog_gen2)
        else:
            data = dotplot.get_dotplot_data(self.blast_result_file, self.gff_gen1, self.gff_gen2)

        data_filtered = [d for d in data if float(d["pident"]) >= threshold]

        if len(data_filtered) == 0:
            messagebox.showwarning("Filter", f"No hit with identity >= {threshold}%.")
            return

        title = f"Dotplot (identity >= {threshold}%): {self.gen1['organism_name']} vs {self.gen2['organism_name']}"

        if mode == "cog" and cog_gen1 is not None:
            cog_colors = self.get_cog_colors()
            colors = [cog_colors.get(d["cog_category"], "black") for d in data_filtered]
            self.display_dotplot_in_interface(data=data_filtered, title=title,
                                              colors=colors, diagonal_blocks=None, show_cog_info=True)
        elif mode == "synteny":
            diagonal_blocks = dotplot.detect_diagonal_blocks(data_filtered, bin_size=50000, min_points=30)
            self.display_dotplot_in_interface(data=data_filtered, title=title,
                                              colors='black', diagonal_blocks=diagonal_blocks, show_cog_info=False)
        else:  # mode "simple"
            self.display_dotplot_in_interface(data=data_filtered, title=title,
                                              colors='black', diagonal_blocks=None, show_cog_info=False)

        self.status_label.config(text=f"Filter applied: {len(data_filtered)}/{len(data)} hits shown (>={threshold}%)")

    #  -----------------------
    #       Export
    #  -----------------------
    def make_default_export_name(self, extension):
        name1 = self.gen1["organism_name"].replace(" ", "_").replace("/","_")
        name2 = self.gen2["organism_name"].replace(" ", "_").replace("/","_")
        evalue = self.evalue_entry.get().replace("-","m").replace("+","p").replace(".","_")

        return f"Summary_{name1}_vs_{name2}_{evalue}.{extension}"

    def save_results(self):
        if self.blast_result_file is None : 
            messagebox.showwarning("Export", "Please run BLAST first.")
            return
        
        #popup de choix de format csv ou pdf
        choice_window = tk.Toplevel(self.fenetre)
        choice_window.title("Save results")
        choice_window.geometry("300x150")
        choice_window.resizable(False, False)

        tk.Label(choice_window, 
                 text="Choose the export format:",
                 font=("Arial",11)).pack(pady=15)
        btn_frame = tk.Frame(choice_window)
        btn_frame.pack()
        
        def export_as_csv():            # si on exporte sous csv
            choice_window.destroy()
            default_name = self.make_default_export_name("csv")
            path = fd.asksaveasfilename(defaultextension=".csv", 
                                        initialfile=default_name,
                                        filetypes=[("CSV files", "*.csv")],
                                        title = "Save as CSV")
            if not path :
                return
            cog_gen1 = getattr(self, 'cog_gen1', None)
            cog_gen2 = getattr(self, 'cog_gen2', None)

            self.status_label.config(text ="Generating CSV...")
            self.fenetre.update_idletasks()

            export_results.export_csv(self.blast_result_file, path, cog_gen1, cog_gen2)
            messagebox.showinfo("CSV export", f"CSV file saved:\n{path}")
            self.status_label.config(text=f"CSV exported: {path}")
        
        def export_as_pdf():            #si on exporte sous pdf
            choice_window.destroy()
            default_name = self.make_default_export_name("pdf")
            path = fd.asksaveasfilename(defaultextension=".pdf", 
                                        initialfile=default_name,
                                        filetypes=[("PDF files", "*.pdf")],
                                        title = "Save as PDF")
            if not path :
                return
            
            stats = genome_stats.compute_blast_stats(self.blast_result_file)
            interpretation = genome_stats.get_interpretation(stats["pident_mean"])
            cog_gen1 = getattr(self, 'cog_gen1', None)
            cog_gen2 = getattr(self, 'cog_gen2', None)

            self.status_label.config(text ="Generating PDF...")
            self.fenetre.update_idletasks()

            export_results.export_pdf(output_path=  path,
                                      name1 = self.gen1["organism_name"],
                                      name2 = self.gen2["organism_name"],
                                      stats = stats, 
                                      interpretation = interpretation,
                                      blast_result_file = self.blast_result_file,
                                      gff_gen1 = self.gff_gen1,
                                      gff_gen2 = self.gff_gen2,
                                      cog_gen1 = cog_gen1, 
                                      cog_gen2 = cog_gen2)
            messagebox.showinfo("PDF export", f"PDF file saved:\n{path}")
            self.status_label.config(text=f"PDF exported: {path}")
        tk.Button(btn_frame, text="CSV", width =10, font=("Arial",10), command=export_as_csv).pack(side=tk.LEFT, padx=10)
        tk.Button(btn_frame, text="PDF", width =10, font=("Arial",10), command=export_as_pdf).pack(side=tk.LEFT, padx=10)
       
interface = Genomes()

