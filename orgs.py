"""This file contains lists of organisms"""

orgs = ["Acinetobacter_ADP1",
        "Pseudomonas_fluorescens", 
        "Acinetobacter_baumannii", 
        "Pseudomonas_fulva", 
        "Azotobacter_vinelandii", 
        "Pseudomonas_mendocina", 
        "Cellvibrio__gilvus", 
        "Pseudomonas_putida", 
        "Cellvibrio_japonicus", 
        "Pseudomonas_stutzeri", 
        "Pseudomonas_syringae", 
        "Moraxella_catarrhalis", 
        "Psychrobacter_arcticus_273", 
        "Pseudomonas_aeruginosa",  
        "Psychrobacter_cryohalolentis", 
        "Pseudomonas_brassicacearum", 
        "Psychrobacter_PRwf", 
        "Pseudomonas_entomophila"]

gammas = ["Haemophilus_influenzae", 
          "Escherichia_coli", 
          "Salmonella_typhimurium", 
          "Shewanella_oneidensis", 
          "Vibrio_cholerae"
          ]

firmicutes = ["Bacillus_subtilis", 
              "Lactobacillus_plantarum", 
              "Lactococcus_lactis", 
              "Bacillus_halodurans", 
              "Enterococcus_faecalis", 
              ]
#
actinos = ["Arthrobacter", 
           "Mycobacterium_smegmatis", 
           "Streptomyces_coelicolor", 
           "Bifidobacterium_breve", 
           "Frankia", 
           "Corynebacterium_diphtheriae"]
misc = ["Helicobacter_pylori_HPAG1",
        ]
pseudos = [org for org in orgs if "Pseudo" in org]
psychros = [org for org in orgs if "Psychro" in org]
all_orgs = orgs + gammas + firmicutes + actinos
new_orgs = gammas + firmicutes
pvp_orgs = pseudos + psychros
last_orgs = pvp_orgs + new_orgs + actinos
ecoli = "Escherichia_coli"
del(org) # namespace leak
