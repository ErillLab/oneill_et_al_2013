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
           "Corynebacterium_diphtheriae"]
misc = ["Helicobacter_pylori_HPAG1",
        ]
enteros = ["Enterococcus_faecalis",
           "Bifidobacterium_breve",
           "Bacteroides_thetaiotaomicron",
           "Fusobacterium_nucleatum",
           "Lactobacillus_gasseri",
           "Proteus_mirabilis",
           "Escherichia_coli",
           "Eubacterium_rectale"]

acinetos = ["Acinetobacter_baumannii",
               "Acinetobacter_ADP1",
               "Acinetobacter_calcoaceticus",
               "Acinetobacter_oleivorans"]

pseudos = [org for org in orgs if "Pseudo" in org]
psychros = [org for org in orgs if "Psychro" in org]
all_orgs = list(set(orgs + gammas + firmicutes + actinos + enteros))
new_orgs = gammas + firmicutes
pvp_orgs = pseudos + psychros
last_orgs = pvp_orgs + new_orgs + actinos
ecoli = "Escherichia_coli"
group_names = ["pseudos","psychros","pvp_orgs","gammas",
               "firmicutes","enteros","actinos","last_orgs"]

#orgs used for expression benchmarking

validation_orgs = ['Thermotoga_maritima_MSB8',
                   'Rhodococcus_jostii_RHA1',
                   'Psychrobacter_arcticus_273_4',
                   'Shewanella_oneidensis_MR_1',
                   'Bacillus_subtilis_168',
                   'Deinococcus_radiodurans_R1',
                   'Bacillus_anthracis__Ames_Ancestor',
                   'Streptomyces_coelicolor_A3_2',
                   'Enterococcus_faecalis_V583',
                   'Synechocystis_PCC_6803',
                   'Pseudomonas_fluorescens_Pf_5',
                   'Neisseria_gonorrhoeae_FA_1090',
                   'Lactococcus_lactis_Il1403',
                   'Streptococcus_pneumoniae_R6',
                   'Pseudomonas_aeruginosa_PAO1',
                   'Yersinia_pestis_CO92',
                   'Propionibacterium_freudenreichii_shermanii_CIRM_BIA1',
                   'Vibrio_cholerae_O1_biovar_El_Tor_N16961',
                   'Thermus_thermophilus_HB8',
                   'Mycoplasma_gallisepticum_R_low',
                   'Pseudomonas_putida_KT2440',
                   'Staphylococcus_aureus_COL',
                   'Neisseria_meningitidis_MC58',
                   'Clostridium_perfringens_str_13',
                   'Haemophilus_influenzae_Rd_KW20',
                   'Streptomyces_avermitilis_MA_4680',
                   'Listeria_monocytogenes_EGD_e',
                   'Clostridium_acetobutylicum_ATCC_824',
                   'Ralstonia_solanacearum_GMI1000',
                   'Chlamydophila_pneumoniae_AR39',
                   'Caulobacter_crescentus_CB15',
                   'Myxococcus_xanthus_DK_1622',
                   'Escherichia_coli',
                   'Mycobacterium_smegmatis']



def org_group(org):
    lookup = {"pseudos":pseudos, #
              "psychros":psychros,
              "gammas":gammas, 
              "firmicutes":firmicutes,
              "actinos":actinos,
              "enteros":enteros,
              # "all_orgs":all_orgs
              }
    for group in ["pseudos", "psychros", "gammas", "firmicutes",
                  "actinos", "enteros", "all_orgs"]:
        if org in lookup[group]:
            return group
    assert False, "couldn't find %s in lookup table" % org


del(org) # namespace leak
