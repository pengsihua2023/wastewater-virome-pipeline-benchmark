# CAMISIM Community Simulation: Species List and Abundance

## Data Source

Sample names in this document follow the convention `HTP-YYYY-MM-DD`, where **HTP** stands
for **Hyperion Treatment Plant** (Hyperion Water Reclamation Plant, Los Angeles, California,
USA) and `YYYY-MM-DD` is the sample collection date. This naming convention originates from
the figshare analysis results dataset accompanying the primary publication.

**Primary publication:**
Grimm, S. L. et al. Deep metatranscriptomic sequencing data of wastewater from Los Angeles,
USA, 2023–2024. *Scientific Data* 13, 158 (2026).
https://doi.org/10.1038/s41597-025-06475-7

**Analysis results dataset (source of sample names and Kraken2 classification data):**
Grimm, S. Output of https://github.com/naobservatory/mgs-workflow/tree/2.5.0, used for Deep
wastewater metatranscriptomic sequencing data, Los Angeles, USA, 2023–2024. figshare (2025).
https://doi.org/10.6084/m9.figshare.28454990.v1

The dataset comprises 20 wastewater samples collected weekly at the Hyperion Water Reclamation
Plant (serving ~4 million residents) between December 3, 2023 and April 14, 2024 (Table 1).
Sample **HTP-2024-01-15** (January 15, 2024) was selected for this simulation as it
represents the peak winter 2023–2024 disease season and showed the highest viral signal among
all 20 samples. Taxonomic classification was performed with Kraken2 v2.1.3 against the
Standard database (2024-06-05 build).

| Sample ID | Collection date | Season |
|-----------|----------------|--------|
| HTP-2023-12-03 | 2023-12-03 | Early winter |
| HTP-2023-12-10 | 2023-12-10 | Early winter |
| HTP-2023-12-17 | 2023-12-17 | Early winter |
| HTP-2023-12-25 | 2023-12-25 | Early winter |
| HTP-2024-01-01 | 2024-01-01 | Winter |
| HTP-2024-01-07 | 2024-01-07 | Winter |
| **HTP-2024-01-15** | **2024-01-15** | **Winter (selected)** |
| HTP-2024-01-21 | 2024-01-21 | Winter |
| HTP-2024-01-28 | 2024-01-28 | Winter |
| HTP-2024-02-04 | 2024-02-04 | Winter |
| HTP-2024-02-11 | 2024-02-11 | Winter |
| HTP-2024-02-19 | 2024-02-19 | Winter |
| HTP-2024-02-25 | 2024-02-25 | Late winter |
| HTP-2024-03-03 | 2024-03-03 | Late winter |
| HTP-2024-03-10 | 2024-03-10 | Late winter |
| HTP-2024-03-17 | 2024-03-17 | Spring |
| HTP-2024-03-25 | 2024-03-25 | Spring |
| HTP-2024-03-31 | 2024-03-31 | Spring |
| HTP-2024-04-07 | 2024-04-07 | Spring |
| HTP-2024-04-14 | 2024-04-14 | Spring |

---

# Automatically Selected Bacteria/Archaea Species: Reads and Abundance

A total of 490 Bacteria/Archaea species were selected by the automated pipeline from sample
HTP-2024-01-15 (Hyperion Water Reclamation Plant, Los Angeles, California, USA; 24-h composite influent; rRNA-depleted library; January 15, 2024, peak of 2023–2024 disease season).
Selection criteria: NCBI rank = `species` and clade reads >= 20 in the sample.
Actual reads (clade) from Kraken2 classification are used directly as
simulation design reads; the two columns are identical for all species.

| No. | Name | NCBI taxid | NCBI rank | Actual reads (HTP-2024-01-15, clade) | Simulation design reads | Relative abundance in simulation (%) |
|-----|------|-----------|-----------|--------------------------------------|------------------------|-------------------------------------|
| 1 | Phocaeicola vulgatus | 821 | S | 2124 | 2124 | 3.270723% |
| 2 | Agathobacter rectalis | 39491 | S | 1909 | 1909 | 2.939647% |
| 3 | Uruburuella suis | 252130 | S | 1855 | 1855 | 2.856493% |
| 4 | Segatella copri | 165179 | S | 1669 | 1669 | 2.570074% |
| 5 | Aliarcobacter cryaerophilus | 28198 | S | 1427 | 1427 | 2.197421% |
| 6 | Sebaldella termitidis | 826 | S | 1202 | 1202 | 1.850946% |
| 7 | Aquaspirillum sp. LM1 | 1938604 | S | 949 | 949 | 1.461354% |
| 8 | uncultured Acetobacteroides sp. | 1760811 | S | 947 | 947 | 1.458274% |
| 9 | Cetobacterium somerae | 188913 | S | 940 | 940 | 1.447495% |
| 10 | Fusobacterium varium | 856 | S | 815 | 815 | 1.255009% |
| 11 | Sneathia vaginalis | 187101 | S | 791 | 791 | 1.218052% |
| 12 | Simiaoa sunii | 2763672 | S | 787 | 787 | 1.211892% |
| 13 | Vitreoscilla massiliensis | 1689272 | S | 780 | 780 | 1.201113% |
| 14 | Fusobacterium canifelinum | 285729 | S | 776 | 776 | 1.194953% |
| 15 | uncultured Propionivibrio sp. | 426737 | S | 763 | 763 | 1.174935% |
| 16 | Haliovirga abyssi | 2996794 | S | 739 | 739 | 1.137977% |
| 17 | [Ruminococcus] torques | 33039 | S | 685 | 685 | 1.054824% |
| 18 | Blautia obeum | 40520 | S | 641 | 641 | 0.987068% |
| 19 | uncultured Ilyobacter sp. | 544433 | S | 611 | 611 | 0.940872% |
| 20 | Bacteroides graminisolvens | 477666 | S | 520 | 520 | 0.800742% |
| 21 | Psychrobacter sp. M13 | 3067275 | S | 464 | 464 | 0.714508% |
| 22 | Streptococcus parasuis | 1501662 | S | 450 | 450 | 0.692950% |
| 23 | Streptobacillus moniliformis | 34105 | S | 444 | 444 | 0.683710% |
| 24 | Bacteroides stercoris | 46506 | S | 437 | 437 | 0.672931% |
| 25 | Faecalibacterium prausnitzii | 853 | S | 428 | 428 | 0.659072% |
| 26 | Staphylococcus epidermidis | 1282 | S | 418 | 418 | 0.643673% |
| 27 | Fusobacterium necrophorum | 859 | S | 394 | 394 | 0.606716% |
| 28 | uncultured Macellibacteroides sp. | 1243007 | S | 383 | 383 | 0.589777% |
| 29 | Bacteroides uniformis | 820 | S | 370 | 370 | 0.569759% |
| 30 | Halomonas sp. JS92-SW72 | 2306583 | S | 356 | 356 | 0.548200% |
| 31 | Bacillus velezensis | 492670 | S | 354 | 354 | 0.545120% |
| 32 | Ruminococcus bicirculans (ex Wegman et al. 2014) | 1160721 | S | 342 | 342 | 0.526642% |
| 33 | Phocaeicola coprophilus | 387090 | S | 321 | 321 | 0.494304% |
| 34 | uncultured Tolumonas sp. | 263765 | S | 321 | 321 | 0.494304% |
| 35 | Bacteriovorax stolpii | 960 | S | 315 | 315 | 0.485065% |
| 36 | Phocaeicola dorei | 357276 | S | 315 | 315 | 0.485065% |
| 37 | Flavobacterium sp. K5-23 | 2746225 | S | 314 | 314 | 0.483525% |
| 38 | Wujia chipingensis | 2763670 | S | 301 | 301 | 0.463506% |
| 39 | Laribacter hongkongensis | 168471 | S | 295 | 295 | 0.454267% |
| 40 | Paralysiella testudinis | 2809020 | S | 289 | 289 | 0.445028% |
| 41 | Ruminococcus champanellensis | 1161942 | S | 278 | 278 | 0.428089% |
| 42 | Pseudomonas sp. CIP-10 | 2892442 | S | 277 | 277 | 0.426549% |
| 43 | Lachnospira eligens | 39485 | S | 276 | 276 | 0.425009% |
| 44 | Tolumonas auensis | 43948 | S | 274 | 274 | 0.421929% |
| 45 | Coprococcus eutactus | 33043 | S | 267 | 267 | 0.411150% |
| 46 | Cloacibacterium caeni | 2004710 | S | 254 | 254 | 0.391132% |
| 47 | Azospira inquinata | 2785627 | S | 247 | 247 | 0.380352% |
| 48 | Fusobacterium gonidiaformans | 849 | S | 243 | 243 | 0.374193% |
| 49 | Thauera aminoaromatica | 164330 | S | 237 | 237 | 0.364954% |
| 50 | Roseburia sp. 831b | 1261635 | S | 236 | 236 | 0.363414% |
| 51 | Flavobacterium channae | 2897181 | S | 233 | 233 | 0.358794% |
| 52 | Stenotrophomonas maltophilia | 40324 | S | 232 | 232 | 0.357254% |
| 53 | Flavobacterium nackdongense | 2547394 | S | 217 | 217 | 0.334156% |
| 54 | Pseudomonas simiae | 321846 | S | 217 | 217 | 0.334156% |
| 55 | Anaerostipes hadrus | 649756 | S | 215 | 215 | 0.331076% |
| 56 | Mediterraneibacter gnavus | 33038 | S | 212 | 212 | 0.326456% |
| 57 | Neisseria zalophi | 640030 | S | 209 | 209 | 0.321837% |
| 58 | Pseudoleptotrichia goodfellowii | 157692 | S | 207 | 207 | 0.318757% |
| 59 | Mycoplasma miroungirhinis | 754516 | S | 206 | 206 | 0.317217% |
| 60 | Roseburia hominis | 301301 | S | 201 | 201 | 0.309518% |
| 61 | Thauera sp. K11 | 2005884 | S | 200 | 200 | 0.307978% |
| 62 | Duodenibacillus massiliensis | 1852381 | S | 194 | 194 | 0.298738% |
| 63 | Acetoanaerobium sticklandii | 1511 | S | 191 | 191 | 0.294119% |
| 64 | Eubacterium ventriosum | 39496 | S | 190 | 190 | 0.292579% |
| 65 | Fusobacterium hwasookii | 1583098 | S | 187 | 187 | 0.287959% |
| 66 | Staphylococcus haemolyticus | 1283 | S | 182 | 182 | 0.280260% |
| 67 | Pseudomonas aeruginosa | 287 | S | 179 | 179 | 0.275640% |
| 68 | Pyramidobacter piscolens | 638849 | S | 178 | 178 | 0.274100% |
| 69 | Candidatus Saccharimonas aalborgensis | 1332188 | S | 176 | 176 | 0.271020% |
| 70 | Coprococcus sp. ART55/1 | 751585 | S | 166 | 166 | 0.255621% |
| 71 | Roseburia intestinalis | 166486 | S | 164 | 164 | 0.252542% |
| 72 | Vescimonas fastidiosa | 2714353 | S | 162 | 162 | 0.249462% |
| 73 | Moraxella osloensis | 34062 | S | 161 | 161 | 0.247922% |
| 74 | Azonexus hydrophilus | 418702 | S | 160 | 160 | 0.246382% |
| 75 | Selenomonas ruminantium | 971 | S | 158 | 158 | 0.243302% |
| 76 | Fibrobacter succinogenes | 833 | S | 156 | 156 | 0.240223% |
| 77 | Escherichia coli | 562 | S | 154 | 154 | 0.237143% |
| 78 | Leptotrichia sp. oral taxon 212 | 712357 | S | 153 | 153 | 0.235603% |
| 79 | Candidatus Megaera polyxenophila | 988779 | S | 148 | 148 | 0.227903% |
| 80 | Anaerolinea thermophila | 167964 | S | 147 | 147 | 0.226364% |
| 81 | Blautia wexlerae | 418240 | S | 146 | 146 | 0.224824% |
| 82 | Epilithonimonas vandammei | 2487072 | S | 144 | 144 | 0.221744% |
| 83 | Paludibacter propionicigenes | 185300 | S | 144 | 144 | 0.221744% |
| 84 | Parabacteroides distasonis | 823 | S | 141 | 141 | 0.217124% |
| 85 | Sutterella wadsworthensis | 40545 | S | 136 | 136 | 0.209425% |
| 86 | Gracilinema caldarium | 215591 | S | 134 | 134 | 0.206345% |
| 87 | Candidatus Mycosynbacter amalyticus | 2665156 | S | 133 | 133 | 0.204805% |
| 88 | Peredibacter starrii | 28202 | S | 133 | 133 | 0.204805% |
| 89 | Fusobacterium nucleatum | 851 | S | 129 | 129 | 0.198646% |
| 90 | Flavobacterium faecale | 1355330 | S | 129 | 129 | 0.198646% |
| 91 | uncultured Desulfobulbus sp. | 239745 | S | 126 | 126 | 0.194026% |
| 92 | Klebsiella pneumoniae | 573 | S | 125 | 125 | 0.192486% |
| 93 | Quatrionicoccus australiensis | 138118 | S | 124 | 124 | 0.190946% |
| 94 | Pelolinea submarina | 913107 | S | 120 | 120 | 0.184787% |
| 95 | Sphaerotilus microaerophilus | 2914710 | S | 119 | 119 | 0.183247% |
| 96 | Aquirhabdus parva | 2283318 | S | 119 | 119 | 0.183247% |
| 97 | Haliangium ochraceum | 80816 | S | 119 | 119 | 0.183247% |
| 98 | Brevefilum fermentans | 1986204 | S | 117 | 117 | 0.180167% |
| 99 | Brucepastera parasyntrophica | 2880008 | S | 116 | 116 | 0.178627% |
| 100 | Salmonella enterica | 28901 | S | 115 | 115 | 0.177087% |
| 101 | Selenobaculum gbiensis | 3054208 | S | 115 | 115 | 0.177087% |
| 102 | Anaerotruncus colihominis | 169435 | S | 115 | 115 | 0.177087% |
| 103 | Ferribacterium limneticum | 76259 | S | 114 | 114 | 0.175547% |
| 104 | Azospira restricta | 404405 | S | 113 | 113 | 0.174007% |
| 105 | Solibaculum mannosilyticum | 2780922 | S | 112 | 112 | 0.172467% |
| 106 | Flavonifractor plautii | 292800 | S | 112 | 112 | 0.172467% |
| 107 | uncultured Bacteroides sp. | 162156 | S | 110 | 110 | 0.169388% |
| 108 | Neisseria bacilliformis | 267212 | S | 110 | 110 | 0.169388% |
| 109 | Conchiformibius steedae | 153493 | S | 110 | 110 | 0.169388% |
| 110 | Prevotella herbatica | 2801997 | S | 110 | 110 | 0.169388% |
| 111 | Dechloromonas sp. A34 | 447588 | S | 109 | 109 | 0.167848% |
| 112 | Odoribacter splanchnicus | 28118 | S | 109 | 109 | 0.167848% |
| 113 | Arcobacter acticola | 1849015 | S | 108 | 108 | 0.166308% |
| 114 | Bacteroides eggerthii | 28111 | S | 106 | 106 | 0.163228% |
| 115 | Dechloromonas denitrificans | 281362 | S | 106 | 106 | 0.163228% |
| 116 | Leptotrichia sp. oral taxon 221 | 712362 | S | 103 | 103 | 0.158609% |
| 117 | Paludicola sp. MB14-C6 | 3070656 | S | 103 | 103 | 0.158609% |
| 118 | Solimonas sp. K1W22B-7 | 2303331 | S | 102 | 102 | 0.157069% |
| 119 | Blautia luti | 89014 | S | 102 | 102 | 0.157069% |
| 120 | Neisseria brasiliensis | 2666100 | S | 102 | 102 | 0.157069% |
| 121 | Anaerobutyricum hallii | 39488 | S | 101 | 101 | 0.155529% |
| 122 | Sterolibacterium denitrificans | 157592 | S | 101 | 101 | 0.155529% |
| 123 | Sphaerotilus sulfidivorans | 639200 | S | 101 | 101 | 0.155529% |
| 124 | Streptomyces sp. T12 | 477697 | S | 101 | 101 | 0.155529% |
| 125 | Psychrobacter sp. van23A | 3064892 | S | 101 | 101 | 0.155529% |
| 126 | Streptococcus oriscaviae | 2781599 | S | 100 | 100 | 0.153989% |
| 127 | Rheinheimera sp. F8 | 1763998 | S | 99 | 99 | 0.152449% |
| 128 | Polyangium aurulentum | 2567896 | S | 99 | 99 | 0.152449% |
| 129 | Pseudomonas putida | 303 | S | 99 | 99 | 0.152449% |
| 130 | Ruminococcus bovis | 2564099 | S | 98 | 98 | 0.150909% |
| 131 | Arcobacter sp. FWKO B | 2593672 | S | 98 | 98 | 0.150909% |
| 132 | Blautia producta | 33035 | S | 97 | 97 | 0.149369% |
| 133 | Aquella oligotrophica | 2067065 | S | 97 | 97 | 0.149369% |
| 134 | Pusillibacter faecalis | 2714358 | S | 96 | 96 | 0.147829% |
| 135 | Blautia liquoris | 2779518 | S | 96 | 96 | 0.147829% |
| 136 | Paraprevotella clara | 454154 | S | 96 | 96 | 0.147829% |
| 137 | Parabacteroides merdae | 46503 | S | 96 | 96 | 0.147829% |
| 138 | Lawsonibacter asaccharolyticus | 2108523 | S | 94 | 94 | 0.144750% |
| 139 | Arcobacter ellisii | 913109 | S | 94 | 94 | 0.144750% |
| 140 | Sulfuritortus calidifontis | 1914471 | S | 93 | 93 | 0.143210% |
| 141 | Aromatoleum bremense | 76115 | S | 93 | 93 | 0.143210% |
| 142 | Chondromyces crocatus | 52 | S | 91 | 91 | 0.140130% |
| 143 | Ethanoligenens harbinense | 253239 | S | 91 | 91 | 0.140130% |
| 144 | Aeromonas salmonicida | 645 | S | 91 | 91 | 0.140130% |
| 145 | Desulfobulbus propionicus | 894 | S | 91 | 91 | 0.140130% |
| 146 | Sutterella faecalis | 2584944 | S | 90 | 90 | 0.138590% |
| 147 | Pseudobutyrivibrio xylanivorans | 185007 | S | 90 | 90 | 0.138590% |
| 148 | Dorea longicatena | 88431 | S | 89 | 89 | 0.137050% |
| 149 | Blautia sp. SC05B48 | 2479767 | S | 87 | 87 | 0.133970% |
| 150 | Vescimonas coprocola | 2714355 | S | 86 | 86 | 0.132430% |
| 151 | Blautia hansenii | 1322 | S | 86 | 86 | 0.132430% |
| 152 | Riemerella anatipestifer | 34085 | S | 84 | 84 | 0.129351% |
| 153 | Azoarcus sp. KH32C | 748247 | S | 83 | 83 | 0.127811% |
| 154 | Acinetobacter baumannii | 470 | S | 83 | 83 | 0.127811% |
| 155 | Ruminococcus albus | 1264 | S | 83 | 83 | 0.127811% |
| 156 | Intestinibacillus sp. NTUH-41-i26 | 3079303 | S | 82 | 82 | 0.126271% |
| 157 | Flavobacterium cyclinae | 2895947 | S | 81 | 81 | 0.124731% |
| 158 | Moraxella canis | 90239 | S | 80 | 80 | 0.123191% |
| 159 | uncultured Carboxylicivirga sp. | 1628156 | S | 80 | 80 | 0.123191% |
| 160 | Moraxella bovoculi | 386891 | S | 80 | 80 | 0.123191% |
| 161 | Haliscomenobacter hydrossis | 2350 | S | 79 | 79 | 0.121651% |
| 162 | Eikenella exigua | 2528037 | S | 78 | 78 | 0.120111% |
| 163 | Coprococcus catus | 116085 | S | 78 | 78 | 0.120111% |
| 164 | Bacillus cereus | 1396 | S | 78 | 78 | 0.120111% |
| 165 | Faecalibacillus intestinalis | 1982626 | S | 77 | 77 | 0.118571% |
| 166 | Ruthenibacterium lactatiformans | 1550024 | S | 76 | 76 | 0.117032% |
| 167 | Oscillibacter hominis | 2763056 | S | 75 | 75 | 0.115492% |
| 168 | Pseudomonas fragi | 296 | S | 75 | 75 | 0.115492% |
| 169 | Arcobacter suis | 1278212 | S | 75 | 75 | 0.115492% |
| 170 | Xanthomonas albilineans | 29447 | S | 74 | 74 | 0.113952% |
| 171 | Verrucomicrobium spinosum | 2736 | S | 72 | 72 | 0.110872% |
| 172 | uncultured Anaeromusa sp. | 673273 | S | 71 | 71 | 0.109332% |
| 173 | Tenuifilum thalassicum | 2590900 | S | 71 | 71 | 0.109332% |
| 174 | Cytophaga hutchinsonii | 985 | S | 71 | 71 | 0.109332% |
| 175 | Capnocytophaga haemolytica | 45243 | S | 70 | 70 | 0.107792% |
| 176 | Oryzomicrobium terrae | 1735038 | S | 70 | 70 | 0.107792% |
| 177 | Chitinimonas arctica | 2594795 | S | 70 | 70 | 0.107792% |
| 178 | Fluviicola taffensis | 191579 | S | 70 | 70 | 0.107792% |
| 179 | Aurantibacillus circumpalustris | 3036359 | S | 69 | 69 | 0.106252% |
| 180 | Streptococcus iniae | 1346 | S | 69 | 69 | 0.106252% |
| 181 | Mariniblastus fucicola | 980251 | S | 69 | 69 | 0.106252% |
| 182 | Bacteroides fragilis | 817 | S | 69 | 69 | 0.106252% |
| 183 | Pseudobdellovibrio exovorus | 453816 | S | 68 | 68 | 0.104712% |
| 184 | Acinetobacter lanii | 2715163 | S | 66 | 66 | 0.101633% |
| 185 | Anaerocolumna cellulosilytica | 433286 | S | 66 | 66 | 0.101633% |
| 186 | Phascolarctobacterium faecium | 33025 | S | 66 | 66 | 0.101633% |
| 187 | Methyloversatilis sp. RAC08 | 1842540 | S | 64 | 64 | 0.098553% |
| 188 | Alkalitalea saponilacus | 889453 | S | 64 | 64 | 0.098553% |
| 189 | Turicimonas muris | 1796652 | S | 64 | 64 | 0.098553% |
| 190 | Nannocystis poenicansa | 2995304 | S | 63 | 63 | 0.097013% |
| 191 | Phototrophicus methaneseepsis | 2710758 | S | 63 | 63 | 0.097013% |
| 192 | Xiashengella succiniciproducens | 2949635 | S | 62 | 62 | 0.095473% |
| 193 | Pseudomonas syringae | 317 | S | 62 | 62 | 0.095473% |
| 194 | Marvinbryantia formatexigens | 168384 | S | 62 | 62 | 0.095473% |
| 195 | Romboutsia ilealis | 1115758 | S | 61 | 61 | 0.093933% |
| 196 | Streptococcus ferus | 1345 | S | 61 | 61 | 0.093933% |
| 197 | Dorea formicigenerans | 39486 | S | 61 | 61 | 0.093933% |
| 198 | Alistipes onderdonkii | 328813 | S | 61 | 61 | 0.093933% |
| 199 | Turneriella parva | 29510 | S | 61 | 61 | 0.093933% |
| 200 | Sutterella megalosphaeroides | 2494234 | S | 60 | 60 | 0.092393% |
| 201 | Aquipluma nitroreducens | 2010828 | S | 59 | 59 | 0.090853% |
| 202 | Intestinimonas butyriciproducens | 1297617 | S | 58 | 58 | 0.089314% |
| 203 | Oxalobacter formigenes | 847 | S | 58 | 58 | 0.089314% |
| 204 | Desulfovibrio desulfuricans | 876 | S | 57 | 57 | 0.087774% |
| 205 | Candidatus Viadribacter manganicus | 1759059 | S | 56 | 56 | 0.086234% |
| 206 | Alistipes communis | 2585118 | S | 56 | 56 | 0.086234% |
| 207 | Ferriphaselus amnicola | 1188319 | S | 56 | 56 | 0.086234% |
| 208 | Bdellovibrio bacteriovorus | 959 | S | 56 | 56 | 0.086234% |
| 209 | Pseudothauera hydrothermalis | 2184083 | S | 55 | 55 | 0.084694% |
| 210 | Segatella hominis | 2518605 | S | 55 | 55 | 0.084694% |
| 211 | Chitinimonas koreensis | 356302 | S | 55 | 55 | 0.084694% |
| 212 | [Ruminococcus] lactaris | 46228 | S | 55 | 55 | 0.084694% |
| 213 | Butyrivibrio crossotus | 45851 | S | 53 | 53 | 0.081614% |
| 214 | Coprococcus comes | 410072 | S | 53 | 53 | 0.081614% |
| 215 | Diaphorobacter limosus | 3036128 | S | 53 | 53 | 0.081614% |
| 216 | Sulfurospirillum cavolei | 366522 | S | 53 | 53 | 0.081614% |
| 217 | Fusobacterium mortiferum | 850 | S | 53 | 53 | 0.081614% |
| 218 | Succinivibrio dextrinosolvens | 83771 | S | 52 | 52 | 0.080074% |
| 219 | Clostridium sp. M62/1 | 411486 | S | 52 | 52 | 0.080074% |
| 220 | Rugosibacter aromaticivorans | 1565605 | S | 52 | 52 | 0.080074% |
| 221 | Pelobacter propionicus | 29543 | S | 52 | 52 | 0.080074% |
| 222 | Sulfuritalea hydrogenivorans | 748811 | S | 52 | 52 | 0.080074% |
| 223 | Candidatus Sulfidibacterium hydrothermale | 2875962 | S | 51 | 51 | 0.078534% |
| 224 | Moraxella ovis | 29433 | S | 51 | 51 | 0.078534% |
| 225 | Clostridium sp. SY8519 | 1042156 | S | 51 | 51 | 0.078534% |
| 226 | Cellvibrio sp. PSBB006 | 1987723 | S | 51 | 51 | 0.078534% |
| 227 | Leptotrichia wadei | 157687 | S | 51 | 51 | 0.078534% |
| 228 | Candidatus Symbiobacter mobilis | 1436290 | S | 51 | 51 | 0.078534% |
| 229 | Desulfomicrobium orale | 132132 | S | 51 | 51 | 0.078534% |
| 230 | Phocaeicola salanitronis | 376805 | S | 51 | 51 | 0.078534% |
| 231 | Flintibacter sp. KGMB00164 | 2610895 | S | 50 | 50 | 0.076994% |
| 232 | Pseudoprevotella muciniphila | 2133944 | S | 50 | 50 | 0.076994% |
| 233 | Aeromonas veronii | 654 | S | 50 | 50 | 0.076994% |
| 234 | Bacteroides caccae | 47678 | S | 50 | 50 | 0.076994% |
| 235 | Alistipes shahii | 328814 | S | 50 | 50 | 0.076994% |
| 236 | Neisseria sp. KEM232 | 655307 | S | 49 | 49 | 0.075455% |
| 237 | Niveibacterium microcysteis | 2811415 | S | 49 | 49 | 0.075455% |
| 238 | Treponema pallidum | 160 | S | 49 | 49 | 0.075455% |
| 239 | Thermocaproicibacter melissae | 2966552 | S | 49 | 49 | 0.075455% |
| 240 | Sediminibacterium sp. TEGAF015 | 575378 | S | 48 | 48 | 0.073915% |
| 241 | Bacteroides ovatus | 28116 | S | 48 | 48 | 0.073915% |
| 242 | Roseimicrobium sp. ORNL1 | 2711231 | S | 47 | 47 | 0.072375% |
| 243 | Bacteroides cellulosilyticus | 246787 | S | 47 | 47 | 0.072375% |
| 244 | Roseburia rectibacter | 2763062 | S | 47 | 47 | 0.072375% |
| 245 | Streptococcus pantholopis | 1811193 | S | 47 | 47 | 0.072375% |
| 246 | Leptotrichia sp. oral taxon 218 | 712361 | S | 47 | 47 | 0.072375% |
| 247 | Treponema primitia | 88058 | S | 47 | 47 | 0.072375% |
| 248 | Moraxella nasovis | 2904121 | S | 46 | 46 | 0.070835% |
| 249 | Aquabacterium olei | 1296669 | S | 46 | 46 | 0.070835% |
| 250 | Collinsella aerofaciens | 74426 | S | 46 | 46 | 0.070835% |
| 251 | Streptococcus dysgalactiae | 1334 | S | 45 | 45 | 0.069295% |
| 252 | Pararhizobium sp. BT-229 | 2986923 | S | 45 | 45 | 0.069295% |
| 253 | Ferrigenium kumadai | 1682490 | S | 45 | 45 | 0.069295% |
| 254 | Acutalibacter muris | 1796620 | S | 45 | 45 | 0.069295% |
| 255 | Bacteroides thetaiotaomicron | 818 | S | 45 | 45 | 0.069295% |
| 256 | Streptococcus anginosus | 1328 | S | 45 | 45 | 0.069295% |
| 257 | Alistipes indistinctus | 626932 | S | 45 | 45 | 0.069295% |
| 258 | Caldilinea aerophila | 133453 | S | 45 | 45 | 0.069295% |
| 259 | Micropruina glycogenica | 75385 | S | 44 | 44 | 0.067755% |
| 260 | Lacrimispora xylanolytica | 29375 | S | 44 | 44 | 0.067755% |
| 261 | Fastidiosipila sanguinis | 236753 | S | 44 | 44 | 0.067755% |
| 262 | Streptococcus thermophilus | 1308 | S | 44 | 44 | 0.067755% |
| 263 | Candidatus Chromulinivorax destructor | 2066483 | S | 43 | 43 | 0.066215% |
| 264 | Leptolyngbya sp. 7M | 2812896 | S | 43 | 43 | 0.066215% |
| 265 | Vandammella animalimorsus | 2029117 | S | 43 | 43 | 0.066215% |
| 266 | Dysgonomonas capnocytophagoides | 45254 | S | 43 | 43 | 0.066215% |
| 267 | Anaeromicropila herbilytica | 2785025 | S | 43 | 43 | 0.066215% |
| 268 | Pseudomonas lundensis | 86185 | S | 42 | 42 | 0.064675% |
| 269 | Acinetobacter sp. Marseille-Q1620 | 2714110 | S | 42 | 42 | 0.064675% |
| 270 | Agathobaculum sp. NTUH-O15-33 | 3079302 | S | 42 | 42 | 0.064675% |
| 271 | Leptospira interrogans | 173 | S | 41 | 41 | 0.063135% |
| 272 | Aquincola tertiaricarbonis | 391953 | S | 40 | 40 | 0.061596% |
| 273 | uncultured Alistipes sp. | 538949 | S | 40 | 40 | 0.061596% |
| 274 | Cellvibrio sp. KY-GH-1 | 2303332 | S | 40 | 40 | 0.061596% |
| 275 | Pseudomonas alcaligenes | 43263 | S | 40 | 40 | 0.061596% |
| 276 | Clostridium perfringens | 1502 | S | 40 | 40 | 0.061596% |
| 277 | Micavibrio aeruginosavorus | 349221 | S | 40 | 40 | 0.061596% |
| 278 | Simonsiella muelleri | 72 | S | 39 | 39 | 0.060056% |
| 279 | Chryseobacterium sp. POL2 | 2713414 | S | 39 | 39 | 0.060056% |
| 280 | Neisseria shayeganii | 607712 | S | 39 | 39 | 0.060056% |
| 281 | Ottowia oryzae | 2109914 | S | 39 | 39 | 0.060056% |
| 282 | Enterobacter sp. T2 | 2707174 | S | 39 | 39 | 0.060056% |
| 283 | Paracidovorax avenae | 80867 | S | 39 | 39 | 0.060056% |
| 284 | Treponema brennaborense | 81028 | S | 39 | 39 | 0.060056% |
| 285 | uncultured Vibrio sp. | 114054 | S | 38 | 38 | 0.058516% |
| 286 | Caproicibacterium lactatifermentans | 2666138 | S | 38 | 38 | 0.058516% |
| 287 | Tepidimonas taiwanensis | 307486 | S | 38 | 38 | 0.058516% |
| 288 | Arcobacter nitrofigilis | 28199 | S | 38 | 38 | 0.058516% |
| 289 | Megamonas funiformis | 437897 | S | 37 | 37 | 0.056976% |
| 290 | Selenomonas sputigena | 69823 | S | 37 | 37 | 0.056976% |
| 291 | Candidatus Paracaedibacter acanthamoebae | 91604 | S | 37 | 37 | 0.056976% |
| 292 | Desulfobulbus oralis | 1986146 | S | 37 | 37 | 0.056976% |
| 293 | Leptotrichia trevisanii | 109328 | S | 37 | 37 | 0.056976% |
| 294 | Niveibacterium sp. SC-1 | 3135646 | S | 37 | 37 | 0.056976% |
| 295 | Tessaracoccus flavescens | 399497 | S | 36 | 36 | 0.055436% |
| 296 | Steroidobacter denitrificans | 465721 | S | 36 | 36 | 0.055436% |
| 297 | Pseudazoarcus pumilus | 2067960 | S | 36 | 36 | 0.055436% |
| 298 | Magnetospirillum gryphiswaldense | 55518 | S | 36 | 36 | 0.055436% |
| 299 | Roseateles amylovorans | 2978473 | S | 36 | 36 | 0.055436% |
| 300 | Limnohabitans sp. 103DPR2 | 1678129 | S | 36 | 36 | 0.055436% |
| 301 | Neisseria lisongii | 2912188 | S | 36 | 36 | 0.055436% |
| 302 | Serratia symbiotica | 138074 | S | 35 | 35 | 0.053896% |
| 303 | Pseudomonas sp. S150 | 2559074 | S | 35 | 35 | 0.053896% |
| 304 | Pseudomonas guguanensis | 1198456 | S | 35 | 35 | 0.053896% |
| 305 | Terricaulis silvestris | 2686094 | S | 35 | 35 | 0.053896% |
| 306 | Labilithrix luteola | 1391654 | S | 35 | 35 | 0.053896% |
| 307 | Qiania dongpingensis | 2763669 | S | 35 | 35 | 0.053896% |
| 308 | Ottowia sp. oral taxon 894 | 1658672 | S | 35 | 35 | 0.053896% |
| 309 | Kineothrix sp. MB12-C1 | 3070215 | S | 35 | 35 | 0.053896% |
| 310 | Treponema phagedenis | 162 | S | 35 | 35 | 0.053896% |
| 311 | Bifidobacterium longum | 216816 | S | 35 | 35 | 0.053896% |
| 312 | Buchnera aphidicola | 9 | S | 35 | 35 | 0.053896% |
| 313 | Bacteroides helcogenes | 290053 | S | 35 | 35 | 0.053896% |
| 314 | Sinimarinibacterium sp. NLF-5-8 | 2698684 | S | 34 | 34 | 0.052356% |
| 315 | Propioniciclava coleopterorum | 2714937 | S | 34 | 34 | 0.052356% |
| 316 | Prevotella sp. E9-3 | 2913621 | S | 34 | 34 | 0.052356% |
| 317 | Sodaliphilus pleomorphus | 2606626 | S | 34 | 34 | 0.052356% |
| 318 | Moraxella haemolytica | 2904119 | S | 34 | 34 | 0.052356% |
| 319 | Pseudomonas fluorescens | 294 | S | 34 | 34 | 0.052356% |
| 320 | Aliarcobacter thereius | 544718 | S | 34 | 34 | 0.052356% |
| 321 | Anaerostipes rhamnosivorans | 1229621 | S | 33 | 33 | 0.050816% |
| 322 | Pseudomonas hefeiensis | 2738125 | S | 33 | 33 | 0.050816% |
| 323 | Oceanirhabdus sp. W0125-5 | 2999116 | S | 33 | 33 | 0.050816% |
| 324 | Cloacibacillus porcorum | 1197717 | S | 33 | 33 | 0.050816% |
| 325 | Dialister succinatiphilus | 487173 | S | 33 | 33 | 0.050816% |
| 326 | [Clostridium] scindens | 29347 | S | 33 | 33 | 0.050816% |
| 327 | Cloacibacterium sp. TD35 | 2976818 | S | 33 | 33 | 0.050816% |
| 328 | Oxalobacter vibrioformis | 933080 | S | 33 | 33 | 0.050816% |
| 329 | Melaminivora suipulveris | 2109913 | S | 33 | 33 | 0.050816% |
| 330 | Leptospira biflexa | 172 | S | 33 | 33 | 0.050816% |
| 331 | Prevotella dentalis | 52227 | S | 33 | 33 | 0.050816% |
| 332 | Roseateles sp. SL47 | 2995138 | S | 32 | 32 | 0.049276% |
| 333 | Changpingibacter yushuensis | 2758440 | S | 32 | 32 | 0.049276% |
| 334 | Denitratisoma sp. DHT3 | 1981880 | S | 32 | 32 | 0.049276% |
| 335 | Methylomusa anaerophila | 1930071 | S | 32 | 32 | 0.049276% |
| 336 | Blautia massiliensis (ex Durand et al. 2017) | 1737424 | S | 32 | 32 | 0.049276% |
| 337 | Thermosulfidibacter takaii | 412593 | S | 32 | 32 | 0.049276% |
| 338 | Comamonas terrigena | 32013 | S | 32 | 32 | 0.049276% |
| 339 | Candidatus Cytomitobacter indipagum | 2601575 | S | 31 | 31 | 0.047737% |
| 340 | Paucibacter aquatile | 2070761 | S | 31 | 31 | 0.047737% |
| 341 | Granulicella sp. WH15 | 2602070 | S | 31 | 31 | 0.047737% |
| 342 | Candidatus Babela massiliensis | 673862 | S | 31 | 31 | 0.047737% |
| 343 | Megasphaera massiliensis | 1232428 | S | 31 | 31 | 0.047737% |
| 344 | uncultured Holophaga sp. | 188475 | S | 31 | 31 | 0.047737% |
| 345 | Acidovorax sp. HDW3 | 2714923 | S | 31 | 31 | 0.047737% |
| 346 | Microvirgula aerodenitrificans | 57480 | S | 31 | 31 | 0.047737% |
| 347 | Acinetobacter johnsonii | 40214 | S | 31 | 31 | 0.047737% |
| 348 | Legionella israelensis | 454 | S | 30 | 30 | 0.046197% |
| 349 | uncultured Hyphomonas sp. | 225298 | S | 30 | 30 | 0.046197% |
| 350 | Chryseolinea soli | 2321403 | S | 30 | 30 | 0.046197% |
| 351 | Christensenella minuta | 626937 | S | 30 | 30 | 0.046197% |
| 352 | Denitratisoma oestradiolicum | 311182 | S | 30 | 30 | 0.046197% |
| 353 | Haemophilus parainfluenzae | 729 | S | 30 | 30 | 0.046197% |
| 354 | Arcobacter venerupis | 1054033 | S | 30 | 30 | 0.046197% |
| 355 | Prevotella sp. Rep29 | 2691580 | S | 30 | 30 | 0.046197% |
| 356 | Morganella morganii | 582 | S | 30 | 30 | 0.046197% |
| 357 | Sulfurospirillum deleyianum | 65553 | S | 30 | 30 | 0.046197% |
| 358 | Nitratidesulfovibrio vulgaris | 881 | S | 29 | 29 | 0.044657% |
| 359 | Bifidobacterium adolescentis | 1680 | S | 29 | 29 | 0.044657% |
| 360 | Staphylococcus aureus | 1280 | S | 29 | 29 | 0.044657% |
| 361 | Thauera humireducens | 1134435 | S | 29 | 29 | 0.044657% |
| 362 | Breznakiella homolactica | 2798577 | S | 29 | 29 | 0.044657% |
| 363 | Maliibacterium massiliense | 2086585 | S | 29 | 29 | 0.044657% |
| 364 | Petroclostridium sp. X23 | 3045146 | S | 29 | 29 | 0.044657% |
| 365 | Arcobacter aquimarinus | 1315211 | S | 29 | 29 | 0.044657% |
| 366 | Tissierella sp. MB52-C2 | 3070999 | S | 29 | 29 | 0.044657% |
| 367 | Fusobacterium animalis | 76859 | S | 29 | 29 | 0.044657% |
| 368 | Treponema succinifaciens | 167 | S | 29 | 29 | 0.044657% |
| 369 | Lachnoclostridium phytofermentans | 66219 | S | 29 | 29 | 0.044657% |
| 370 | Blautia hydrogenotrophica | 53443 | S | 29 | 29 | 0.044657% |
| 371 | Monoglobus pectinilyticus | 1981510 | S | 28 | 28 | 0.043117% |
| 372 | Nitratidesulfovibrio liaohensis | 2604158 | S | 28 | 28 | 0.043117% |
| 373 | Leptotrichia sp. oral taxon 847 | 1785996 | S | 28 | 28 | 0.043117% |
| 374 | Ruminococcus gauvreauii | 438033 | S | 28 | 28 | 0.043117% |
| 375 | Halobacteriovorax sp. BALOs_7 | 2109558 | S | 28 | 28 | 0.043117% |
| 376 | Treponema bryantii | 163 | S | 28 | 28 | 0.043117% |
| 377 | Rhodoferax sp. BAB1 | 2741720 | S | 28 | 28 | 0.043117% |
| 378 | Sulfurospirillum barnesii | 44674 | S | 28 | 28 | 0.043117% |
| 379 | Faecalibacterium duncaniae | 411483 | S | 27 | 27 | 0.041577% |
| 380 | Silvimonas soli | 2980100 | S | 27 | 27 | 0.041577% |
| 381 | Thauera sp. JM12B12 | 3142262 | S | 27 | 27 | 0.041577% |
| 382 | Thomasclavelia ramosa | 1547 | S | 27 | 27 | 0.041577% |
| 383 | Butyrivibrio fibrisolvens | 831 | S | 27 | 27 | 0.041577% |
| 384 | Aeromonas simiae | 218936 | S | 27 | 27 | 0.041577% |
| 385 | Alistipes finegoldii | 214856 | S | 27 | 27 | 0.041577% |
| 386 | Ramlibacter tataouinensis | 94132 | S | 27 | 27 | 0.041577% |
| 387 | Mageeibacillus indolicus | 884684 | S | 27 | 27 | 0.041577% |
| 388 | Faecalibaculum rodentium | 1702221 | S | 26 | 26 | 0.040037% |
| 389 | Comamonas thiooxydans | 363952 | S | 26 | 26 | 0.040037% |
| 390 | Rheinheimera sp. D18 | 2545632 | S | 26 | 26 | 0.040037% |
| 391 | Acinetobacter bouvetii | 202951 | S | 26 | 26 | 0.040037% |
| 392 | Opitutus sp. GAS368 | 1882749 | S | 26 | 26 | 0.040037% |
| 393 | Enterocloster bolteae | 208479 | S | 26 | 26 | 0.040037% |
| 394 | Pedobacter sp. SL55 | 2995161 | S | 26 | 26 | 0.040037% |
| 395 | Faecalibacterium sp. I3-3-33 | 2929492 | S | 26 | 26 | 0.040037% |
| 396 | Bacterioplanes sanyensis | 1249553 | S | 26 | 26 | 0.040037% |
| 397 | Megalodesulfovibrio gigas | 879 | S | 26 | 26 | 0.040037% |
| 398 | Helicobacter hepaticus | 32025 | S | 26 | 26 | 0.040037% |
| 399 | Hirschia baltica | 2724 | S | 26 | 26 | 0.040037% |
| 400 | Crenobacter cavernae | 2290923 | S | 25 | 25 | 0.038497% |
| 401 | Fontisphaera persica | 2974023 | S | 25 | 25 | 0.038497% |
| 402 | Chordicoccus furentiruminis | 2709410 | S | 25 | 25 | 0.038497% |
| 403 | Pajaroellobacter abortibovis | 1882918 | S | 25 | 25 | 0.038497% |
| 404 | Flagellatimonas centrodinii | 2806210 | S | 25 | 25 | 0.038497% |
| 405 | Streptococcus suis | 1307 | S | 25 | 25 | 0.038497% |
| 406 | Subdoligranulum variabile | 214851 | S | 24 | 24 | 0.036957% |
| 407 | Dechloromonas sp. TW-R-39-2 | 2654218 | S | 24 | 24 | 0.036957% |
| 408 | Candidatus Pseudomonas adelgestsugas | 1302376 | S | 24 | 24 | 0.036957% |
| 409 | Treponema parvum | 138851 | S | 24 | 24 | 0.036957% |
| 410 | Pectinatus frisingensis | 865 | S | 24 | 24 | 0.036957% |
| 411 | Helicobacter pylori | 210 | S | 24 | 24 | 0.036957% |
| 412 | Gemmobacter fulvus | 2840474 | S | 24 | 24 | 0.036957% |
| 413 | Caproiciproducens sp. CPB-2 | 3030017 | S | 24 | 24 | 0.036957% |
| 414 | Ottowia sp. SB7-C50 | 3081231 | S | 24 | 24 | 0.036957% |
| 415 | Larkinella insperata | 332158 | S | 24 | 24 | 0.036957% |
| 416 | Verminephrobacter eiseniae | 364317 | S | 24 | 24 | 0.036957% |
| 417 | Methylobacterium sp. FF17 | 2984843 | S | 23 | 23 | 0.035417% |
| 418 | Acinetobacter sp. WCHA45 | 2004644 | S | 23 | 23 | 0.035417% |
| 419 | Aeromonas sp. CU5 | 2033033 | S | 23 | 23 | 0.035417% |
| 420 | Enterobacter ludwigii | 299767 | S | 23 | 23 | 0.035417% |
| 421 | Amniculibacterium sp. G2-70 | 2767188 | S | 23 | 23 | 0.035417% |
| 422 | Candidatus Methanoplasma termitum | 1577791 | S | 23 | 23 | 0.035417% |
| 423 | Sandaracinus amylolyticus | 927083 | S | 23 | 23 | 0.035417% |
| 424 | Candidatus Promineifilum breve | 1806508 | S | 23 | 23 | 0.035417% |
| 425 | Faecalibacterium sp. I3-3-89 | 2929493 | S | 23 | 23 | 0.035417% |
| 426 | Romboutsia sp. 13368 | 2708053 | S | 23 | 23 | 0.035417% |
| 427 | Anaeropeptidivorans aminofermentans | 2934315 | S | 23 | 23 | 0.035417% |
| 428 | Agarilytica rhodophyticola | 1737490 | S | 23 | 23 | 0.035417% |
| 429 | Hoylesella buccalis | 28127 | S | 23 | 23 | 0.035417% |
| 430 | Clostridioides difficile | 1496 | S | 23 | 23 | 0.035417% |
| 431 | Hungatella hathewayi | 154046 | S | 23 | 23 | 0.035417% |
| 432 | Paracoccus aminophilus | 34003 | S | 23 | 23 | 0.035417% |
| 433 | Synechococcus sp. JA-3-3Ab | 321327 | S | 22 | 22 | 0.033878% |
| 434 | Fuerstiella marisgermanici | 1891926 | S | 22 | 22 | 0.033878% |
| 435 | Humisphaera borealis | 2807512 | S | 22 | 22 | 0.033878% |
| 436 | Streptococcus merionis | 400065 | S | 22 | 22 | 0.033878% |
| 437 | Lactococcus raffinolactis | 1366 | S | 22 | 22 | 0.033878% |
| 438 | [Clostridium] innocuum | 1522 | S | 22 | 22 | 0.033878% |
| 439 | Bdellovibrio sp. 22V | 3044166 | S | 22 | 22 | 0.033878% |
| 440 | Wansuia hejianensis | 2763667 | S | 22 | 22 | 0.033878% |
| 441 | Beggiatoa leptomitoformis | 288004 | S | 22 | 22 | 0.033878% |
| 442 | Caproicibacterium sp. BJN0003 | 2994078 | S | 22 | 22 | 0.033878% |
| 443 | Ereboglobus luteus | 1796921 | S | 22 | 22 | 0.033878% |
| 444 | Petrimonas mucosa | 1642646 | S | 22 | 22 | 0.033878% |
| 445 | Sorangium cellulosum | 56 | S | 22 | 22 | 0.033878% |
| 446 | Lachnoclostridium phocaeense | 1871021 | S | 22 | 22 | 0.033878% |
| 447 | Jeongeupia sp. HS-3 | 1009682 | S | 22 | 22 | 0.033878% |
| 448 | Fusobacterium polymorphum | 76857 | S | 22 | 22 | 0.033878% |
| 449 | Candidatus Karelsulcia muelleri | 336810 | S | 22 | 22 | 0.033878% |
| 450 | Neisseria gonorrhoeae | 485 | S | 22 | 22 | 0.033878% |
| 451 | Enterocloster asparagiformis | 333367 | S | 22 | 22 | 0.033878% |
| 452 | Roseburia sp. 499 | 1261634 | S | 21 | 21 | 0.032338% |
| 453 | Capsulimonas corticalis | 2219043 | S | 21 | 21 | 0.032338% |
| 454 | Stomatobaculum sp. F0698 | 3059030 | S | 21 | 21 | 0.032338% |
| 455 | Oscillibacter sp. PEA192 | 2109687 | S | 21 | 21 | 0.032338% |
| 456 | Tichowtungia aerotolerans | 2697043 | S | 21 | 21 | 0.032338% |
| 457 | Aliarcobacter cibarius | 255507 | S | 21 | 21 | 0.032338% |
| 458 | uncultured Marinifilum sp. | 1127472 | S | 21 | 21 | 0.032338% |
| 459 | Faecalibacterium sp. IP-1-18 | 2929488 | S | 21 | 21 | 0.032338% |
| 460 | Azoarcus olearius | 418699 | S | 21 | 21 | 0.032338% |
| 461 | Porphyromonas gingivalis | 837 | S | 21 | 21 | 0.032338% |
| 462 | Haematospirillum jordaniae | 1549855 | S | 21 | 21 | 0.032338% |
| 463 | Desulfoluna limicola | 2810562 | S | 21 | 21 | 0.032338% |
| 464 | Chakrabartyella piscis | 2918914 | S | 21 | 21 | 0.032338% |
| 465 | Alistipes senegalensis | 1288121 | S | 21 | 21 | 0.032338% |
| 466 | uncultured Caproiciproducens sp. | 1738647 | S | 21 | 21 | 0.032338% |
| 467 | Methylogaea oryzae | 1295382 | S | 20 | 20 | 0.030798% |
| 468 | Zobellella denitrificans | 347534 | S | 20 | 20 | 0.030798% |
| 469 | Heliorestis convoluta | 356322 | S | 20 | 20 | 0.030798% |
| 470 | Bacteroides zoogleoformans | 28119 | S | 20 | 20 | 0.030798% |
| 471 | Piscinibacter gummiphilus | 946333 | S | 20 | 20 | 0.030798% |
| 472 | Clostridium butyricum | 1492 | S | 20 | 20 | 0.030798% |
| 473 | Pseudobacter ginsenosidimutans | 661488 | S | 20 | 20 | 0.030798% |
| 474 | Hankyongella ginsenosidimutans | 1763828 | S | 20 | 20 | 0.030798% |
| 475 | Spirosoma pollinicola | 2057025 | S | 20 | 20 | 0.030798% |
| 476 | Treponema peruense | 2787628 | S | 20 | 20 | 0.030798% |
| 477 | Campylobacter gracilis | 824 | S | 20 | 20 | 0.030798% |
| 478 | Casimicrobium huifangae | 2591109 | S | 20 | 20 | 0.030798% |
| 479 | Leptotrichia hofstadii | 157688 | S | 20 | 20 | 0.030798% |
| 480 | Candidatus Portiera aleyrodidarum | 91844 | S | 20 | 20 | 0.030798% |
| 481 | Indioceanicola profundi | 2220096 | S | 20 | 20 | 0.030798% |
| 482 | Proteiniphilum propionicum | 2829812 | S | 20 | 20 | 0.030798% |
| 483 | Horticoccus luteus | 2862869 | S | 20 | 20 | 0.030798% |
| 484 | Sulfuriroseicoccus oceanibius | 2707525 | S | 20 | 20 | 0.030798% |
| 485 | Tannerella serpentiformis | 712710 | S | 20 | 20 | 0.030798% |
| 486 | Vibrio harveyi | 669 | S | 20 | 20 | 0.030798% |
| 487 | Photobacterium damselae | 38293 | S | 20 | 20 | 0.030798% |
| 488 | Vibrio cholerae | 666 | S | 20 | 20 | 0.030798% |
| 489 | Pseudomonas chlororaphis | 587753 | S | 20 | 20 | 0.030798% |
| 490 | Anaerotignum propionicum | 28446 | S | 20 | 20 | 0.030798% |

---

# Automatically Selected Bacteriophage Species: Reads and Abundance

A total of 67 Bacteriophage species were selected by the automated pipeline from sample
HTP-2024-01-15 (Hyperion Water Reclamation Plant, Los Angeles, California, USA; 24-h composite influent; rRNA-depleted library; January 15, 2024, peak of 2023–2024 disease season).
Selection criteria: NCBI rank = `species` and clade reads >= 20 in the sample.
Actual reads (clade) from Kraken2 classification are used directly as
simulation design reads; the two columns are identical for all species.

| No. | Name | NCBI taxid | NCBI rank | Actual reads (HTP-2024-01-15, clade) | Simulation design reads | Relative abundance in simulation (%) |
|-----|------|-----------|-----------|--------------------------------------|------------------------|-------------------------------------|
| 1 | ssRNA phage SRR6960799_40 | 2786599 | S1 | 1625 | 1625 | 1.272763% |
| 2 | ssRNA phage SRR6960803_7 | 2786623 | S1 | 1181 | 1181 | 0.925005% |
| 3 | ssRNA phage SRR6049586_2 | 2786481 | S1 | 663 | 663 | 0.519287% |
| 4 | ssRNA phage SRR7976357_7 | 2786747 | S1 | 561 | 561 | 0.439397% |
| 5 | ssRNA phage SRR7976300_2 | 2786651 | S1 | 474 | 474 | 0.371255% |
| 6 | ssRNA phage SRR5466727_7 | 2786436 | S1 | 444 | 444 | 0.347758% |
| 7 | Acinetobacter phage AP205 | 154784 | S1 | 378 | 378 | 0.296064% |
| 8 | ssRNA phage SRR5466727_2 | 2786431 | S1 | 311 | 311 | 0.243587% |
| 9 | ssRNA phage SRR6960803_4 | 2786620 | S1 | 269 | 269 | 0.210691% |
| 10 | ssRNA phage SRR6960507_10 | 2786508 | S1 | 269 | 269 | 0.210691% |
| 11 | ssRNA phage SRR5466369_1 | 2786404 | S1 | 267 | 267 | 0.209125% |
| 12 | ssRNA phage SRR6960509_17 | 2786528 | S1 | 232 | 232 | 0.181711% |
| 13 | ssRNA phage SRR7976300_8 | 2786657 | S1 | 221 | 221 | 0.173096% |
| 14 | ssRNA phage SRR5466727_4 | 2786433 | S1 | 194 | 194 | 0.151948% |
| 15 | ssRNA phage Zoerhiza.1_5 | 2786783 | S1 | 155 | 155 | 0.121402% |
| 16 | ssRNA phage SRR5466337_3 | 2786389 | S1 | 149 | 149 | 0.116703% |
| 17 | ssRNA phage SRR7976299_5 | 2786645 | S1 | 147 | 147 | 0.115136% |
| 18 | ssRNA phage SRR7976357_8 | 2786748 | S1 | 139 | 139 | 0.108870% |
| 19 | ssRNA phage SRR6960509_10 | 2786521 | S1 | 119 | 119 | 0.093205% |
| 20 | ssRNA phage Esthiorhiza.1_1 | 2786011 | S1 | 102 | 102 | 0.079890% |
| 21 | ssRNA phage SRR6960799_25 | 2786582 | S1 | 101 | 101 | 0.079107% |
| 22 | ssRNA phage SRR7976310_16 | 2786678 | S1 | 98 | 98 | 0.076757% |
| 23 | ssRNA phage SRR5467137_1 | 2786475 | S1 | 86 | 86 | 0.067359% |
| 24 | ssRNA phage SRR7976325_7 | 2786722 | S1 | 84 | 84 | 0.065792% |
| 25 | ssRNA phage SRR6960549_3 | 2786540 | S1 | 81 | 81 | 0.063442% |
| 26 | ssRNA phage Esthiorhiza.1_2 | 2786017 | S1 | 75 | 75 | 0.058743% |
| 27 | ssRNA phage SRR7976299_8 | 2786648 | S1 | 73 | 73 | 0.057176% |
| 28 | Enterobacteria phage GA | 12018 | S1 | 73 | 73 | 0.057176% |
| 29 | ssRNA phage ESE006 | 2785996 | S1 | 72 | 72 | 0.056393% |
| 30 | ssRNA phage Zoerhiza.1_25 | 2786767 | S1 | 68 | 68 | 0.053260% |
| 31 | ssRNA phage SRR7976310_7 | 2786685 | S1 | 67 | 67 | 0.052477% |
| 32 | ssRNA phage SRR5467091_15 | 2786465 | S1 | 65 | 65 | 0.050911% |
| 33 | ssRNA phage SRR5466364_4 | 2786398 | S1 | 65 | 65 | 0.050911% |
| 34 | ssRNA phage SRR6960509_1 | 2786520 | S1 | 65 | 65 | 0.050911% |
| 35 | ssRNA phage Gephyllon.1_24 | 2786124 | S1 | 64 | 64 | 0.050127% |
| 36 | ssRNA phage SRR7976310_10 | 2786672 | S1 | 62 | 62 | 0.048561% |
| 37 | ssRNA phage SRR5466725_16 | 2786414 | S1 | 60 | 60 | 0.046994% |
| 38 | ssRNA phage Zoerhiza.4_12 | 2786828 | S1 | 60 | 60 | 0.046994% |
| 39 | ssRNA phage ESE029 | 2786005 | S1 | 57 | 57 | 0.044645% |
| 40 | ssRNA phage SRR6960799_27 | 2786584 | S1 | 55 | 55 | 0.043078% |
| 41 | ssRNA phage SRR6960799_30 | 2786588 | S1 | 52 | 52 | 0.040728% |
| 42 | ssRNA phage SRR6960797_1 | 2786561 | S1 | 50 | 50 | 0.039162% |
| 43 | ssRNA phage Gephyllon.2_8 | 2786147 | S1 | 49 | 49 | 0.038379% |
| 44 | ssRNA phage SRR7976357_6 | 2786746 | S1 | 47 | 47 | 0.036812% |
| 45 | ssRNA phage SRR7976327_1 | 2786730 | S1 | 46 | 46 | 0.036029% |
| 46 | ssRNA phage SRR6960799_17 | 2786573 | S1 | 43 | 43 | 0.033679% |
| 47 | ssRNA phage SRR6960799_10 | 2786566 | S1 | 41 | 41 | 0.032113% |
| 48 | ssRNA phage SRR6960509_7 | 2786534 | S1 | 38 | 38 | 0.029763% |
| 49 | ssRNA phage SRR7976325_29 | 2786717 | S1 | 35 | 35 | 0.027413% |
| 50 | ssRNA phage SRR6960799_15 | 2786571 | S1 | 35 | 35 | 0.027413% |
| 51 | ssRNA phage SRR6255513_1 | 2786496 | S1 | 34 | 34 | 0.026630% |
| 52 | ssRNA phage SRR5466365_2 | 2786401 | S1 | 33 | 33 | 0.025847% |
| 53 | ssRNA phage SRR6253161_3 | 2786490 | S1 | 32 | 32 | 0.025064% |
| 54 | ssRNA phage SRR5467091_3 | 2786468 | S1 | 32 | 32 | 0.025064% |
| 55 | ssRNA phage Zoerhiza.2_22 | 2786802 | S1 | 31 | 31 | 0.024280% |
| 56 | ssRNA phage SRR6960799_36 | 2786594 | S1 | 28 | 28 | 0.021931% |
| 57 | ssRNA phage SRR7976326_4 | 2786728 | S1 | 28 | 28 | 0.021931% |
| 58 | ssRNA phage Gerhypos.3_16 | 2786289 | S1 | 28 | 28 | 0.021931% |
| 59 | Caulobacter phage phiCb5 | 767473 | S1 | 27 | 27 | 0.021147% |
| 60 | Enterobacteria phage Hgal1 | 1206300 | S1 | 26 | 26 | 0.020364% |
| 61 | ssRNA phage SRR6960551_9 | 2786560 | S1 | 24 | 24 | 0.018798% |
| 62 | ssRNA phage SRR7976357_2 | 2786742 | S1 | 22 | 22 | 0.017231% |
| 63 | ssRNA phage SRR5467090_10 | 2786448 | S1 | 21 | 21 | 0.016448% |
| 64 | ssRNA phage AVE016 | 2785988 | S1 | 21 | 21 | 0.016448% |
| 65 | ssRNA phage Gerhypos.1_12 | 2786193 | S1 | 20 | 20 | 0.015665% |
| 66 | ssRNA phage Esthiorhiza.2_21 | 2786038 | S1 | 20 | 20 | 0.015665% |
| 67 | ssRNA phage Gephyllon.4_11 | 2786169 | S1 | 20 | 20 | 0.015665% |

---

# Automatically Selected Virus Species: Reads and Abundance

A total of 138 viral species were selected by the automated pipeline from sample
HTP-2024-01-15 (Hyperion Water Reclamation Plant, Los Angeles, California, USA; 24-h composite influent; rRNA-depleted library; January 15, 2024, peak of 2023–2024 disease season).
Selection criteria: NCBI rank = `species` and clade reads >= 20 in the sample.
Actual reads (clade) from Kraken2 classification are used directly as
simulation design reads; the two columns are identical for all 138 species.

| No. | Virus | NCBI taxid | NCBI rank | Actual reads (HTP-2024-01-15, clade) | Simulation design reads | Relative abundance in simulation (%) |
|-----|-------|-----------|-----------|--------------------------------------|------------------------|-------------------------------------|
| 1 | Mamastrovirus 1 | 1239565 | species | 559347 | 559347 | 2.355158% |
| 2 | Norwalk virus | 11983 | species | 83135 | 83135 | 0.350044% |
| 3 | Sapporo virus | 95342 | species | 70569 | 70569 | 0.297134% |
| 4 | Mamastrovirus 6 | 1239570 | species | 57851 | 57851 | 0.243584% |
| 5 | Mamastrovirus 9 | 1239573 | species | 21096 | 21096 | 0.088826% |
| 6 | Astrovirus MLB2 | 683172 | species | 19441 | 19441 | 0.081857% |
| 7 | HMO Astrovirus A | 682382 | species | 18930 | 18930 | 0.079706% |
| 8 | Husavirus sp. | 1930509 | species | 14941 | 14941 | 0.062910% |
| 9 | Kobuvirus aichi | 72149 | species | 13393 | 13393 | 0.056392% |
| 10 | Orthopicobirnavirus hominis | 2956252 | species | 8992 | 8992 | 0.037861% |
| 11 | Porcine picobirnavirus | 442302 | species | 8939 | 8939 | 0.037638% |
| 12 | Canine astrovirus | 1157338 | species | 7194 | 7194 | 0.030291% |
| 13 | Salivirus A | 1330524 | species | 6118 | 6118 | 0.025760% |
| 14 | Enterovirus C | 138950 | species | 5222 | 5222 | 0.021987% |
| 15 | Hudisavirus sp. | 2021738 | species | 4786 | 4786 | 0.020152% |
| 16 | Astrovirus VA3 | 683174 | species | 3648 | 3648 | 0.015360% |
| 17 | Genomoviridae sp. | 2202565 | species | 3278 | 3278 | 0.013802% |
| 18 | Porcine torovirus | 237020 | species | 3140 | 3140 | 0.013221% |
| 19 | Astrovirus MLB3 | 1247114 | species | 1891 | 1891 | 0.007962% |
| 20 | Alphacoronavirus 1 | 693997 | species | 1860 | 1860 | 0.007832% |
| 21 | Human feces pecovirus | 1820160 | species | 1810 | 1810 | 0.007621% |
| 22 | Picobirnavirus sp. | 1907787 | species | 1632 | 1632 | 0.006872% |
| 23 | Picornavirales Tottori-HG1 | 1795648 | species | 1545 | 1545 | 0.006505% |
| 24 | Rotavirus A | 28875 | species | 1509 | 1509 | 0.006354% |
| 25 | Bovine torovirus | 74501 | species | 1456 | 1456 | 0.006131% |
| 26 | Rocahepevirus ratti | 1678145 | species | 1365 | 1365 | 0.005747% |
| 27 | Feline astrovirus 2 | 1389204 | species | 1356 | 1356 | 0.005710% |
| 28 | Picobirnavirus Equ3 | 1673646 | species | 1191 | 1191 | 0.005015% |
| 29 | Mamastrovirus 2 | 1239566 | species | 1104 | 1104 | 0.004648% |
| 30 | Goose astrovirus | 1349999 | species | 983 | 983 | 0.004139% |
| 31 | Severe acute respiratory syndrome-related coronavirus | 694009 | species | 934 | 934 | 0.003933% |
| 32 | Enterovirus A | 138948 | species | 844 | 844 | 0.003554% |
| 33 | Sicinivirus A | 1755588 | species | 726 | 726 | 0.003057% |
| 34 | Cressdnaviricota sp. | 2748378 | species | 719 | 719 | 0.003027% |
| 35 | Cosavirus A | 1330491 | species | 703 | 703 | 0.002960% |
| 36 | Otarine picobirnavirus | 1187973 | species | 660 | 660 | 0.002779% |
| 37 | Enterovirus B | 138949 | species | 637 | 637 | 0.002682% |
| 38 | Avian coronavirus | 694014 | species | 623 | 623 | 0.002623% |
| 39 | Cardiovirus B | 1821750 | species | 622 | 622 | 0.002619% |
| 40 | Canine kobuvirus | 1836608 | species | 619 | 619 | 0.002606% |
| 41 | Chicken picobirnavirus | 1930304 | species | 587 | 587 | 0.002472% |
| 42 | Primate norovirus | 1903318 | species | 584 | 584 | 0.002459% |
| 43 | Megrivirus C | 2170006 | species | 577 | 577 | 0.002429% |
| 44 | Parechovirus A | 1803956 | species | 553 | 553 | 0.002328% |
| 45 | Mischivirus D | 2560598 | species | 521 | 521 | 0.002194% |
| 46 | Rotavirus C | 36427 | species | 510 | 510 | 0.002147% |
| 47 | Astrovirus VA4 | 1247113 | species | 490 | 490 | 0.002063% |
| 48 | Rhinovirus C | 463676 | species | 459 | 459 | 0.001933% |
| 49 | Human fecal virus Jorvi2 | 2017081 | species | 441 | 441 | 0.001857% |
| 50 | Rhinovirus A | 147711 | species | 427 | 427 | 0.001798% |
| 51 | Mouse astrovirus M-52/USA/2008 | 1074210 | species | 412 | 412 | 0.001735% |
| 52 | Bovine astrovirus B18/HK | 1027245 | species | 381 | 381 | 0.001604% |
| 53 | Rosavirus B | 1902501 | species | 345 | 345 | 0.001453% |
| 54 | Circular ssDNA virus sp. | 2805939 | species | 315 | 315 | 0.001326% |
| 55 | Cardiovirus C | 1821751 | species | 305 | 305 | 0.001284% |
| 56 | Bocaparvovirus carnivoran2 | 3052027 | species | 295 | 295 | 0.001242% |
| 57 | Bocaparvovirus primate2 | 3052041 | species | 272 | 272 | 0.001145% |
| 58 | Human cosavirus | 1233383 | species | 269 | 269 | 0.001133% |
| 59 | Kobuvirus bejaponia | 194965 | species | 261 | 261 | 0.001099% |
| 60 | Murine astrovirus | 1141625 | species | 240 | 240 | 0.001011% |
| 61 | Avastrovirus 2 | 1239438 | species | 212 | 212 | 0.000893% |
| 62 | Gemykibivirus humas5 | 2844585 | species | 203 | 203 | 0.000855% |
| 63 | Canine vesivirus | 1176484 | species | 195 | 195 | 0.000821% |
| 64 | Mamastrovirus 3 | 1239567 | species | 185 | 185 | 0.000779% |
| 65 | Mamastrovirus 13 | 1239577 | species | 184 | 184 | 0.000775% |
| 66 | Protoparvovirus carnivoran1 | 3052707 | species | 178 | 178 | 0.000749% |
| 67 | Cardiovirus D | 2734523 | species | 166 | 166 | 0.000699% |
| 68 | Chicken picornavirus 4 | 1534547 | species | 165 | 165 | 0.000695% |
| 69 | Human smacovirus 1 | 1595998 | species | 155 | 155 | 0.000653% |
| 70 | Betacoronavirus 1 | 694003 | species | 150 | 150 | 0.000632% |
| 71 | Dependoparvovirus primate1 | 1511891 | species | 146 | 146 | 0.000615% |
| 72 | Rhinovirus B | 147712 | species | 138 | 138 | 0.000581% |
| 73 | Human mastadenovirus F | 130309 | species | 137 | 137 | 0.000577% |
| 74 | Kobuvirus cebes | 1298633 | species | 136 | 136 | 0.000573% |
| 75 | Tioga picorna-like virus 1 | 2029771 | species | 133 | 133 | 0.000560% |
| 76 | Giant panda associated gemycircularvirus | 2016461 | species | 133 | 133 | 0.000560% |
| 77 | Astrovirus wild boar/WBAstV-1/2011/HUN | 1137931 | species | 132 | 132 | 0.000556% |
| 78 | Melegrivirus A | 1330070 | species | 128 | 128 | 0.000539% |
| 79 | Hunnivirus A | 1431460 | species | 124 | 124 | 0.000522% |
| 80 | Bocaparvovirus primate1 | 3052040 | species | 115 | 115 | 0.000484% |
| 81 | Goat torovirus | 1837217 | species | 115 | 115 | 0.000484% |
| 82 | Chicken associated smacovirus | 1932006 | species | 109 | 109 | 0.000459% |
| 83 | Cosavirus D | 2003650 | species | 100 | 100 | 0.000421% |
| 84 | Nordland virus | 2560615 | species | 95 | 95 | 0.000400% |
| 85 | Bocaparvovirus carnivoran1 | 3052026 | species | 94 | 94 | 0.000396% |
| 86 | Avian leukosis virus | 11864 | species | 91 | 91 | 0.000383% |
| 87 | Chicken picornavirus 1 | 1534544 | species | 86 | 86 | 0.000362% |
| 88 | Rodent astrovirus | 1914442 | species | 80 | 80 | 0.000337% |
| 89 | Human coronavirus 229E | 11137 | species | 78 | 78 | 0.000328% |
| 90 | Gallivirus Pf-CHK1/GV | 1766557 | species | 69 | 69 | 0.000291% |
| 91 | Ferret coronavirus | 1264898 | species | 66 | 66 | 0.000278% |
| 92 | Mamastrovirus 10 | 1239574 | species | 63 | 63 | 0.000265% |
| 93 | Cosavirus E | 2003651 | species | 57 | 57 | 0.000240% |
| 94 | Qinghai Himalayan marmot astrovirus 1 | 1961665 | species | 50 | 50 | 0.000211% |
| 95 | Bastrovirus BAS-3 | 2025398 | species | 47 | 47 | 0.000198% |
| 96 | Goose dicistrovirus | 1776109 | species | 47 | 47 | 0.000198% |
| 97 | Porcine stool-associated circular virus 4 | 1475062 | species | 47 | 47 | 0.000198% |
| 98 | Enterovirus D | 138951 | species | 46 | 46 | 0.000194% |
| 99 | Meles meles circovirus-like virus | 1128073 | species | 46 | 46 | 0.000194% |
| 100 | Calicivirus chicken/V0021/Bayern/2004 | 889389 | species | 45 | 45 | 0.000189% |
| 101 | Tremovirus A | 70796 | species | 45 | 45 | 0.000189% |
| 102 | Hepatovirus A | 12092 | species | 44 | 44 | 0.000185% |
| 103 | Bavaria virus | 2560343 | species | 43 | 43 | 0.000181% |
| 104 | Cyclovirus humana | 3047957 | species | 40 | 40 | 0.000168% |
| 105 | Rabovirus A | 1603962 | species | 40 | 40 | 0.000168% |
| 106 | Burkina Faso astrovirus | 1518574 | species | 40 | 40 | 0.000168% |
| 107 | Gemykibivirus haeme2 | 2844576 | species | 38 | 38 | 0.000160% |
| 108 | Canarypox virus | 44088 | species | 37 | 37 | 0.000156% |
| 109 | Blattambidensovirus incertum1 | 3059782 | species | 37 | 37 | 0.000156% |
| 110 | Cripavirus NB-1/2011/HUN | 1541211 | species | 37 | 37 | 0.000156% |
| 111 | Rabbit astrovirus TN/2208/2010 | 1161132 | species | 37 | 37 | 0.000156% |
| 112 | Betapolyomavirus hominis | 1891762 | species | 36 | 36 | 0.000152% |
| 113 | Gemycircularvirus monas1 | 2844484 | species | 35 | 35 | 0.000147% |
| 114 | Gemycircularvirus porci1 | 2844505 | species | 34 | 34 | 0.000143% |
| 115 | Sakobuvirus A | 1659771 | species | 34 | 34 | 0.000143% |
| 116 | Anativirus biphaco | 2734519 | species | 33 | 33 | 0.000139% |
| 117 | Porprismacovirus humas3 | 2845899 | species | 30 | 30 | 0.000126% |
| 118 | Gemykibivirus canfam1 | 2844564 | species | 29 | 29 | 0.000122% |
| 119 | Posavirus sp. | 1930510 | species | 29 | 29 | 0.000122% |
| 120 | Avisivirus C | 1985983 | species | 28 | 28 | 0.000118% |
| 121 | Posavirus 1 | 1105380 | species | 28 | 28 | 0.000118% |
| 122 | Sripuvirus almpiwar | 1985693 | species | 27 | 27 | 0.000114% |
| 123 | Picornavirales Bu-3 | 1813613 | species | 27 | 27 | 0.000114% |
| 124 | Gemykibivirus cynas1 | 2844567 | species | 26 | 26 | 0.000109% |
| 125 | Human associated porprismacovirus 2 | 2170114 | species | 26 | 26 | 0.000109% |
| 126 | Bovine astrovirus B76-2/HK | 1027247 | species | 26 | 26 | 0.000109% |
| 127 | Chicken stool-associated circular virus | 1930302 | species | 25 | 25 | 0.000105% |
| 128 | Cadicivirus A | 1330068 | species | 25 | 25 | 0.000105% |
| 129 | Qinghai Himalayan marmot astrovirus 2 | 1961666 | species | 23 | 23 | 0.000097% |
| 130 | Avisivirus Pf-CHK1/AsV | 1766559 | species | 23 | 23 | 0.000097% |
| 131 | Pteropine orthoreovirus | 1459003 | species | 23 | 23 | 0.000097% |
| 132 | Aalivirus apekidu | 2169685 | species | 22 | 22 | 0.000093% |
| 133 | Porcine astrovirus 3 | 1254425 | species | 22 | 22 | 0.000093% |
| 134 | Astrovirus Er/SZAL6/HUN/2011 | 1671382 | species | 21 | 21 | 0.000088% |
| 135 | Rous sarcoma virus | 11886 | species | 21 | 21 | 0.000088% |
| 136 | Inpeasmacovirus humas1 | 2845010 | species | 20 | 20 | 0.000084% |
| 137 | Gemycircularvirus lamas1 | 2844465 | species | 20 | 20 | 0.000084% |
| 138 | Betapolyomavirus secuhominis | 1891763 | species | 20 | 20 | 0.000084% |

---

# Manually Added Viruses: Selection Rationale

Thirteen viral species were manually added to the simulation community
(`species_list_final.tsv`) after the automated pipeline produced the
primary species list from sample HTP-2024-01-15 (Hyperion Water Reclamation Plant, Los Angeles, California, USA; 24-h composite influent; rRNA-depleted library; January 15, 2024, peak of 2023–2024 disease season). These additions fell into three
categories based on the reason they were absent from the automated list.

---

## Category 1: Present in sample but excluded by taxonomic rank filter

The automated pipeline retained only taxa classified at the NCBI
`species` rank. Two high-abundance viruses were excluded because NCBI
assigns them a `no rank` status — they are strains beneath a broader
species node rather than species themselves.

| Virus | NCBI taxid | NCBI rank | Actual reads (HTP-2024-01-15, direct) | Simulation design reads | Relative abundance in simulation (%) | Parent species node |
|-------|-----------|-----------|--------------------------------------|------------------------|-------------------------------------|---------------------|
| Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) | 2697049 | no rank | 875 | 875 | 0.003684% | *Severe acute respiratory syndrome-related coronavirus* (taxid 694009, species) |
| Human coronavirus OC43 | 31631 | no rank | 47 | 47 | 0.000198% | *Betacoronavirus 1* (taxid 694003, species) |

Simulation design reads match the actual Kraken2-classified reads from
sample HTP-2024-01-15, as these viruses were genuinely abundant in the
sample.

---

## Category 2: Present in sample but excluded by read-count filter

Taxa with fewer than 20 directly assigned reads were removed from the
automated list. Several respiratory viruses had detectable but very low
signals in the sample, and were additionally classified at non-species
ranks (serotype, no rank). Simulation design reads were manually
assigned to reflect expected winter epidemiological signal rather than
the near-zero sample counts.

| Virus | NCBI taxid | NCBI rank | Actual reads (HTP-2024-01-15, direct) | Actual reads (clade) | Simulation design reads | Relative abundance in simulation (%) |
|-------|-----------|-----------|--------------------------------------|---------------------|------------------------|-------------------------------------|
| Influenza A virus (H3N2) | 119210 | serotype | 0 | 1 | 80 | 0.000337% |
| Influenza A virus (H1N1) | 114727 | serotype | 0 | 6 | 60 | 0.000253% |
| Influenza B virus | 2955465 | species | 0 | 1 | 30 | 0.000126% |
| Human orthopneumovirus A (RSV A) | 11250 | no rank | 0 | 3 | 50 | 0.000211% |
| Human parainfluenza virus 3 (PIV3) | 11216 | no rank | 1 | 1 | 20 | 0.000084% |
| Monkeypox virus | 10244 | species | 6 | 6 | 6 | 0.000025% |

---

## Category 3: Not detected in sample — added for epidemiological completeness

Five viruses had no detectable reads in HTP-2024-01-15 under any taxon
node. They were included based on their public-health or environmental
surveillance relevance for a winter wastewater metagenome simulation.
Simulation design reads are entirely manually assigned.

| Virus | NCBI taxid | Actual reads (HTP-2024-01-15) | Simulation design reads | Relative abundance in simulation (%) | Rationale |
|-------|-----------|------------------------------|------------------------|--------------------------------------|-----------|
| Influenza A virus H5N1 (A/goose/Guangdong/1/1996) | 93838 | 0 (not detected) | 4 | 0.000017% | Early spillover sentinel: represents low-level zoonotic spillover signal relevant to HPAI H5N1 surveillance scenarios |
| Human metapneumovirus (hMPV) | 162145 | 0 (not detected) | 35 | 0.000147% | Winter peak respiratory virus; co-circulates with RSV during cold months |
| Human astrovirus 1 | 37108 | 0 (not detected) | 25 | 0.000105% | Ubiquitous enteric virus commonly detected in municipal wastewater; included as a representative gastrointestinal pathogen |
| Human papillomavirus 16 (HPV 16) | 333760 | 0 (not detected) | 8 | 0.000034% | Representative of low-abundance DNA viruses detectable in wastewater; included to broaden community diversity |
| Poliovirus 1 | 12110 | 0 (not detected) | 5 | 0.000021% | Environmental surveillance sentinel; routinely monitored in wastewater-based epidemiology programs globally |

---

## Abundance assignment

Relative abundances for all 13 viruses were set proportionally within
the 4% viral fraction allocated to the Viruses domain. The following
principles guided the simulation design read assignments:

- **SARS-CoV-2 and OC43**: design reads equal actual sample reads, as
  both viruses were genuinely abundant in HTP-2024-01-15
- **Influenza A (H3N2 > H1N1)**: reflects typical winter co-circulation
  dominance hierarchy
- **Influenza B**: lower than Influenza A, consistent with surveillance
  data from January 2024
- **RSV A > hMPV > PIV3**: winter peak ordering based on seasonal
  incidence patterns
- **Monkeypox**: design reads match the 6 actual sample reads
- **Astrovirus 1**: modest endemic signal appropriate for a common
  enteric virus
- **HPV 16 and Poliovirus 1**: minimal abundance, representing rare or
  sentinel-level detection

All abundances were re-normalized to sum to 1.0 across the full
708-genome community after manual additions.
