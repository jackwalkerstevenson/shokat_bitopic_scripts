# manual colors for plotting
# specify treatments to rename in plots. "old_name" = "new_name"
display_names_treatments <- c(
  "netgr_Asc30" = "net growth rate: asciminib 30 nM",
  "netgr_Asc10" = "net growth rate: asciminib 10 nM",
  "netgr_PL1" = "net growth rate: PonatiLink-2 1 nM",
  "netgr_PL5" = "net growth rate: PonatiLink-2 5 nM",
  "netgr_Pon0.1" = "net growth rate: ponatinib 0.1 nM",
  "netgr_Pon0.4" = "net growth rate: ponatinib 0.4 nM",
  "ponatinib + asciminib" = "pona+asc",
  "dasatinib + asciminib" = "dasa+asc",
  "PonatiLink-1-24" = "PonatiLink-1",
  "PonatiLink-2-7-10" = "PonatiLink-2"
)
color_map_treatments <- c(
  "ponatinib" = "#A40000", # apparent darkred from chimera
  "asciminib" = "#FF5300", # apparent orangered from chimera
  "ponatinib + asciminib" = "gold2",
  "dasatinib" = "mediumpurple",
  "PonatiLink-1-24" = "forestgreen",
  "PonatiLink-2-7-10" = "skyblue3"
)
shape_map_treatments <- c(
  "ponatinib" = "circle",
  "asciminib" = "triangle",
  "ponatinib + asciminib" = "square",
  "dasatinib" = "triangle open",
  "PonatiLink-1-24" = "circle cross",
  "PonatiLink-2-7-10" = "square open" # consistent across plots
)
# specify targets to rename in plots. "old_name" = "new_name"
display_names_targets <- c(
  "K562 T315I" = "K562 BE-T315I",
  "K562 pUltra control" = "control",
  "K562 pUltra BCR-ABL1 wt" = "wt",
  "K562 pUltra BCR-ABL1 E255V" = "E255V",
  "K562 pUltra BCR-ABL1 T315I" = "T315I"
)
color_map_targets <- c(
  # ABL1 cell lines-----------------------------------------------
  "K562 pUltra control" = "gray",
  "K562 pUltra BCR-ABL1 wt" = "forestgreen",
  "K562 pUltra BCR-ABL1 E255V" = "palegreen2",
  "K562 pUltra BCR-ABL1 T315I" = "skyblue2",
  # selectscreen kinases-----------------------------------------------
  "ABL1" = "forestgreen",
  "ABL1 T315I" = "skyblue2",
  # other----------------------------------------------------------------
  "HUVEC" = "royalblue",
  "HEK 293T" = "forestgreen"
)
shape_map_targets <- c(
  # ABL1 cell lines-----------------------------------------------
  "K562 pUltra control" = "circle open",
  "K562 pUltra BCR-ABL1 wt" = "circle",
  "K562 pUltra BCR-ABL1 E255V" = "triangle",
  "K562 pUltra BCR-ABL1 T315I" = "square",
  # selectscreen kinases-----------------------------------------------
  "ABL1" = "circle",
  "ABL1 T315I" = "square",
  # other----------------------------------------------------------------
  "HUVEC" = "circle",
  "HEK 293T" = "circle"
)
