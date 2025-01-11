# manual colors for plotting
# specify treatments to rename in plots. "old_name" = "new_name"
display_names_treatments <- c(
  "ponatinib + asciminib" = "pona+asc",
  "RMC-5552" = "RMC-5552",
  "PonatiLink-1-12" = "PonatiLink-1-PEG12",
  "PonatiLink-1-16" = "PonatiLink-1-PEG16",
  "PonatiLink-1-20" = "PonatiLink-1-PEG20",
  # "PonatiLink-1-24" = "PonatiLink-1-PEG24",
  "PonatiLink-1-24" = "PonatiLink-1",
  "PonatiLink-1-28" = "PonatiLink-1-PEG28",
  "PonatiLink-2-7" = "PonatiLink-2-PEG7",
  "PonatiLink-2-7-4" = "PonatiLink-2-PEG11",
  "PonatiLink-2-7-6" = "PonatiLink-2-PEG13",
  "PonatiLink-2-7-8" = "PonatiLink-2-PEG15",
  # "PonatiLink-2-7-10" = "PonatiLink-2-PEG17",
  "PonatiLink-2-7-10" = "PonatiLink-2",
  "PL-2-7-10" = "PonatiLink-2",
  "PonatiLink-2-7-12" = "PonatiLink-2-PEG19",
  "PonatiLink-2-7-16" = "PonatiLink-2-PEG23",
  "DasatiLink-4" = "DasatiLink-1-PEG6",
  "DasatiLink-3" = "DasatiLink-1-PEG8",
  "DasatiLink-2" = "DasatiLink-1-PEG10",
  # "DasatiLink-1" = "DasatiLink-1-PEG12",
  "DasatiLink-1" = "DasatiLink-1"
)
color_map_treatments <- c(
  "ponatinib" = "#A40000", # apparent darkred from chimera
  "asciminib" = "#FF5300", # apparent orangered from chimera
  "ponatinib + asciminib" = "gold2",
  "imatinib" = "palegreen2",
  "dasatinib" = "mediumpurple",
  "dasatinib + asciminib" = "plum2",
  "RMC-5552" = "olivedrab3",
  "DasatiLink-1" = "black",
  "DasatiLink-2" = "green4",
  "DasatiLink-3" = "palegreen3",
  "DasatiLink-4" = "palegreen",
  "PonatiLink-1-12" = "plum",
  "PonatiLink-1-16" = "paleturquoise3",
  "PonatiLink-1-20" = "royalblue",
  "PonatiLink-1-24" = "forestgreen",
  "PonatiLink-1-28" = "olivedrab3",
  "PonatiLink-2-7" = "burlywood2",
  "PonatiLink-2-7-4" = "darkolivegreen1",
  "PonatiLink-2-7-6" = "palegreen2",
  "PonatiLink-2-7-8" = "forestgreen",
  "PonatiLink-2-7-10" = "skyblue3",
  "PL-2-7-10" = "skyblue3",
  "PonatiLink-2" = "skyblue3",
  "PonatiLink-2-7-12" = "blue2",
  "PonatiLink-2-7-16" = "purple3"
)
shape_map_treatments <- c(
  "ponatinib" = "circle",
  "asciminib" = "triangle",
  "ponatinib + asciminib" = "square",
  "imatinib" = "triangle down open",
  "dasatinib" = "triangle open",
  "dasatinib + asciminib" = "triangle down open",
  "RMC-5552" = "diamond",
  "DasatiLink-4" = "plus",
  "DasatiLink-3" = "diamond plus",
  "DasatiLink-2" = "circle cross",
  "DasatiLink-1" = "square cross",
  "PonatiLink-1-12" = "plus",
  "PonatiLink-1-16" = "square cross",
  "PonatiLink-1-20" = "diamond plus",
  "PonatiLink-1-24" = "circle cross",
  "PonatiLink-1-28" = "asterisk",
  "PonatiLink-2-7" = "plus",
  "PonatiLink-2-7-4" = "circle open",
  "PonatiLink-2-7-6" = "triangle open",
  "PonatiLink-2-7-8" = "triangle down open",
  "PonatiLink-2-7-10" = "square open", # consistent across plots
  "PL-2-7-10" = "square open", # consistent across plots
  "PonatiLink-2" = "square open", # consistent across plots
  "PonatiLink-2-7-12" = "diamond open",
  "PonatiLink-2-7-16" = "cross"
)
# specify targets to rename in plots. "old_name" = "new_name"
display_names_targets <- c(
  "K562 T315I" = "K562 BE-T315I",
  "K562 CRISPRi nt" = "K562 CRISPRi non-targeting",
  "K562 pUltra control" = "control",
  # "K562 pUltra BCR-ABL1 wt" = "K562 pU BCR::ABL1 wt",
  "K562 pUltra BCR-ABL1 wt" = "wt",
  "K562 pUltra BCR-ABL1 E255V" = "E255V",
  "K562 pUltra BCR-ABL1 T315I" = "T315I",
  "K562 pUltra BCR-ABL1 E255V T315I" = "E255V T315I",
  "K562 pUltra BCR-ABL1 T315M" = "T315M",
  "K562 pUltra BCR-ABL1 V468F" = "V468F",
  "K562 pUltra BCR-ABL1 E255V V468F" = "E255V V468F",
  "K562 pUltra BCR-ABL1 T315I V468F" = "T315I V468F",
  "K562 pUltra BCR-ABL1 E255V T315I V468F" = "E255V T315I V468F",
  "K562 pUltra BCR-ABL1 T315M V468F" = "T315M V468F"
)
color_map_targets <- c(
  # ABL1 cell lines-----------------------------------------------
  "K562 wt" = "forestgreen",
  "K562 T315I" = "royalblue",
  "K562 pUltra control" = "gray",
  "K562 pUltra BCR-ABL1 wt" = "forestgreen",
  "K562 pUltra BCR-ABL1 E255V" = "palegreen2",
  "K562 pUltra BCR-ABL1 T315I" = "skyblue2",
  "K562 pUltra BCR-ABL1 E255V T315I" = "royalblue",
  "K562 pUltra BCR-ABL1 T315M" = "purple",
  "K562 pUltra BCR-ABL1 V468F" = "gold",
  "K562 pUltra BCR-ABL1 E255V V468F" = "orchid1",
  "K562 pUltra BCR-ABL1 T315I V468F" = "plum1",
  "K562 pUltra BCR-ABL1 E255V T315I V468F" = "darkorange",
  "K562 pUltra BCR-ABL1 T315M V468F" = "firebrick3"
)

shape_map_targets <- c(
  # ABL1 cell lines-----------------------------------------------
  "K562 wt" = "circle",
  "K562 T315I" = "square",
  "K562 pUltra control" = "circle open",
  "K562 pUltra BCR-ABL1 wt" = "circle",
  "K562 pUltra BCR-ABL1 E255V" = "triangle",
  "K562 pUltra BCR-ABL1 T315I" = "square",
  "K562 pUltra BCR-ABL1 E255V T315I" = "diamond",
  "K562 pUltra BCR-ABL1 T315M" = "cross",
  "K562 pUltra BCR-ABL1 V468F" = "circle cross",
  "K562 pUltra BCR-ABL1 E255V V468F" = "triangle open",
  "K562 pUltra BCR-ABL1 T315I V468F" = "square open",
  "K562 pUltra BCR-ABL1 E255V T315I V468F" = "diamond open",
  "K562 pUltra BCR-ABL1 T315M V468F" = "asterisk"
)