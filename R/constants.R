
#' order for matching
#' @export
h3_clade_grouping_heirarchy <- c(
  "3C.2A1",
  "3C.2A1A",
  "3C.2A1B",
  "3C.2A2",
  "3C.2A",
  "3C.3A",
  "3C.3"
)

#' order for plotting
#' @export
h3_clade_order <- c(
  "3C.3",
  "3C.3A",
  "3C.2A",
  "3C.2A1",
  "3C.2A1A",
  "3C.2A1B",
  "3C.2A2",
  "Other",
  "Clade missing"
)


#' @export
h3_clusters <- function(){
  c(
    "HK68",
    "EN72",
    "VI75",
    "TX77",
    "BK79",
    "SI87",
    "BE89",
    "BE92",
    "WU95",
    "SY97",
    "FU02",
    "CA04",
    "WI05",
    "PE09",
    "SW13",
    "HK14"
  )
}

#' @export
h3_cluster_colors <- c(
  "HK68" = "#a208bd",
  "EN72" = "#00ffe1",
  "VI75" = "#f9d004",
  "TX77" = "#ab4c00",
  "BK79" = "#00ff00",
  "SI87" = "#0000ff",
  "BE89" = "#ff0000",
  "BE92" = "#f894f8",
  "WU95" = "#37802b",
  "SY97" = "#00afff",
  "FU02" = "#ffd700",
  "CA04" = "green",
  "WI05" = "blue",
  "PE09" = "purple",
  "SW13" = "grey80",
  "HK14" = "grey80"
)


#' @export
h3_clade_colors <- c(
  "3C.3" = '#7570b3',
  "3C.3A" = '#e7298a',
  "3C.2A" = '#0B691B',
  "3C.2A1" = '#7D2B07',
  "3C.2A1A" = '#E2D11F',
  "3C.2A1B" = '#F4B329',
  "3C.2A2" = '#1FD16E',
  "Other" = 'black',
  "Clade missing" = 'grey50'
)



#' @export
h3_cluster_subs <- list(
  HK68 = c("155T"),
  EN72 = c("155Y"),
  VI75 = c("189K"),
  TX77 = c("158E", "193N"),
  BK79 = c("156E"),
  SI87 = c("155H", "159Y", "189R"),
  BE89 = c("145K"),
  BE92 = c("156K"),
  WU95 = c("145K"),
  SY97 = c("156Q", "158K"),
  FU02 = c("156H"),
  CA04 = c("145N"),
  WI05 = c("193F"),
  PE09 = c("158N", "189K")
)

#' @export
clusters.subsMatrix <- function(){
  output <- do.call(rbind, lapply(h3_cluster_subs, function(x){
    position <- substr(x, 1, 3)
    aa       <- substr(x, 4, 4)
    result   <- rep("", 7)
    result[match(position, b7_positions)] <- aa
    result
  }))
  colnames(output) <- b7_positions
  output
}


#' @export
clusters.colour <- function(clusters){
  h3_cluster_colors[as.character(clusters)]
}

#' @export
clusters.number <- function(clusters){
  clusters[clusters == "BR07"] <- "WI05"
  match(
    clusters,
    h3_clusters()
  )
}

#' @export
clusters.toFactor <- function(clusters, excluded_clusters = NULL){

  cluster_levels <- h3_clusters()
  factor(clusters, levels = cluster_levels[!cluster_levels %in% excluded_clusters])

}


#' @export
b7_positions <- c(145, 155, 156, 158, 159, 189, 193)

#' @export
accessory_residues <- c(133, 143)

#' @export
h3_2004_variable_positions <- c(
  1, 2, 3, 5, 6, 7, 8, 10, 18, 21, 25, 29, 31, 33, 34, 45, 47, 49, 50, 53, 54, 55, 57, 58, 59,
  62, 63, 67, 75, 78, 79, 80, 81, 82, 83, 88, 92, 94, 96, 101, 103, 106, 112, 114, 115, 121, 122,
  124, 126, 128, 129, 131, 132, 133, 135, 137, 138, 139, 140, 142, 143, 144, 145, 146, 148, 155,
  156, 157, 158, 159, 160, 163, 164, 165, 167, 169, 171, 172, 173, 174, 175, 176, 182, 185, 186,
  188, 189, 190, 192, 193, 194, 196, 197, 198, 199, 201, 202, 203, 204, 205, 207, 208, 209, 213,
  214, 216, 217, 219, 220, 222, 225, 226, 227, 229, 230, 231, 233, 235, 236, 242, 244, 246, 247,
  248, 259, 260, 261, 262, 264, 267, 268, 269, 271, 272, 273, 275, 276, 278, 279, 280, 285, 286,
  289, 291, 292, 299, 304, 307, 309, 310, 312, 313, 323, 325, 326
)

siteA <- c(122, 124, 126, 130, 131, 132, 133, 135, 138, 142, 145)
siteB <- c(128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 190, 193, 194, 196, 197, 198)
siteC <- c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312)
siteD <- c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248)
RBS <- c(98, 134, 135, 136, 137, 138, 153, 183, 190, 194, 224, 225, 226, 227, 228)

siteA_Annette <- c(siteA, 137, 140, 143, 144, 146, 150, 152, 168)
siteB_Annette <- c(siteB, 187, 188, 189, 192)
siteC_Annette <- c(siteC, 276)
siteD_Annette <- siteD
siteE_Annette <- c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)

#' @export
h3_sites <- list(
  RBS = RBS,
  A = siteA,
  B = siteB,
  C = siteC,
  D = siteD,
  A_Annette = siteA_Annette,
  B_Annette = siteB_Annette,
  C_Annette = siteC_Annette,
  D_Annette = siteD_Annette,
  E_Annette = siteE_Annette
)



#' @export
aa_characteristics <- tibble::tribble(
  ~aa, ~character,
  "R", "Positive",
  "K", "Positive",
  "H", "Positive",
  "D", "Negative",
  "E", "Negative",
  "S", "Polar uncharged",
  "T", "Polar uncharged",
  "N", "Polar uncharged",
  "Q", "Polar uncharged",
  "C", "Special",
  "G", "Special",
  "P", "Special",
  "A", "Hydrophobic",
  "V", "Hydrophobic",
  "I", "Hydrophobic",
  "L", "Hydrophobic",
  "M", "Hydrophobic",
  "F", "Hydrophobic",
  "Y", "Hydrophobic",
  "W", "Hydrophobic"
)





#' @export
aa.values <- function(){
  c(
    c("G","A","S","T"),
    c("C","V","I","L","P","F","Y","M","W"),
    c("N","Q","H"),
    c("D","E"),
    c("K","R"),
    c("*", "X")
  )
}

#' @export
aa.colors <- function(){
  aa     <- aa.values()
  aacols <- as.vector(coloraa(rbind(aa)))
  names(aacols) <- aa
  aacols
}

#' @export
coloraa <- function(x, alpha = 0.5, simplify = T){
  if(is.vector(x)) x <- rbind(x)
  if(ncol(x) == 0) return(x)
  aacolors   <- x
  aacolors[] <- ""
  aacolors[x %in% c("G","A","S","T")                     ] <- "orange"
  aacolors[x %in% c("C","V","I","L","P","F","Y","M","W") ] <- "green"
  aacolors[x %in% c("N","Q","H","Z")                     ] <- "magenta"
  aacolors[x %in% c("D","E","B")                         ] <- "red"
  aacolors[x %in% c("K","R")                             ] <- "blue"
  aacolors[x %in% c("*", "X", "-", ".", "")              ] <- "grey50"
  if(sum(aacolors == "") > 0){
    stop(sprintf("aas '%s' not recognised.", paste(unique(x[aacolors == ""]), collapse = ", ")))
  }
  aacolors[] <- adjustcolor(aacolors, alpha.f = alpha)
  if (simplify){
    if (length(aacolors) == 1) return(aacolors[[1]])
  }

    aacolors
}

#'@export
aacolors_al <- c(
  "D" = "#F7CF08",
  "E" = "#CD8D00",
  "R" = "#B2A1FF",
  "H" = "#7571BC",
  "K" = "#440FAE",
  "N" = "#18E0F7",
  "Q" = "#0080FF",
  "S" = "#0F52BA",
  "T" = "navy",
  "C" = "#00F000",
  "G" = "#4CBB17",
  "P" = "#1A4F30",
  "A" = "#EB1919",
  "I" = "#F93BBB",
  "L" = "#FF2976",
  "M" = "#FF5A3D",
  "F" = "#C34822",
  "W" = "#A30000",
  "Y" = "#941651",
  "V" = "#681226"
)

#' All known two-letter place name abbreviations
#' @export
place_abvs <- c(
  "-A" = "SAUDI-ARABIA",
  "-B" = "SUPHAN-BURI",
  "-D" = "NORDRHEIN-WESTFALEN",
  "-G" = "FRENCH-GUIANA",
  "-H" = "SENDI-H",
  "-I" = "RHODE-ISLAND",
  "-K" = "SA-KAEO",
  "-L" = "SRI-LANKA",
  "-M" = "PUERTO-MONTT",
  "-N" = "SAKON-NAKHON",
  "-O" = "DISTRICT-OF-COLUMBIA",
  "-P" = "ST-PETERSBURG",
  "-S" = "SAN-SEBASTIAN",
  "-T" = "SAN-ANTONIO",
  "-V" = "MECKLENBURG-VORPO.",
  "AA" = "ANNARBOR",
  "AB" = "ALABAMA",
  "AC" = "AICHI",
  "AD" = "ALBANIA",
  "AE" = "ALGERIA",
  "AF" = "SOUTH-AFRICA",
  "AG" = "ARGENTINA",
  "AH" = "AUSTRIA",
  "AI" = "ASTURIAS",
  "AJ" = "SAKAI",
  "AK" = "AKITA",
  "AL" = "AUSTRALIA",
  "AM" = "AMSTERDAM",
  "AN" = "ANHUI",
  "AO" = "AOMORI",
  "AP" = "CHIANGMAI",
  "AQ" = "ANG-THONG",
  "AR" = "ARIZONA",
  "AS" = "ALICESPRINGS",
  "AT" = "ATLANTA",
  "AU" = "AUCKLAND",
  "AV" = "SANTIAGAO",
  "AX" = "ALASKA",
  "AY" = "AYTTHAYA",
  "AZ" = "AYATTHAYA",
  "BA" = "BANGKOK",
  "BB" = "BILBAO",
  "BC" = "BUCHAREST",
  "BD" = "BADEN-WURTTEMBERG",
  "BE" = "BEIJING",
  "BF" = "BALEARES",
  "BG" = "BELGIUM",
  "BH" = "BELGRADE",
  "BI" = "BILTHOVEN",
  "BJ" = "BEIIJNG",
  "BK" = "BANGKOKI",
  "BL" = "BELEM",
  "BM" = "BREMEN",
  "BN" = "BERLIN",
  "BO" = "BARCELONA",
  "BP" = "BURIRUM",
  "BQ" = "PATTANI",
  "BR" = "BRISBANE",
  "BS" = "BRASOV",
  "BT" = "BAGKOK",
  "BU" = "BUSAN",
  "BV" = "PHARE",
  "BW" = "PHICHIT",
  "BX" = "BURIRAM",
  "BY" = "BAYERN",
  "BZ" = "BRAZIL",
  "CA" = "CANBERRA",
  "CB" = "CHIBA",
  "CC" = "CHRISTCHURCH",
  "CD" = "COLLINDALE",
  "CE" = "CAEN",
  "CF" = "CALIFORNIA",
  "CG" = "CHUNGNAM",
  "CH" = "CHINA",
  "CI" = "CAEN",
  "CJ" = "CHEJU",
  "CK" = "CHEONBUK",
  "CL" = "CHILE",
  "CM" = "CHIANGMAI",
  "CN" = "CONNECTICUT",
  "CO" = "COLORADO",
  "CP" = "CLUJ",
  "CQ" = "CANARIAS",
  "CR" = "CHIANGRAI",
  "CS" = "CALARASI",
  "CT" = "CHITA",
  "CU" = "CHANTABURI",
  "CV" = "CARTAGENA",
  "CW" = "CHANGWON",
  "CZ" = "CZECHOSLOVAKIA",
  "DA" = "DAEGU",
  "DB" = "NY",
  "DC" = "CA",
  "DD" = "DUNDIN",
  "DE" = "DELFT",
  "DF" = "KYUNGGI",
  "DG" = "PAU",
  "DH" = "MILAN",
  "DI" = "SALAMANCA",
  "DJ" = "DAEJEON",
  "DK" = "DAKAR",
  "DL" = "CASABLANCA",
  "DM" = "MEKNES",
  "DN" = "DENMARK",
  "DO" = "MARRAKECH",
  "DP" = "AGADIR",
  "DQ" = "ATHENS",
  "DR" = "SANTANDER",
  "DS" = "DESGENETTES",
  "DU" = "DUNEDIN",
  "DV" = "DEVA",
  "DW" = "DARWIN",
  "DX" = "TENNESEE",
  "DY" = "BRANDYS",
  "DZ" = "TRANG",
  "EA" = "SIENA",
  "EB" = "ANGTHONG",
  "EC" = "ECUADOR",
  "ED" = "NIEDERSACHSEN",
  "EE" = "NEPAL",
  "EF" = "KALASIN",
  "EG" = "EGYPT",
  "EH" = "EHIME",
  "EI" = "EINDHOVEN",
  "EJ" = "PHITSANULOK",
  "EK" = "EKATERINBURG",
  "EL" = "EL-SALVADOR",
  "EM" = "PRACHUAPKHIRIKKAN",
  "EN" = "ENGLAND",
  "EO" = "CHACHOENGSAO",
  "EP" = "NEW-HAMPSHIRE",
  "EQ" = "EQUADOR",
  "ER" = "TOMSK",
  "ES" = "ENSCHEDE",
  "ET" = "EXTREMAD",
  "EU" = "EUSKADI",
  "EV" = "PONTEVEORA",
  "EX" = "EXTREMADURA",
  "EY" = "KYOTO",
  "EZ" = "EKATEZINBURG",
  "FA" = "FATICK",
  "FB" = "FUKUI",
  "FC" = "FUKUOKA-C",
  "FD" = "FES",
  "FE" = "FIRENZE",
  "FF" = "ANNECY",
  "FG" = "MACU",
  "FH" = "FOSHAN",
  "FI" = "FINLAND",
  "FJ" = "FIJI",
  "FK" = "FUOKA",
  "FL" = "FLORIDA",
  "FM" = "SUKHBAATAR",
  "FN" = "SAINSHAND",
  "FO" = "FLORENCE",
  "FP" = "CAMBODIA",
  "FQ" = "HOLLAND",
  "FR" = "FRANCE",
  "FS" = "FUKUSHIMA",
  "FT" = "DUBAI",
  "FU" = "FUJIAN",
  "FV" = "GANSUCHENGGUAN",
  "FW" = "HONDURAS",
  "FX" = "KALININGRAD",
  "FY" = "ASTRAKHAN",
  "FZ" = "FUZHOU",
  "GA" = "GUAM",
  "GB" = "GENOA",
  "GC" = "GREECE",
  "GD" = "GUANGDONG",
  "GE" = "GENEVA",
  "GF" = "GUILDFORD",
  "GG" = "GEORGIA",
  "GH" = "GIRONA",
  "GI" = "GIFU",
  "GJ" = "GUADALAJARA",
  "GK" = "GOTEBURG",
  "GL" = "GUALDALAGARA",
  "GM" = "GUMMA",
  "GN" = "GRANADA",
  "GO" = "GOTENBORG",
  "GP" = "GIFU-C",
  "GQ" = "GALICIA",
  "GR" = "GRONINGEN",
  "GS" = "GANSU",
  "GT" = "GOTEBORG",
  "GU" = "GUIZHOU",
  "GV" = "SEGOVIA",
  "GW" = "GUNMA",
  "GX" = "GUANGXI",
  "GY" = "GERMANY",
  "GZ" = "GUANGZHOU",
  "HA" = "HAWAII",
  "HB" = "HARBIN",
  "HC" = "SACHSEN",
  "HD" = "HOKKAIDO",
  "HE" = "HEBEI",
  "HF" = "HIROSHIMA",
  "HG" = "HUNGARY",
  "HH" = "SHIZUOKA",
  "HI" = "HUBEI",
  "HJ" = "HYOGO",
  "HK" = "HONGKONG",
  "HL" = "HIME",
  "HM" = "HENAN",
  "HN" = "HUNAN",
  "HO" = "HOUSTON",
  "HP" = "HONG-KONGK",
  "HQ" = "HAINAN",
  "HS" = "HIROSHIMA-C",
  "HT" = "HAMAMATU-C",
  "HU" = "HUNTINGTON",
  "HV" = "HANNOVER",
  "HW" = "NAKHON-SAWAN",
  "HX" = "SHAANXI",
  "HY" = "CHANTHABURI",
  "HZ" = "SHIZUOKA-C",
  "IA" = "IASI",
  "IB" = "IBARAKI",
  "IC" = "INCHEON",
  "ID" = "INDONESIA",
  "IE" = "ICELAND",
  "IF" = "IRAN",
  "IH" = "IDAHO",
  "II" = "INDIA",
  "IJ" = "BEIJINGXUANWU",
  "IK" = "ISHIKAWA",
  "IL" = "ILLINOIS",
  "IM" = "SHIMANE",
  "IN" = "INDIANA",
  "IO" = "SOLOMON-ISLANDS",
  "IP" = "MURMANSK",
  "IQ" = "STPETERSBURG",
  "IR" = "IRELAND",
  "IS" = "ISRAEL",
  "IT" = "ITALY",
  "IU" = "ISTANBUL",
  "IV" = "INVERNESS",
  "IW" = "IWATE",
  "IX" = "BEIJINXUANWU",
  "IY" = "CHAIYAPUM",
  "IZ" = "HAMBURG",
  "JA" = "JIANGXI",
  "JB" = "BRANDENBURG",
  "JC" = "BURSA",
  "JD" = "EDIRNE",
  "JE" = "ANTALYA",
  "JF" = "KHOVD",
  "JG" = "NEWYORK",
  "JH" = "JHB",
  "JI" = "JILIN",
  "JJ" = "JEJU",
  "JK" = "NOVGOROD",
  "JL" = "TURKEY",
  "JM" = "SUDAN",
  "JN" = "AFGHANISTAN",
  "JO" = "JOHANNESBURG",
  "JP" = "JAPAN",
  "JQ" = "NIGERIA",
  "JR" = "KENTUCKY",
  "JS" = "JIANGSU",
  "JT" = "TALCAHUANO",
  "JU" = "UTAH",
  "JV" = "MASSACHUSETTS",
  "JW" = "KANSAS",
  "JX" = "JIANGXIDONGHU",
  "JY" = "MONTANA",
  "JZ" = "IOWA",
  "KA" = "KASAULI",
  "KB" = "KAGAWA",
  "KC" = "KOCHI",
  "KD" = "KANAGAWA",
  "KE" = "KUMAMOTO-C",
  "KF" = "KUM",
  "KG" = "KAGOSHIMA",
  "KH" = "KHABAROVSK",
  "KI" = "KYONGGI",
  "KJ" = "KWANGJU",
  "KK" = "KYONGBUK",
  "KL" = "KHON-KAEN",
  "KM" = "KUMAMOTO",
  "KN" = "KANGWON",
  "KO" = "KOBE",
  "KP" = "NAKHON-PATHOM",
  "KQ" = "KAMPHAENG-PHET",
  "KR" = "KOREA",
  "KS" = "KAWASAKI",
  "KT" = "KYOTO-C",
  "KU" = "KITAKYUSYU",
  "KV" = "KANCHANABURI",
  "KW" = "KWANGJIN",
  "KX" = "KYUNGNAM",
  "KY" = "KITAKYUSHU",
  "KZ" = "KYONGNAM",
  "LA" = "LAUSANNE",
  "LB" = "LOPBURI",
  "LC" = "LOP-BURI",
  "LD" = "LYON-TRS",
  "LE" = "LENINGRAD",
  "LF" = "LOBBURI",
  "LG" = "LIAONING",
  "LH" = "LYON-CHU",
  "LI" = "LINNKOPING",
  "LJ" = "LAMPANG",
  "LK" = "LAMOANG",
  "LL" = "CASTILLA",
  "LM" = "LIMOGES",
  "LN" = "LEON",
  "LO" = "LOSANGELES",
  "LP" = "LIPETZK",
  "LQ" = "LOEI",
  "LR" = "LA-REUNION",
  "LS" = "LISBON",
  "LT" = "LATVIA",
  "LU" = "LOUISIANA",
  "LV" = "SLOVENIA",
  "LW" = "SCHLESWIG-HOLSTEIN",
  "LX" = "LEICESTERSHIRE",
  "LY" = "LYON",
  "LZ" = "LINCOLN",
  "MA" = "MADRID",
  "MB" = "MOROCCO",
  "MC" = "MICHIGAN",
  "MD" = "MAD",
  "ME" = "MEMPHIS",
  "MF" = "CLERMONTFERRAND",
  "MG" = "MIYAGI",
  "MH" = "MAE-HONG-SORN",
  "MI" = "MISSISSIPPI",
  "MJ" = "MIE",
  "MK" = "MIYAZAKI",
  "ML" = "MALY",
  "MM" = "MALMO",
  "MN" = "MAE-HONG-SON",
  "MO" = "MISSOURI",
  "MP" = "MONTPELLIER",
  "MQ" = "MACAU",
  "MR" = "MADAGASCAR",
  "MS" = "MINNESOTA",
  "MT" = "SAMUT-SAKHON",
  "MU" = "SAMUT-PRAKAN",
  "MV" = "MONGOLIA",
  "MW" = "MOSCOW",
  "MX" = "MEXICO",
  "MY" = "MAYOCLINIC",
  "MZ" = "MASSACHUSETS",
  "NA" = "NANCHANG",
  "NB" = "NEBRASKA",
  "NC" = "NEWCASTLE",
  "ND" = "NORTH-DAKOTA",
  "NE" = "NICE",
  "NF" = "NRW",
  "NG" = "NIIGATA",
  "NH" = "NARA",
  "NI" = "NIJMEGEN",
  "NJ" = "NEW-JERSEY",
  "NK" = "NAGANO",
  "NL" = "NETHERLANDS",
  "NM" = "NAGOYA",
  "NN" = "NINGBO",
  "NO" = "NORTH-CAROLINA",
  "NP" = "NIGATA",
  "NQ" = "NAPAL",
  "NR" = "NOVOSIBIRSK",
  "NS" = "NAGASAKI",
  "NT" = "NIIGATA-C",
  "NU" = "NONTHABURI",
  "NV" = "NEVADA",
  "NW" = "NEW-CALEDONIA",
  "NX" = "NINGXIA",
  "NY" = "NEW-YORK",
  "NZ" = "NAKHON-RATCHASIMA",
  "OA" = "OSAKA-C",
  "OB" = "NOVY-BYDZOV",
  "OC" = "OTAGA",
  "OD" = "NORDRHEIN",
  "OE" = "ORENSE",
  "OF" = "OREL",
  "OG" = "OREGON",
  "OH" = "OHIO",
  "OI" = "OITA",
  "OJ" = "BOLIVIA",
  "OK" = "OKLAHOMA",
  "OL" = "CHOIBALSAN",
  "OM" = "OMSK",
  "ON" = "OKINAWA",
  "OO" = "NAKONNAYOK",
  "OP" = "GUADELOUPE",
  "OQ" = "CHONGQING",
  "OR" = "ORADEA",
  "OS" = "OSLO",
  "OT" = "OTAGO",
  "OU" = "OUJDA",
  "OV" = "OVIEDO",
  "OW" = "NORWAY",
  "OX" = "MAINE",
  "OY" = "OKAYAMA",
  "OZ" = "RUSSIA",
  "PA" = "PARIS",
  "PB" = "PARMA",
  "PC" = "PORT-CHALMERS",
  "PD" = "POLAND",
  "PE" = "PERTH",
  "PF" = "PHATTHALUNG",
  "PG" = "PERUGIA",
  "PH" = "PHILIPPINES",
  "PI" = "POITIERS",
  "PJ" = "PRAJIANBURI",
  "PK" = "PHUKET",
  "PL" = "PILSEN",
  "PM" = "PANAMA",
  "PN" = "PUSAN",
  "PO" = "PUERTO-RICO",
  "PP" = "PRACHUAP-KHIRI-KHAN",
  "PQ" = "PATHUMTHANI",
  "PR" = "PRAGUE",
  "PS" = "PENNSYLVANIA",
  "PT" = "PATHUM-THANI",
  "PU" = "PERU",
  "PV" = "PHILLIPINES",
  "PW" = "PRACHINBURI",
  "PX" = "PHETCHABUN",
  "PY" = "PARAGUAY",
  "PZ" = "PRABCHINBURI",
  "QA" = "QUANZHOU",
  "QB" = "BANGLADESH",
  "QC" = "CORDOBA",
  "QD" = "ARKANSAS",
  "QE" = "MARTINIQUE",
  "QF" = "COLOMBIA",
  "QG" = "JAMAICA",
  "QH" = "POL",
  "QI" = "KIEV",
  "QJ" = "ODESSA",
  "QK" = "PALENCIA",
  "QL" = "NOUACKCHOTT",
  "QM" = "HESSEN",
  "QN" = "GERONA",
  "QO" = "MAURITIUS",
  "QP" = "FUKUOKAC",
  "QQ" = "CONGQING",
  "QR" = "GIFUC",
  "QS" = "HIROSHIMAC",
  "QT" = "QINGDAO",
  "QU" = "QUEENSLAND",
  "QV" = "KUMAMOTOC",
  "QW" = "NIIGATAC",
  "QX" = "SENDAIH",
  "QY" = "MARSEILLE",
  "QZ" = "BORDEAUX",
  "RA" = "ROME",
  "RB" = "ROI-ET",
  "RC" = "PRACHINMURI",
  "RD" = "ROTTERDAM",
  "RE" = "REUNION",
  "RF" = "RABAT",
  "RG" = "ARAGON",
  "RH" = "RHEINLAND-PFALZ",
  "RI" = "SORIA",
  "RJ" = "RIO-DE-JANEIRO",
  "RK" = "SAMUTPRAKHAN",
  "RL" = "BRATISLAVA",
  "RM" = "ROMANIA",
  "RN" = "SURIN",
  "RO" = "ROMA",
  "RP" = "PRACHUAPKHIRIKHAN",
  "RQ" = "ROSTOV-ON-DON",
  "RR" = "SARABURI",
  "RS" = "BURGOS",
  "RT" = "RATCHABURI",
  "RU" = "RU",
  "RV" = "ARVAIKHEER",
  "RW" = "NARATHIWAT",
  "RX" = "ROSTOVDON",
  "RY" = "KRASNOYARSK",
  "RZ" = "TRABZON",
  "SA" = "SOUTH-AUSTRALIA",
  "SB" = "SAGA",
  "SC" = "SOUTH-CAROLINA",
  "SD" = "SHANGDONG",
  "SE" = "SENDAI",
  "SF" = "SOFIA",
  "SG" = "SHIGA",
  "SH" = "SHANGHAI",
  "SI" = "SICHUAN",
  "SJ" = "SPAIN",
  "SK" = "ST.-ETIENNE",
  "SL" = "SCOTLAND",
  "SM" = "SAMARA",
  "SN" = "ST-ETIENNE",
  "SO" = "SOUTH-DAKOTA",
  "SP" = "SINGAPORE",
  "SQ" = "SACHSEN-ANHALT",
  "SR" = "SAPPORO",
  "SS" = "ST.-PETERSBURG",
  "ST" = "STOCKHOLM",
  "SU" = "SEOUL",
  "SV" = "SHANTOU",
  "SW" = "SW",
  "SX" = "SOPHIA",
  "SY" = "SYDNEY",
  "SZ" = "SANTIAGO",
  "TA" = "TAIWAN",
  "TB" = "THURINGEN",
  "TC" = "TOCHIGI",
  "TD" = "TRENTO",
  "TE" = "TEXAS",
  "TF" = "TARRAGONA",
  "TG" = "TONGA",
  "TH" = "THESSALONIKI",
  "TI" = "TILBURG",
  "TJ" = "TIANJIN",
  "TK" = "TAK",
  "TL" = "THAILAND",
  "TM" = "TASMANIA",
  "TN" = "TOULON",
  "TO" = "TOULOUSE",
  "TP" = "TOYAMA",
  "TQ" = "TOKUSHIMA",
  "TR" = "TEHRAN",
  "TS" = "TRIESTE",
  "TT" = "TOTTORI",
  "TU" = "TULA",
  "TV" = "TOWNSVILLE",
  "TW" = "TENNESSEE",
  "TX" = "TX",
  "TY" = "TOKYO",
  "TZ" = "TANGER",
  "UA" = "U.K.",
  "UB" = "UTTARADIT",
  "UC" = "UD",
  "UD" = "UDORN",
  "UE" = "UBON-RATCHATHANI",
  "UF" = "UBONRATCHATHANI",
  "UG" = "URUAGUAY",
  "UH" = "SUPHANBURI",
  "UI" = "UNITED-KINGDOM",
  "UJ" = "UTHAI-THANI",
  "UK" = "UK",
  "UL" = "ULSAN",
  "UM" = "UMEA",
  "UN" = "UNITEDKINGDOM",
  "UO" = "SUKHOTHAI",
  "UP" = "SAMUTPRAKAN",
  "UQ" = "UMEA",
  "UR" = "URUGUAY",
  "US" = "USSR/RUSSIA",
  "UT" = "UTRECHT",
  "UU" = "ULAN-UDE",
  "UV" = "SUCEAVA",
  "UW" = "ULAANBAATAR",
  "UX" = "BENELUX",
  "UY" = "GUYANE",
  "UZ" = "UKRAINE",
  "VA" = "VALENCIA",
  "VB" = "VINA-DEL-MAR",
  "VC" = "COSTARICA",
  "VD" = "VALLADOLID",
  "VE" = "VIETNAM",
  "VF" = "COOK-ISLAND",
  "VG" = "VIRGINIA",
  "VH" = "JEONBUK",
  "VI" = "VICTORIA",
  "VJ" = "GUATEMALA",
  "VK" = "SLOVAKIA",
  "VL" = "VOLGOGRAD",
  "VM" = "VLADIMIR",
  "VN" = "VIENNA",
  "VO" = "VORONEZH",
  "VP" = "KENYA",
  "VQ" = "SING",
  "VR" = "STAVROPOL",
  "VT" = "VERMONT",
  "VV" = "AVILA",
  "VY" = "IVORY-COAST",
  "VZ" = "VENEZUELA",
  "WA" = "WASHINGTON",
  "WB" = "WY--",
  "WE" = "WELLINGTON",
  "WG" = "GWANGJU",
  "WH" = "WUZHOU",
  "WK" = "WAIKATO",
  "WM" = "WAKAYAMA",
  "WN" = "WISCONSIN",
  "WO" = "GANGWON",
  "WR" = "WARSAW",
  "WS" = "WEST-VIRGINIA",
  "WU" = "WUHAN",
  "WV" = "WESTVIRGINIA",
  "WW" = "NEWMARKET",
  "WX" = "NEW-MEXICO",
  "WY" = "WYOMING",
  "WZ" = "SWITZERLAND",
  "XI" = "XINJIANGCHANGJI",
  "XJ" = "XINJIANG-HUTUBI",
  "XM" = "XIAMEN",
  "XN" = "XINJIANGHUTUBI",
  "XX" = "CX-ROUSSE",
  "YA" = "YAMAGA",
  "YB" = "YUNAN",
  "YD" = "BYDGOSZCZ",
  "YE" = "GYEONGNAM",
  "YG" = "YAMAGUCHI",
  "YI" = "YAMANESHI",
  "YK" = "YOKOSUKA",
  "YL" = "MARYLAND",
  "YM" = "YAMANASHI",
  "YN" = "KYUNGBUK",
  "YO" = "YOKOHAMA",
  "YR" = "YARYSLAVL",
  "YS" = "YAROSLAVL",
  "YT" = "YAMAGATA",
  "YU" = "YUNNAN",
  "YY" = "NAKHONNAYOK",
  "YZ" = "YUNGNAM",
  "ZA" = "ZAMBIA",
  "ZB" = "AZERBAIJAN",
  "ZE" = "ZAGREB",
  "ZG" = "ZARAGOSSA",
  "ZH" = "ZHEJIANG",
  "ZI" = "IZMIR",
  "ZL" = "ZLIN",
  "ZM" = "ZAMORA",
  "ZN" = "RYAZAN",
  "ZR" = "ZARAGOZA",
  "ZU" = "SHIZUOKAC",
  "ZX" = "CHUNGBUK",
  "ZY" = "SAITAMA",
  "OT" = "GOETEBORG"
)
