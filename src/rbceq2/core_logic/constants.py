from enum import Enum, auto

# Define version
VERSION = "2.4.0"
DB_VERSION = "2.4.0"


COMMON_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

HOM_REF_DUMMY_QUAL = ""
HOM_REF_DUMMY_QUAL += "./.:"  # GT (Genotype ie  0/1, 1|0, 0/0, 1/1.)
HOM_REF_DUMMY_QUAL += "1,29:"  # AD (Allelic Depth)
HOM_REF_DUMMY_QUAL += "30:"  # GQ (Genotype Quality)
HOM_REF_DUMMY_QUAL += "30:"  # DP (Read Depth):
HOM_REF_DUMMY_QUAL += "1"  # PS (Phase Set)

LOW_WEIGHT = 1_000


class AlleleState:
    CO = "co_existing"
    NORMAL = "pairs"
    RAW = "raw"
    # POS = "possible" #everything now a filter
    FILT = "filtered"


class PhenoType(Enum):
    alphanumeric = auto()
    numeric = auto()


TWO = 2


class BgName(Enum):
    CROM = auto()
    FY = auto()
    KN = auto()
    SC = auto()
    VEL = auto()
    DO = auto()
    ABCC4 = auto()
    JK = auto()
    FUT2 = auto()
    FUT3 = auto()
    LU = auto()
    ABCB6 = auto()
    GLOB = auto()
    AUG = auto()
    CO = auto()
    KEL = auto()
    GBGT1 = auto()
    XK = auto()
    ABO = auto()
    GYPA = auto()
    GYPB = auto()
    CD59 = auto()
    IN = auto()
    RAPH = auto()
    JMH = auto()
    ER = auto()
    DI = auto()
    SID = auto()
    CTL2 = auto()
    FUT1 = auto()
    KLF = auto()
    LW = auto()
    MAM = auto()
    OK = auto()
    GE = auto()
    KANNO = auto()
    A4GALT = auto()
    ABCG2 = auto()
    PIGG = auto()
    GCNT2 = auto()
    RHAG = auto()
    YT = auto()
    GIL = auto()
    GATA1 = auto()
    XG = auto()
    ABCC1 = auto()
    RHD = auto()
    RHCE = auto()
    C4A = auto()
    C4B = auto()
    ATP11C = auto()
    HPA1 = auto()
    HPA10 = auto()
    HPA11 = auto()
    HPA14 = auto()
    HPA16 = auto()
    HPA17 = auto()
    HPA19 = auto()
    HPA2 = auto()
    HPA20 = auto()
    HPA21 = auto()
    HPA22 = auto()
    HPA23 = auto()
    HPA24 = auto()
    HPA26 = auto()
    HPA27 = auto()
    HPA28 = auto()
    HPA29 = auto()
    HPA3 = auto()
    HPA30 = auto()
    HPA32 = auto()
    HPA33 = auto()
    HPA34 = auto()
    HPA35 = auto()
    HPA4 = auto()
    HPA6 = auto()
    HPA7 = auto()
    HPA8 = auto()
    HPA9 = auto()
    HPA12 = auto()
    HPA31 = auto()
    HPA13 = auto()
    HPA18 = auto()
    HPA25 = auto()
    HPA5 = auto()
    HPA15 = auto()
    CD99 = auto()

    @classmethod
    def from_string(cls, value: str):
        value = value.upper()
        if "KLF" in value:
            value = "KLF"
        try:
            return cls[value]
        except KeyError:
            raise ValueError(f"'{value}' is not a valid {cls.__name__}")


ANTITHETICAL = {
    PhenoType.numeric: {
        BgName.LU: {
            "1": ("2",),
            "2": ("1",),
            "6": ("9",),
            "9": ("6",),
            "8": ("14",),
            "14": ("8",),
            "18": ("19",),
            "19": ("18",),
        },
        BgName.JK: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.DI: {
            "1": ("2",),
            "2": ("1",),
            "3": ("4",),
            "4": ("3",),
            "9": ("22",),
            "22": ("9",),
            "11": ("12",),
            "12": ("11",),
            "15": ("16",),
            "16": ("15",),
            "17": ("18",),
            "18": ("17",),
        },
        BgName.ER: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.DO: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.SC: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.FY: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.KN: {
            "1": ("2",),
            "2": ("1",),
            "3": ("6",),
            "6": ("3",),
            "4": ("7",),
            "7": ("4",),
            "9": ("10",),
            "10": ("9",),
            "11": ("12",),
            "12": ("11",),
        },
        BgName.KEL: {
            "1": ("2",),
            "2": ("1",),
            "3": ("4", "21"),  # all 3 possible
            "4": ("3", "21"),
            "21": ("4", "3"),
            "6": ("7",),
            "7": ("6",),
            "11": ("17",),
            "17": ("11",),
            "14": ("24",),
            "24": ("14",),
            "31": ("38",),
            "38": ("31",),
            "37": ("39",),
            "39": ("37",),
            "40": ("41",),
            "41": ("40",),
        },
        # GPT Added blood groups
        BgName.GYPA: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.GYPB: {  # are these actaully antithetical?
            "3": ("4",),
            "4": ("3",),
        },
        BgName.YT: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.CO: {  # are these actaully antithetical?
            "1": ("2",),
            "2": ("1",),
        },
        BgName.LW: {
            "5": ("7",),
            "7": ("5",),
        },
        BgName.IN: {
            "1": ("2",),
            "2": ("1",),
        },
        BgName.CROM: {
            "2": ("3", "4"),
            "3": ("2", "4"),
            "4": ("2", "3"),
        },
        BgName.RHCE: {
            "2": ("4",),
            "4": ("2",),
            "3": ("5",),
            "5": ("3",),
            "8": ("9"),
            "9": ("8"),
            "26": ("55",),
            "55": ("26",),
            "32": ("46",),
            "46": ("32",),
            "43": ("58",),
            "58": ("43",),
            "48": ("57",),
            "57": ("48",),
        },
    },
    PhenoType.alphanumeric: {  # TODO RHCE - need to sort out ant names with Eileen
        BgName.LU: {
            "Lu(a)": ("Lu(b)",),
            "Lu(b)": ("Lu(a)",),
            "Lu6": ("Lu9",),
            "Lu9": ("Lu6",),
            "Lu8": ("Lu14",),
            "Lu14": ("Lu8",),
            "Au(a)": ("Au(b)",),
            "Au(b)": ("Au(a)",),
        },
        BgName.JK: {
            "Jk(a)": ("Jk(b)",),
            "Jk(b)": ("Jk(a)",),
        },
        BgName.DI: {
            "Di(a)": ("Di(b)",),
            "Di(b)": ("Di(a)",),
            "Wr(a)": ("Wr(b)",),
            "Wr(b)": ("Wr(a)",),
            "Wu": ("DISK",),
            "DISK": ("Wu",),
            "Moa": ("Hga",),
            "Hga": ("Moa",),
            "BOW": ("NFLD",),
            "NFLD": ("BOW",),
            "Jna": ("KREP",),
            "KREP": ("Jna",),
        },
        BgName.ER: {
            "Er(a)": ("Er(b)",),
            "ER(b)": ("Er(a)",),
        },
        BgName.DO: {
            "Do(a)": ("Do(b)",),
            "Do(b)": ("Do(a)",),
        },
        BgName.SC: {
            "Sc1": ("Sc2",),
            "Sc2": ("Sc1",),
        },
        BgName.FY: {
            "Fy(a)": ("Fy(b)",),
            "Fy(b)": ("Fy(a)",),
        },
        BgName.CROM: {
            "Tc(a)": ("Tc(b)", "Tc(c)"),  # only 2 possible
            "Tc(b)": ("Tc(a)", "Tc(c)"),
            "Tc(c)": ("Tc(a)", "Tc(b)"),
        },
        BgName.KN: {
            "Kn(a)": ("Kn(b)",),
            "Kn(b)": ("Kn(a)",),
            "McC(a)": ("McC(b)",),
            "McC(b)": ("McC(a)",),
            "Sl1": ("Vil",),
            "Vil": ("Sl1",),
            "KCAM": ("KDAS",),
            "KDAS": ("KCAM",),
            "DACY": ("YCAD",),
            "YCAD": ("DACY",),
        },
        BgName.KEL: {
            "K": ("k",),
            "k": ("K",),
            "Kp(a)": ("Kp(b)", "Kp(c)"),
            "Kp(b)": ("Kp(a)", "Kp(c)"),
            "Kp(c)": ("Kp(a)", "Kp(b)"),
            "Js(a)": ("Js(b)",),
            "Js(b)": ("Js(a)",),
            "K11": ("K17",),
            "K17": ("K11",),
            "K14": ("K24",),
            "K24": ("K14",),
            "KYO": ("KYOR",),
            "KYOR": ("KYO",),
            "KHUL": ("KEAL",),
            "KEAL": ("KHUL",),
            "KHIZ": ("KHOZ",),
            "KHOZ": ("KHIZ",),
        },
        # Added blood groups
        BgName.GYPA: {
            "M": ("N",),
            "N": ("M",),
        },
        BgName.GYPB: {
            "S": ("s",),
            "s": ("S",),
        },
        BgName.RHCE: {
            "C": ("c",),
            "c": ("C",),
            "E": ("e",),
            "e": ("E",),
            "Cw": ("Cx",),
            "Cx": ("Cw",),
            "c-like": ("LOCR",),
            "LOCR": ("c-like",),
            "Rh32": ("Sec",),
            "Sec": ("Rh32",),
            "Crawford": ("CELO",),
            "CELO": ("Crawford",),
            "JAL": ("CEST",),
            "CEST": ("JAL",),
        },
        BgName.YT: {
            "Yt(a)": ("Yt(b)",),
            "Yt(b)": ("Yt(a)",),
        },
        BgName.CO: {
            "Co(a)": ("Co(b)",),
            "Co(b)": ("Co(a)",),
        },
        BgName.LW: {
            "LW(a)": ("LW(b)",),
            "LW(b)": ("LW(a)",),
        },
        BgName.IN: {
            "In(a)": ("In(b)",),
            "In(b)": ("In(a)",),
        },
    },
}

LANE = {  # definition of lane has expanded to include any position that is implicity ref
    "chr1": {
        "1000": "no_ALT",  # unit test 159175354
        "159205564": "G_A",  # 38 FY
        "159175354": "G_A",  # 37
        "207331122": "G_C",  # 38 CROM
        "207504467": "G_C",  # 37
        "207609424": "G_A",  # 38 KN
        "207782769": "G_A",  # 37
        "207609571": "A_T",  # 38 KN
        "207782916": "A_T",  # 37
        "25317062": "T_C",  # 38 RHD
        "25643553": "T_C",  # 37
        "25390874": "C_G",  # 38 RHCE
        "25717365": "C_G",  # 37
        "25408711": "G_A",  # 38 RHCE
        "25735202": "G_A",  # 37
        "25420739": "G_C",  # 38 RHCE
        "25747230": "G_C",  # 37
        "42830851": "G_A",  # 38 SC
        "43296522": "G_A",  # 37
    },
    "chr11": {
        "35176644": "G_C",  # 38 IN
        "35198191": "G_C",
    },  # 37
    "chr12": {
        "14840505": "C_T",  # 38 DO
        "14993439": "C_T",
    },  # 37
    "chr16": {
        "88716069": "C_T",  # 38 ER
        "88782477": "C_T",
    },  # 37
    "chr17": {
        "44251253": "G_A",  # 38 DI
        "42328621": "G_A",  # 37
        "42453065": "A_C",  # HPAs have same coords
        "42453072": "G_T",
        "42453084": "C_T",
        "42453291": "C_G",
        "42453713": "C_A",
        "42455875": "G_A",
        "42457790": "C_T",
        "42462694": "T_G",
        "45351803": "C_T",
        "45360730": "T_C",
        "45360817": "G_A",
        "45360903": "C_T",
        "45361934": "A_C",
        "45361944": "C_T",
        "45361968": "A_G",
        "45363673": "C_T",
        "45369541": "C_G",
        "45369617": "A_G",
        "45369758": "G_A",
        "45369788": "G_A",
        "45376801": "G_T",
        "45376892": "TAAG_T",
        "45377872": "C_T",
        "45377890": "G_A",
        "45377906": "G_A",
        "45377914": "C_T",
        "45631953": "G_A",
        "4836381": "C_T",
    },
    "chr18": {
        "45739554": "G_A",  # 38 JK
        "43319519": "G_A",  # 37
    },
    "chr19": {
        "44812188": "G_A",  # 38 LU
        "45315445": "G_A",  # 37
        "5844526": "no_ALT",  # 38 FUT3
        "5844537": "no_ALT",  # 37
        "10287311": "A_G",  # 38 LW
        "10397987": "A_G",  # 37
        "5844638": "no_ALT",  # 38 FUT3
        "5844649": "no_ALT",  # 37
        "10631494": "no_ALT",  # 38 CTL2
        "10742170": "no_ALT",  # 37
    },
    "chr22": {
        "19711485": "G_A",
        "42717787": "C_A",  # 38 A4GALT
        "43113793": "C_A",  # 37
    },
    "chr3": {
        "161086379": "C_T",  # 38 GLOB
        "160804167": "C_T",  # 37
        "128780950": "C_T",
    },
    "chr4": {
        "144120555": "T_C",  # 38 GYPA
        "145041708": "T_C",  # 37
        "144120567": "A_G",  # 38 GYPA
        "145041720": "A_G",  # 37
        "143997559": "C_G",  # 38 GYPB
        "144918712": "C_G",  # 37
        "143999443": "G_A",  # 38 GYPB
        "144920596": "G_A",  # 37
        "144120554": "C_A",  # 38 GYPA
        "145041707": "C_A",  # 37
    },
    "chr5": {
        "52366090": "G_T",
        "52358757": "G_A",
        "52369001": "C_T",
        "52382870": "C_T",
    },
    "chr6": {
        "44232918": "G_A",  # 38 AUG
        "44200655": "G_A",  # 37
        "10586805": "no_ALT",  # 38 GCNT2
        "10587038": "no_ALT",  # 37
    },
    "chr7": {
        "142957921": "G_A",  # 38 KEL
        "142655008": "G_A",  # 37
        "100893176": "G_T",  # 38 YT
        "100490797": "G_T",
    },  # 37
    "chr8": {
        "74493432": "C_A",
    },
    "chr9": {
        "133257521": "T_TC",  # 38 ABO
        "136132908": "T_TC",  # 37
        "136133506": "no_ALT",  # 37 (only)
        "136135237": "no_ALT",  # 37 (only)
        "136136770": "no_ALT",  # 37 (only)
        "136135238": "no_ALT",
    },  # 37 (only)
}

RHD_ANT_MAP = {
    "RHD*01W.1": "Type1",
    "RHD*01W.1.1": "Type1.1",
    "RHD*01W.1.2": "Type1.2",
    "RHD*01W.1.3": "Type1.3",
    "RHD*01W.10": "Type10",
    "RHD*01W.10.1": "Type10.1",
    "RHD*01W.10.2": "Type10.2",
    "RHD*01W.100": "Type100",
    "RHD*01W.101": "Type101",
    "RHD*01W.102": "Type102",
    "RHD*01W.103": "Type103",
    "RHD*01W.104": "Type104",
    "RHD*01W.105": "Type105",
    "RHD*01W.106": "Type106",
    "RHD*01W.107": "Type107",
    "RHD*01W.108": "Type108",
    "RHD*01W.109": "Type109",
    "RHD*01W.110": "Type110",
    "RHD*01W.111": "Type111",
    "RHD*01W.112": "Type112",
    "RHD*01W.113": "Type113",
    "RHD*01W.114": "Type114",
    "RHD*01W.115": "Type115",
    "RHD*01W.116": "Type116",
    "RHD*01W.117": "Type117",
    "RHD*01W.118": "Type118",
    "RHD*01W.119": "Type119",
    "RHD*01W.12": "Type12",
    "RHD*01W.120": "Type120",
    "RHD*01W.121": "Type121",
    "RHD*01W.122": "Type122",
    "RHD*01W.123": "Type123",
    "RHD*01W.124": "Type124",
    "RHD*01W.125": "Type125",
    "RHD*01W.126": "Type126",
    "RHD*01W.127": "Type127",
    "RHD*01W.128": "Type128",
    "RHD*01W.129": "Type129",
    "RHD*01W.13": "Type13",
    "RHD*01W.130": "Type130",
    "RHD*01W.131": "Type131",
    "RHD*01W.132": "Type132",
    "RHD*01W.133": "Type133",
    "RHD*01W.134": "Type134",
    "RHD*01W.135": "Type135",
    "RHD*01W.136": "Type136",
    "RHD*01W.137": "Type137",
    "RHD*01W.138": "Type138",
    "RHD*01W.139": "Type139",
    "RHD*01W.14": "Type14",
    "RHD*01W.140": "Type140",
    "RHD*01W.141": "Type141",
    "RHD*01W.142": "Type142",
    "RHD*01W.143": "Type143",
    "RHD*01W.144": "Type144",
    "RHD*01W.145": "Type145",
    "RHD*01W.146": "Type146",
    "RHD*01W.147": "Type147",
    "RHD*01W.148": "Type148",
    "RHD*01W.149": "Type149",
    "RHD*01W.150": "Type150",
    "RHD*01W.151": "Type151",
    "RHD*01W.152": "Type152",
    "RHD*01W.153": "Type153",
    "RHD*01W.154": "Type154",
    "RHD*01W.155": "Type155",
    "RHD*01W.156": "Type156",
    "RHD*01W.157": "Type157",
    "RHD*01W.158": "Type158",
    "RHD*01W.159": "Type159",
    "RHD*01W.16": "Type16",
    "RHD*01W.160": "Type160",
    "RHD*01W.161": "Type161",
    "RHD*01W.162": "Type162",
    "RHD*01W.163": "Type163",
    "RHD*01W.164": "Type164",
    "RHD*01W.165": "Type165",
    "RHD*01W.17": "Type17",
    "RHD*01W.18": "Type18",
    "RHD*01W.19": "Type19",
    "RHD*01W.2": "Type2",
    "RHD*01W.2.1": "Type2.1",
    "RHD*01W.2.2": "Type2.2",
    "RHD*01W.20": "Type20",
    "RHD*01W.22": "Type22",
    "RHD*01W.23": "Type23",
    "RHD*01W.24": "Type24",
    "RHD*01W.25": "Type25",
    "RHD*01W.26": "Type26",
    "RHD*01W.27": "Type27",
    "RHD*01W.28": "Type28",
    "RHD*01W.29": "Type29",
    "RHD*01W.3": "Type3",
    "RHD*01W.3.1": "Type3.1",
    "RHD*01W.30": "Type30",
    "RHD*01W.31": "Type31",
    "RHD*01W.32": "Type32",
    "RHD*01W.33": "Type33",
    "RHD*01W.34": "Type34",
    "RHD*01W.35": "Type35",
    "RHD*01W.36": "Type36",
    "RHD*01W.37": "Type37",
    "RHD*01W.38": "Type38",
    "RHD*01W.39": "Type39",
    "RHD*01W.40": "Type40",
    "RHD*01W.41": "Type41",
    "RHD*01W.41.0.1": "Type41.0.1",
    "RHD*01W.42": "Type42",
    "RHD*01W.43": "Type43",
    "RHD*01W.44": "Type44",
    "RHD*01W.45": "Type45",
    "RHD*01W.45.1": "Type45.1",
    "RHD*01W.45.2": "Type45.2",
    "RHD*01W.46": "Type46",
    "RHD*01W.47": "Type47",
    "RHD*01W.48": "Type48",
    "RHD*01W.49": "Type49",
    "RHD*01W.5": "Type5",
    "RHD*01W.50": "Type50",
    "RHD*01W.51": "Type51",
    "RHD*01W.52": "Type52",
    "RHD*01W.53": "Type53",
    "RHD*01W.54": "Type54",
    "RHD*01W.55": "Type55",
    "RHD*01W.56": "Type56",
    "RHD*01W.57": "Type57",
    "RHD*01W.58": "Type58",
    "RHD*01W.59": "Type59",
    "RHD*01W.6": "Type6",
    "RHD*01W.60": "Type60",
    "RHD*01W.61": "Type61",
    "RHD*01W.62": "Type62",
    "RHD*01W.63": "Type63",
    "RHD*01W.64": "Type64",
    "RHD*01W.65": "Type65",
    "RHD*01W.66": "Type66",
    "RHD*01W.67": "Type67",
    "RHD*01W.68": "Type68",
    "RHD*01W.69": "Type69",
    "RHD*01W.7": "Type7",
    "RHD*01W.70": "Type70",
    "RHD*01W.71": "Type71",
    "RHD*01W.72": "Type72",
    "RHD*01W.73": "Type73",
    "RHD*01W.74": "Type74",
    "RHD*01W.75": "Type75",
    "RHD*01W.76": "Type76",
    "RHD*01W.77": "Type77",
    "RHD*01W.78": "Type78",
    "RHD*01W.79": "Type79",
    "RHD*01W.8": "Type8",
    "RHD*01W.80": "Type80",
    "RHD*01W.81": "Type81",
    "RHD*01W.82": "Type82",
    "RHD*01W.83": "Type83",
    "RHD*01W.84": "Type84",
    "RHD*01W.85": "Type85",
    "RHD*01W.86": "Type86",
    "RHD*01W.87": "Type87",
    "RHD*01W.88": "Type88",
    "RHD*01W.89": "Type89",
    "RHD*01W.9": "Type9",
    "RHD*01W.90": "Type90",
    "RHD*01W.91": "Type91",
    "RHD*01W.92": "Type92",
    "RHD*01W.93": "Type93",
    "RHD*01W.94": "Type94",
    "RHD*01W.95": "Type95",
    "RHD*01W.96": "Type96",
    "RHD*01W.97": "Type97",
    "RHD*01W.98": "Type98",
    "RHD*01W.99": "Type99",

    "RHD*02": "DII",
    "RHD*03.01": "DIIIa",
    "RHD*03.02": "DIIIb",
    "RHD*03.03": "DIIIc",
    "RHD*03.04": "DIIItype4",
    "RHD*03.04.02": "Nottested",
    "RHD*03.06": "DIIItype6",
    "RHD*03.07": "DIIItype7",
    "RHD*03.08": "DIIItype8",
    "RHD*03.09": "DIIItype9",
    "RHD*03N.01": "D-",
    "RHD*03N.02": "D-",
    "RHD*04.01": "DIVa",
    "RHD*04.01.02": "DIVa_like",
    "RHD*04.03": "DIVtype3",
    "RHD*04.04": "DIVtype4",
    "RHD*04.05": "DIVtype5",
    "RHD*04.06": "DIVb",
    "RHD*05.01": "DVtype1",
    "RHD*05.02": "DVtype2",
    "RHD*05.03": "DVtype3",
    "RHD*05.04": "DVtype4",
    "RHD*05.05": "DVtype5",
    "RHD*05.06": "DVtype6",
    "RHD*05.07": "DVtype7",
    "RHD*05.08": "DVtype8",
    "RHD*05.09": "DVtype9",
    "RHD*05.10": "DVtype10",
    "RHD*06.01": "DVItype1",
    "RHD*06.02": "DVItype2",
    "RHD*06.03.01": "DVItype3",
    "RHD*06.03.02": "DVItype3.2",
    "RHD*06.04": "DVItype4",
    "RHD*07.01": "DVII",
    "RHD*07.02": "DVIItype2",
    "RHD*08.01": "DFV",
    "RHD*08N.01": "D-",
    "RHD*09.01": "DAR(T203A)",
    "RHD*09.01.00": "DAR1(weakD4.2)",
    "RHD*09.01.01": "DAR1.1(weakD4.2.1)",
    "RHD*09.01.02": "DAR1.2(weakD4.2.2)",
    "RHD*09.01.03": "DAR1.3(weakD4.2.3)",
    "RHD*09.02": "DAR2(DARE)",
    "RHD*09.02.01": "DAR2.1",
    "RHD*09.03": "DAR3(weakPartial_D4.0.1)",
    "RHD*09.03.01": "DAR3.1(weakPartial_D4.0)",
    "RHD*09.04": "DAR4(weakD4.1)",
    "RHD*09.05": "DAR5(weakD4.3orDel)",
    "RHD*09.06": "DAR6",
    "RHD*10.00": "DAU0",
    "RHD*10.00.01": "DAU0.01",
    "RHD*10.00.02": "DAU0.02",
    "RHD*10.01": "DAU1",
    "RHD*10.02": "DAU2",
    "RHD*10.03": "DAU3",
    "RHD*10.04": "DAU4",
    "RHD*10.05": "DAU5",
    "RHD*10.05.01": "DAU5.1",
    "RHD*10.06": "DAU6",
    "RHD*10.07": "DAU7",
    "RHD*10.08": "DAU8",
    "RHD*10.09": "DAU9",
    "RHD*10.10": "DAU10",
    "RHD*10.11": "DAU11",
    "RHD*10.12": "DAU12",
    "RHD*10.13": "DAU13",
    "RHD*10.14": "DAU14",
    "RHD*10.15": "RHD(M1V,T379M)",
    "RHD*11": "weakPartial_11orDel",
    "RHD*12.01": "DOL1",
    "RHD*12.02": "DOL2",
    "RHD*12.03": "DOL3",
    "RHD*12.04": "DOL4",
    "RHD*13.01": "DBS1",
    "RHD*13.02": "DBS2",
    "RHD*14.01": "DBT1",
    "RHD*14.02": "DBT2",
    "RHD*15": "WeakPartial_Type15",
    "RHD*16.01": "DCS1",
    "RHD*16.02": "DCS2",
    "RHD*16.03": "DCS3",
    "RHD*17.01": "DFR1",
    "RHD*17.02": "DFR2",
    "RHD*17.03": "DFR3",
    "RHD*17.04": "DFR4",
    "RHD*17.05": "DFR5",
    "RHD*18": "DFW",
    "RHD*19": "DHMi",
    "RHD*20": "DHO",
    "RHD*21": "D+weak_partial",
    "RHD*22": "DHR",
    "RHD*23": "DMH",
    "RHD*24": "DNAK",
    "RHD*25": "DNB",
    "RHD*26": "DNU",
    "RHD*27": "DDE",
    "RHD*28": "DFL",
    "RHD*29": "DYU",
    "RHD*30": "DQC",
    "RHD*31": "DQC",
    "RHD*32": "DVL2",
    "RHD*33": "DWI(DWLLE)",
    "RHD*34": "DIM(DIleM)",
    "RHD*35": "DMA",
    "RHD*36": "DLO",
    "RHD*37": "DUC2",
    "RHD*38": "DNT",
    "RHD*39": "RHD(S103P),G-",
    "RHD*40": "D-SPM",
    "RHD*41": "DBU",
    "RHD*42": "DCC",
    "RHD*43": "DDN",
    "RHD*44": "DHQ",
    "RHD*45": "DKK",
    "RHD*46": "DLX",
    "RHD*47": "DMI",
    "RHD*47.01": "DMI_1.1",
    "RHD*48": "DNS",
    "RHD*49": "DWN",
    "RHD*50": "RHD(A354T)",
    "RHD*51": "RHD(del44L)",
    "RHD*52": "RHD(F223S)",
    "RHD*53": "RHD(IVS2_2delA)",
    "RHD*54": "RHD(IVS4_2A>C)",
    "RHD*55": "RHD(L81P)",
    "RHD*56": "DBA",
    "RHD*57": "weakPartial_type57",
    "RHD*58": "D-CE(7)-D",
    "RHD*59": "RHD(F175L)",
    "RHD*60": "D+weak_partial",
    "RHD*61": "D+weak_partial",
    "RHD*62": "DNT(V270G)",
    "RHD*63": "D+partial",
    "RHD*64": "D+partial",
    "RHD*65": "D+partial",
    "RHD*66": "D+partial",
}
