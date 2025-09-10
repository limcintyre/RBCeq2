from enum import Enum, auto

# Define version
VERSION = "2.3.4"
DB_VERSION = "2.3.3"


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
            "8": ("9", "52"),
            "9": ("8", "52"),
            "52": ("8", "9"),
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
            "Cw": ("Cx", "BARC"),
            "Cx": ("Cw", "BARC"),
            "BARC": ("Cw", "Cx"),
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

LANE = {
    "chr1": {
        "1000": "no_ALT",  # unit test
        "159205564": "G_A",
        "207331122": "G_C",
        "207609424": "G_A",
        "207609571": "A_T",
        "25317062": "T_C",
        "25390874": "C_G",
        "25408711": "G_A",
        "25420739": "G_C",
        "42830851": "G_A",
    },
    "chr11": {"35176644": "G_C"},
    "chr12": {"14840505": "C_T"},
    "chr16": {"88716069": "C_T"},
    "chr17": {
        "42453065": "A_C",
        "42453072": "G_T",
        "42453084": "C_T",
        "42453291": "C_G",
        "42453713": "C_A",
        "42455875": "G_A",
        "42457790": "C_T",
        "42462694": "T_G",
        "44251253": "G_A",
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
        "45739554": "G_A",
    },
    "chr19": {
        "44812188": "G_A",
        "5844526": "no_ALT",
        "10287311": "A_G",
        "5844638": "no_ALT",
        "10631494": "no_ALT",
    },
    "chr22": {
        "19711485": "G_A",
        "42717787": "C_A",
    },
    "chr3": {
        "161086379": "C_T",
        "128780950": "C_T",
    },
    "chr4": {
        "144120555": "T_C",
        "144120567": "A_G",
        "143997559": "C_G",
        "143999443": "G_A",
        "144120554": "C_A",
    },
    "chr5": {
        "52366090": "G_T",
        "52358757": "G_A",
        "52369001": "C_T",
        "52382870": "C_T",
    },
    "chr6": {
        "44232918": "G_A",
        "10586805": "no_ALT",
    },
    "chr7": {"142957921": "G_A", "100893176": "G_T"},
    "chr8": {
        "74493432": "C_A",
    },
    "chr9": {"133257521": "T_TC"},
}
