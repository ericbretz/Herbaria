import os

# script path and where the data is cached
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(os.path.dirname(SCRIPT_DIR), '02_data')
CACHE_FILE = os.path.join(SCRIPT_DIR, 'plot_data_cache.pkl')

SAMPLES = {
    'WA02'  : 'WA02_paired', 'WA03'  : 'WA03_paired', 'WA07'  : 'WA07_paired',
    'WA08'  : 'WA08_paired', 'WA09'  : 'WA09_paired', 'WA11'  : 'WA11_paired',
    'WA12'  : 'WA12_paired', 'WA13'  : 'WA13_paired', 'WA14'  : 'WA14_paired',
    'WA18'  : 'WA18_paired', 'WA22'  : 'WA22_paired', 'WA25'  : 'WA25_paired',
    'WA26'  : 'WA26_paired', 'DAL192': 'DAL192',      'DAL193': 'DAL193',
    'DAL195': 'DAL195',      'DAL212': 'DAL212',      'DAL218': 'DAL218',
    'DAL224': 'DAL224',      'DAL227': 'DAL227'
}

SAMPLE_NAMES = {
    'WA02'  : 'Darlingtonia1-Silica',  'WA03'  : 'Ipomopsis-Silica',
    'WA07'  : 'Leucothoe-Silica',      'WA08'  : 'Darlingtonia2-Silica',
    'WA09'  : 'Primula-Silica',        'WA11'  : 'Pyrola-Silica',
    'WA12'  : 'Navarretia-Silica',     'WA13'  : 'Collomia-Silica',
    'WA14'  : 'Vaccinium-Silica',      'WA18'  : 'Cyclamen-Silica',
    'WA22'  : 'Antistrophe-Herbarium', 'WA25'  : 'Jacquinia-Herbarium',
    'WA26'  : 'Clavija-Herbarium',     'DAL192': 'Collomia-Fresh',
    'DAL193': 'Navarretia-Fresh',      'DAL195': 'Pyrola-Fresh',
    'DAL212': 'Vaccinium-Fresh',       'DAL218': 'Ipomopsis-Fresh',
    'DAL224': 'Leucothoe-Fresh',       'DAL227': 'Darlingtonia1-Fresh'
}

SAMPLE_TYPES           = {'WA22': 'herbaria', 'WA13': 'silica', 'DAL192': 'fresh'}
REPRESENTATIVE_SAMPLES = ['WA22', 'DAL192', 'WA13']

COLORS = {
    'herbaria': '#BA6F6D',
    'fresh'   : '#C5CF9D',
    'silica'  : '#80B0B8',
    'busco': {
        'S': '#56B4E9',
        'D': '#3492C7',
        'F': '#F0E442',
        'M': '#F04442',
    },
    'transrate': {
        'sCnuc_Harmonic': '#264653',
        'sCcov_Harmonic': '#C74F41',
        'sCord_Harmonic': '#E9C46A',
        'sCseg_Harmonic': '#4D9267',
    },
}

BUSCO_CATEGORIES  = ['S', 'D', 'F', 'M']
BUSCO_LABELS      = ['Complete\nSingle-copy', 'Complete\nDuplicated', 'Fragmented', 'Missing']

TRANSRATE_COLS    = ['sCnuc_Harmonic', 'sCcov_Harmonic', 'sCord_Harmonic', 'sCseg_Harmonic']
TRANSRATE_LABELS  = ['sCnuc', 'sCcov', 'sCord', 'sCseg']

TYPE_ORDER        = ('fresh', 'silica', 'herbaria')
TYPE_DISPLAY      = {'fresh': 'Fresh', 'silica': 'Silica', 'herbaria': 'Herbarium'}

REP_SAMPLE_LABELS = {
    'DAL192': 'Collomia\nFresh',
    'WA13':   'Collomia\nSilica',
    'WA22':   'Antistrophe\nHerbarium',
}

CONCAT_PLOT_ORDER = ['deamination', 'inserts', 'transrate', 'busco']
CONCAT_SAMPLES    = ['DAL192', 'WA13', 'WA22']

ROW_TITLES = {
    'deamination': 'Deamination',
    'inserts':     'Insert Length',
    'transrate':   'TransRate2',
    'busco':       'BUSCO',
}

FIGSIZE = {
    'busco':       (8,  6), # 4 plots = 8
    'deamination': (10, 6), # 3 plots = 10
    'inserts':     (10, 6),
    'transrate':   (8,  6),
}

FIGSIZE_REP              = (38.48, 13.2) # 52.76 / 4
FIGSIZE_TRANSRATE_SCORES = (16,    12)
FIGSIZE_CONCAT           = (38.48, 52.76) # tried to make this proportional to an 8x11 page

TITLE_SUFFIX = {
    'busco':       ' BUSCO Scores',
    'deamination': ' Deamination Pattern',
    'inserts':     ' Insert Length Distribution',
    'transrate':   ' Transrate2 Scores',
}


# Font info 
SCALE = 2.0 # the font proportions i needed to fine tune, this is for the sizes

# initially representative and individual plots had different sizes,
# but then realized that the plots would look different
# keeping here incase need to revert back for whatever reason
FONTS_REP = {
    'title'   : 20*SCALE, 'axis_label'    : 22*SCALE, 'tick_label': 22*SCALE,
    'legend'  : 22*SCALE, 'bar_annotation': 18*SCALE, 'error'     : 18*SCALE,
    'labelpad': 10,       'title_pad'     : 20,
}

FONTS_IND = {
    'title'   : 16, 'axis_label'    : 12, 'tick_label': 11,
    'legend'  : 10, 'bar_annotation': 12, 'error'     : 12,
    'labelpad': 6,  'title_pad'     : 10,
}

FONTS_CONCAT = {
    'title'   : 20*SCALE, 'axis_label'    : 18*SCALE, 'tick_label': 16*SCALE,
    'legend'  : 16*SCALE, 'bar_annotation': 18*SCALE, 'error'     : 20*SCALE,
    'labelpad': 25,       'title_pad'     : 20,
}

FONTS_SLIDE = {
    'title'   : 7,   'axis_label'    : 5.5, 'tick_label': 5,
    'legend'  : 4.5, 'bar_annotation': 5,   'error'     : 5,
    'labelpad': 2,   'title_pad'     : 3,
}

FONT_ROW_LABEL = 48 * SCALE
FONT_ROW_TITLE = 24 * SCALE

# as close as i could get to powerpoint dimensions
SLIDE_FIGSIZE = (13.33, 7.5)
SLIDE_NCOLS   = 5
SLIDE_NROWS   = 4
