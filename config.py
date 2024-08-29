import os

class Config(object):
    # Security
    SECRET_KEY = os.environ.get('SECRET_KEY') or b'6\xe9\xda\xead\x81\xf7\x8d\xbbH\x87\xe8m\xdd3%'
    
    # Data Storage
    IMMUNOLYSER_DATA = os.environ.get('IMMUNOLYSER_DATA')
    
    # Celery Configuration
    CELERY_BROKER_URL = 'redis://localhost:6379/0'
    CELERY_RESULT_BACKEND = 'db+sqlite:///results.sqlite'  # SQLite
    CELERY_DEFAULT_QUEUE = 'celery'  # Ensure all tasks are routed to 'celery' queue
    
    # Demo Task
    DEMO_TASK_ID = os.environ.get('DEMO_TASK_ID')

    # Debugging
    DEBUG = True

    # Application Settings
    PIN = '123'
    
    # Job Input Limits
    SAMPLE_NAME_MAX_LENGTH = 30
    MAX_SAMPLES = 5
    MAX_TOTAL_PEPTIDES = 100000

    # Dropdown Option Sets
    OPTION_SETS = {
        'Species': [
            {'value': 'human', 'label': 'Human'},
            {'value': 'mouse', 'label': 'Mouse (H-2)'}
        ],
        'MHCClass': [
            {'value': 'mhc1', 'label': 'MHC I'},
            {'value': 'mhc2', 'label': 'MHC II'}
        ],
        'BindingTool': [
            {'value': 'NetMHCpan_4_1', 'label': 'NetMHCpan 4.1'},
            {'value': 'MixMHCpred_2_1', 'label': 'MixMHCpred 2.1'},
            {'value': 'NetMHCIIpan_4_0', 'label': 'NetMHCIIpan 4.0'},
            {'value': 'MixMHC2pred_2_0', 'label': 'MixMHC2pred 2.0'}
        ]
    }
