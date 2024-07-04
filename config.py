import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or b'6\xe9\xda\xead\x81\xf7\x8d\xbbH\x87\xe8m\xdd3%'

    # Location to store all the data
    IMMUNOLYSER_DATA = os.environ.get('IMMUNOLYSER_DATA')

    # Task id of demo job
    DEMO_TASK_ID = os.environ.get('DEMO_TASK_ID')

    CELERY_BROKER_URL='redis://localhost:6379/0'
    CELERY_RESULT_BACKEND='db+sqlite:///results.sqlite'  # SQLite
    CELERY_DEFAULT_QUEUE='celery'  # Ensure all tasks are routed to 'celery' queue
    DEBUG = True
    PIN = '123'