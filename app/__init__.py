from flask import Flask, render_template
from flask_restplus import Api
from config import Config
from celery import Celery
import os

def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL']
    )
    celery.conf.update(app.config)
    
    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)
    celery.Task = ContextTask
    return celery

api = Api()

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html", index=True)
    
# app.config.from_object(Config)
app.config.update(
    SECRET_KEY = os.environ.get('SECRET_KEY') or b'6\xe9\xda\xead\x81\xf7\x8d\xbbH\x87\xe8m\xdd3%',

    # Location to store all the data
    IMMUNOLYSER_DATA = os.environ.get('IMMUNOLYSER_DATA'),

    # Task id of demo job
    DEMO_TASK_ID = os.environ.get('DEMO_TASK_ID'),
    CELERY_BROKER_URL='redis://localhost:6379/0',
    CELERY_RESULT_BACKEND='db+sqlite:///results.sqlite',  # SQLite
    CELERY_DEFAULT_QUEUE='celery',  # Ensure all tasks are routed to 'celery' queue
    DEBUG = True,
    PIN = '123'
)

celery = make_celery(app)

@app.route('/process/<name>')
def process(name):
    task = reverse_name.delay(name)
    return f'Request to process {name} has been received. Task ID is {task.id}'

@celery.task(name='myapp.reverse_name', bind=True)
def reverse_name(self, name):
    try:
        # Simulate an error condition
        # raise ValueError("Input must be a string")
        
        result = name[::-1]
        return result
    except Exception as e:
        # Manually raise an exception to simulate failure
        raise e

api.init_app(app)

from app import routes

app.run(debug=True)