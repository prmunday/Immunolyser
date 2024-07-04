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
    
app.config.from_object(Config)

celery = make_celery(app)

api.init_app(app)

from app import routes

# app.run(debug=True)