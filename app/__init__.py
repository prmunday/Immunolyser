from flask import Flask
from flask_restplus import Api
from config import Config

api = Api()

app = Flask(__name__)
app.config.from_object(Config)

api.init_app(app)

from app import routes

app.run()