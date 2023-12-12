from flask import Flask, render_template
from flask_restplus import Api
from config import Config

api = Api()

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html", index=True)
    
app.config.from_object(Config)

api.init_app(app)

from app import routes

app.run(host='0.0.0.0', port=5000)