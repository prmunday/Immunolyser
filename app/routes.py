from app import app
from flask import render_template, request, redirect, url_for
from app.forms import InitialiserForm
import pandas as pd
import os
from werkzeug.utils import secure_filename
from app.utils import create_plot, filterPeaksFile

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    return render_template("index.html", index=True )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    form = InitialiserForm()
    if form.validate_on_submit():

        f = form.file_input.data
        filename = secure_filename(f.filename)
        f.save(os.path.join('data', filename))

        data = pd.read_csv('data/{}'.format(filename))
        data = filterPeaksFile(data, maxLen=30)
        bar = create_plot(data, normalized = False)

        return render_template('analytics.html', plot=bar)
    return render_template("initialiser.html", form=form)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html")
