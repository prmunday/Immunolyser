from app import app
from flask import render_template, request, redirect, url_for
from app.forms import InitialiserForm
import plotly
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import json
import os
from werkzeug.utils import secure_filename

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    return render_template("index.html", index=True )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    form = InitialiserForm()
    if form.validate_on_submit():

        col_name   = form.col_name.data
        f = form.file_input.data
        filename = secure_filename(f.filename)
        f.save(os.path.join('data', filename))

        bar = create_plot(filename, col_name)

        return render_template('analytics.html', plot=bar)
    return render_template("initialiser.html", form=form)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html")

def create_plot(filename, col_name):

    print(filename, col_name)
    data = pd.read_csv('data/VMM1_1st_nil_DT9_peptide.csv')


    data = [
        go.Histogram(
            x=data['Length']
        )
    ]

    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON