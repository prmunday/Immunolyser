from app import app
from flask import render_template, request, redirect, url_for
from app.forms import InitialiserForm
import pandas as pd
import os
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile
from app.sample import Sample

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    return render_template("index.html", index=True )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    form = InitialiserForm()
    if form.validate_on_submit():
        
        sampleName = form.col_name.name
        
        # Have to take this input from user
        maxLen = 30
        minLen = 1
        # filename = secure_filename(f.filename)
        # f.save(os.path.join('data', filename))
        files_filenames = []
        for smaple in form.file_input.data:
            file_filename = secure_filename(smaple.filename)
            smaple.save(os.path.join('data', file_filename))
            files_filenames.append(file_filename)
        print(files_filenames)

        sample = Sample(sampleName)

        for replicate in files_filenames:
            sample.addReplicate(replicate, pd.read_csv(os.path.join('data', replicate)))

        # Have to later add the user input for length
        sample.data = filterPeaksFile(sample.data, minLen=minLen, maxLen=maxLen)
        bar = plot_lenght_distribution({sample.name:sample})

        return render_template('analytics.html', plot=bar)
    return render_template("initialiser.html", form=form)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html")
