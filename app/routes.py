from app import app
from flask import render_template, request, redirect, url_for
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile
from app.sample import Sample
import time

# Experiment ID
EXPERIMENT_ID = 0

@app.route("/")
@app.route("/index")
@app.route("/home")
def index():
    return render_template("index.html", index=True )

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():
    form = ParentForm()

    if form.add_sample.data:
        form.children.append_entry()
        return render_template("initialiser.html", form=form, initialiser=True)

    if form.validate_on_submit():

        experiment_id = getExperminetId()
        

        samples = []
        for sample in form.children.data:
            samples.append(sample)
    #     sample_one_name = form.sample_one_name.data
    #     sample_two_name = form.sample_two_name.data
    
        # Have to take this input from user
        maxLen = 30
        minLen = 1

        # Creating directory for the experiment
        # Create directory
        dirName = os.path.join('data', str(experiment_id)+'-ID-'+time.strftime("%Y-%m-%d"))
        try:
            # Create target Directory
            os.makedirs(dirName)
            print("Directory " , dirName ,  " Created ") 
        except FileExistsError:
            print("Directory " , dirName ,  " already exists")

        # Creating sub directories to store sample data
        for sample in samples:
            try:
                os.mkdir(os.path.join(dirName, sample['sample_name']))
                print("Directory Created") 
            except FileExistsError:
                print("Directory already exists")     

        data = {}

        # Saving the data and loading into the dictionary
        for sample in samples:
            data[sample['sample_name']] = list()
            for replicate in sample['sample']:
                file_filename = secure_filename(replicate.filename)
                replicate.save(os.path.join(dirName, sample['sample_name'], file_filename))
                # sample_one.append(file_filename)
                data[sample['sample_name']].append(file_filename)

    #     # Storing in classes
    #     # sample1 = Sample(sample_one_name)
    #     # sample2 = Sample(sample_two_name)

    #     # Or in variables
    #     sample1 = {}
    #     sample2 = {}

        sample_data = {}

        for sample_name, file_names in data.items():
            sample_data[sample_name] = dict()
            for replicate in file_names:
                sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

        # Have to later add the user input for length
        for sample_name, sample in sample_data.items():
            sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


        bar_percent = plot_lenght_distribution(sample_data, hist='percent')
        bar_density = plot_lenght_distribution(sample_data, hist='density')

        return render_template('analytics.html', peptide_percent=bar_percent, peptide_density=bar_density, analytics=True)
    return render_template("initialiser.html", form=form, initialiser=True)

@app.route("/analytics")
def analytics():
    return render_template("analytics.html", analytics=True)

# Method to manage experiment ID
def getExperminetId():
    global EXPERIMENT_ID
    EXPERIMENT_ID = EXPERIMENT_ID+1

    return EXPERIMENT_ID