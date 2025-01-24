#!/bin/bash

# Activate the virtual environment
source lenv/bin/activate

# Set Flask environment variables
export FLASK_APP=firstdemo.py
export FLASK_ENV=development
export IMMUNOLYSER_DATA=${IMMUNOLYSER_DATA}

# Run the appropriate process
if [ "$1" = "flask" ]; then
    flask run --host=0.0.0.0 --port=5000
elif [ "$1" = "celery" ]; then
    celery -A app.celery worker --loglevel=info --concurrency=1
else
    exec "$@"
fi
