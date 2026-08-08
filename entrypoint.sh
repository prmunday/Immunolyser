#!/bin/bash
set -e

# Activate virtual environment
source lenv/bin/activate

# app/static/images is where job output (prediction CSVs, seqlogos, gibbs-
# cluster/MHC-TP results, allele_compatibility_matrix.csv, ...) gets written
# during processing, and where the web UI reads it back from. On bare-metal
# this is a symlink to /pvol (the persistent data volume) so it's the same
# filesystem everywhere. flask_app and celery_worker are SEPARATE containers
# though — without this symlink, each has its own private, ephemeral
# app/static/images/, so the worker's output is invisible to the web
# container serving report pages. Both containers already mount /pvol
# (see docker-compose.yml), so replicate the same symlink here.
if [ ! -L app/static/images ]; then
    rm -rf app/static/images
    ln -s /pvol app/static/images
    echo "Linked app/static/images -> /pvol"
fi

# Seed the bind-mounted ref_data volume with the committed defaults (mouse
# Class I motifs, human/mouse .db, warmed numba_cache) the first time it's
# empty — the bind mount in docker-compose.yml hides whatever git cloned
# there, so an empty host dir would otherwise start with nothing at all.
REF_DATA_DIR="app/tools/HLA-PepClust/data/ref_data"
if [ -d "/app/.hlapepclust-ref-data-seed" ] && [ -z "$(ls -A "$REF_DATA_DIR" 2>/dev/null)" ]; then
    echo "Seeding $REF_DATA_DIR from image defaults (mouse Class I + numba cache) ..."
    cp -r /app/.hlapepclust-ref-data-seed/. "$REF_DATA_DIR"/
fi

# DB init (Python-based — keeps schema in sync with job_registry.py / email_registry.py)
DB_FILE="${IMMUNOLYSER_DATA}/results.sqlite"
if [ "$1" = "flask" ]; then
    if [ ! -f "$DB_FILE" ]; then
        echo "Initialising SQLite database at $DB_FILE ..."
        python3 -c "
from app.job_registry import init_job_registry
from app.email_registry import init_email_registry
init_job_registry()
init_email_registry()
print('Database initialised.')
"
    else
        echo "Database already exists at $DB_FILE"
    fi

    # Warn about missing licensed tools (non-fatal)
    for tool in \
        "app/tools/netMHCpan-4.2/netMHCpan:netMHCpan 4.2" \
        "app/tools/netMHCIIpan-4.3/netMHCIIpan:netMHCIIpan 4.3"; do
        path="${tool%%:*}"
        name="${tool##*:}"
        if [ ! -f "$path" ]; then
            echo "WARNING: $name not found at $path — predictions using this tool will fail."
            echo "         See app/tools/README.md for download instructions."
        fi
    done
fi

# Run the service
if [ "$1" = "flask" ]; then
    exec gunicorn --workers 2 --bind 0.0.0.0:5000 --timeout 120 firstdemo:app
elif [ "$1" = "celery" ]; then
    exec celery -A app.celery worker --loglevel=info --concurrency=1
elif [ "$1" = "beat" ]; then
    # --schedule points the persistent schedule file (last-run timestamps, used
    # to avoid re-firing a task on restart) at /pvol instead of the container's
    # ephemeral filesystem, so it survives container recreation.
    exec celery -A app.celery beat --loglevel=info --schedule=/pvol/celerybeat-schedule
else
    exec "$@"
fi
