#!/bin/bash
set -e

# Activate virtual environment
source lenv/bin/activate

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
else
    exec "$@"
fi
