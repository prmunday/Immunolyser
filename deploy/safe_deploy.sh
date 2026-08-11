#!/bin/bash
# Safe deploy wrapper for the Docker stack.
#
# `docker compose build && docker compose up -d` (or any recreate of
# celery_worker/flask_app) kills whatever job is currently running on the
# worker — SIGTERM/SIGKILL takes out the Celery process and any subprocess
# it spawned (NetMHCpan, GibbsCluster, ...) mid-write, with no error logged.
# The job is just left stuck at its last status forever. This happened for
# real on 2026-08-10 (job 842149c6, mid-GibbsCluster, killed by a routine
# deploy). See docker-compose.yml's celery_worker.stop_grace_period comment
# for the same story from the container-shutdown side.
#
# This script checks for active/reserved Celery jobs first and requires
# explicit confirmation (or --force) before proceeding, instead of silently
# taking them out. Run this INSTEAD of raw `docker compose build/up`
# commands for any change that touches celery_worker or flask_app.
#
# Usage:
#   ./deploy/safe_deploy.sh                # checks, prompts if jobs active
#   ./deploy/safe_deploy.sh --force        # skip the check entirely
#   ./deploy/safe_deploy.sh --check-only   # just report active jobs, exit

set -euo pipefail
cd "$(dirname "$0")/.."

FORCE=false
CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=true ;;
        --check-only) CHECK_ONLY=true ;;
    esac
done

echo "=== Checking for active Celery jobs on celery_worker ==="
ACTIVE_OUTPUT=$(docker exec celery_worker bash -c 'source lenv/bin/activate && celery -A app.celery inspect active 2>&1' || true)
echo "$ACTIVE_OUTPUT"

# "- empty -" appears per-node when nothing is active; anything scheduling an
# app.submit_job task means real prediction work is in flight (the frequent
# dispatch_pending_jobs housekeeping task doesn't count as "a job").
if echo "$ACTIVE_OUTPUT" | grep -q "app.submit_job"; then
    echo
    echo "!!! A prediction job is currently running on celery_worker !!!"
    echo "Recreating celery_worker or flask_app now will kill it with no error logged —"
    echo "it'll be stuck at its last status in job_registry forever."
    echo

    if $CHECK_ONLY; then
        exit 1
    fi

    if ! $FORCE; then
        read -r -p "Proceed anyway and kill the running job? [y/N] " reply
        if [[ ! "$reply" =~ ^[Yy]$ ]]; then
            echo "Aborted. Wait for the job to finish, or re-run with --force to override."
            exit 1
        fi
    fi
else
    echo "No active prediction jobs. Safe to proceed."
    if $CHECK_ONLY; then
        exit 0
    fi
fi

echo
echo "=== Deploying ==="
git pull
docker compose build
docker compose up -d
echo
echo "=== Done. Current status: ==="
docker compose ps
