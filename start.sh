echo Starting Gunicorn.
exec gunicorn holmes_core.wsgi:application \
    --bind 0.0.0.0:8000 \
    --workers 3
