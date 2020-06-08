FROM alpine:3.11

RUN apk add \
        python3 postgresql-client py3-flask py3-sqlalchemy py3-gunicorn \
        py3-prometheus-client py3-psycopg2 && \
    pip3 install boto3 # boto3 not in 3.11

# Run Python in UTF-8 mode
ENV FLASK_APP=flask_app
ENV FLASK_ENV=production
ENV LANG=C.UTF-8
ENV PYTHONUNBUFFERED=1

WORKDIR /opt/
COPY serratus-scheduler/flask_app/ ./flask_app

ENTRYPOINT ["gunicorn", "flask_app:create_app()", "--bind", "0.0.0.0:8000"]
