FROM alpine:3.11

RUN apk add \
        python3 postgresql-client py3-flask py3-sqlalchemy py3-gunicorn \
        py3-prometheus-client py3-psycopg2 py3-gevent && \
    pip3 install boto3 psycogreen # boto3 not in 3.11

# Run Python in UTF-8 mode
ENV FLASK_APP=flask_app
ENV FLASK_ENV=production
ENV LANG=C.UTF-8
ENV PYTHONUNBUFFERED=1

RUN mkdir /opt/metrics
ENV prometheus_multiproc_dir=/opt/metrics

WORKDIR /opt/
COPY serratus-scheduler/flask_app/ ./flask_app
COPY serratus-scheduler/gunicorn.conf.py /etc/gunicorn.conf.py

ENTRYPOINT ["gunicorn", "flask_app:create_app()", "--bind", "0.0.0.0:8000", "--config", "/etc/gunicorn.conf.py"]
