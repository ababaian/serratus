FROM ubuntu

RUN apt-get update -yq
RUN apt-get install -yqq python3-pip

RUN pip3 install flask sqlalchemy gunicorn prometheus_client boto3 psycopg2-binary

# Run Python in UTF-8 mode
ENV FLASK_APP=flask_app
ENV FLASK_ENV=production
ENV LANG=C.UTF-8
ENV PYTHONUNBUFFERED=1

WORKDIR /opt/
ADD serratus-scheduler/flask_app/ ./flask_app

ENTRYPOINT ["gunicorn", "flask_app:create_app()", "--bind", "0.0.0.0:8000"]
