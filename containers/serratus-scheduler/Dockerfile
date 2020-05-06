FROM ubuntu

RUN apt-get update -yq
RUN apt-get install -yqq python3-pip

RUN pip3 install flask sqlalchemy gunicorn prometheus_client boto3

# Run Python in UTF-8 mode
ENV FLASK_APP=flask_app
ENV FLASK_ENV=production
ENV LANG=C.UTF-8
ENV PYTHONUNBUFFERED=1

WORKDIR /opt/
ADD serratus-scheduler/flask_app/ ./flask_app

# If something fails, `flask init-db` prints a usage message, not the error.
# So run "flask --help" in that case which does print the error
RUN flask init-db || (flask --help; false)

# Do NOT enable multithreading here, as the underlying app is
# thread-unsafe.
CMD gunicorn 'flask_app:create_app()' --bind '0.0.0.0:8000'
#CMD flask run -p 8000

