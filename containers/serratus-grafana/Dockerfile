FROM grafana/grafana
COPY serratus-grafana/grafana.ini /etc/grafana/grafana.ini
COPY serratus-grafana/datasources.yaml /etc/grafana/provisioning/datasources/
COPY serratus-grafana/dashboards.yaml /etc/grafana/provisioning/dashboards/
COPY serratus-grafana/dashboards /var/lib/grafana/dashboards/
