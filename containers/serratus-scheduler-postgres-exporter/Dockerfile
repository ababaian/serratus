FROM wrouesnel/postgres_exporter:latest
COPY serratus-scheduler-postgres-exporter/queries.yaml /queries.yaml
CMD ["--extend.query-path", "/queries.yaml", "--disable-settings-metrics"]
