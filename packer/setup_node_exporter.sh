#!/bin/bash
set -eux

# Grab the binary from prometheus github
curl -sL https://github.com/prometheus/node_exporter/releases/download/v1.0.0-rc.0/node_exporter-1.0.0-rc.0.linux-amd64.tar.gz -o node_exporter.tgz
tar -xvzf node_exporter.tgz

# Install it
sudo cp node_exporter-1.0.0-rc.0.linux-amd64/node_exporter /usr/bin/

# Create a service, which runs automatically in the background.
sudo useradd -s /bin/false prometheus
sudo tee /etc/systemd/system/node_exporter.service << END
[Unit]
Description=Prometheus Node Exporter
Wants=network-online.target
After=network-online.target

[Service]
User=prometheus
ExecStart=/usr/bin/node_exporter

[Install]
WantedBy=default.target
END
sudo systemctl enable node_exporter.service
