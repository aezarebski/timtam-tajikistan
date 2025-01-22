#!/bin/bash

# Define variables
BEAST_BASE_URL="https://github.com/CompEvol/beast2/releases/download/v2.7.3"
BEAST_TAR_FILE="BEAST.v2.7.3.Linux.x86.tgz"
TRACER_BASE_URL="https://github.com/beast-dev/tracer/releases/download/v1.7.2"
TRACER_TAR_FILE="Tracer_v1.7.2.tgz"
DEST_DIR="lib"

# Create the lib directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Download the BEAST tarball
curl -L "${BEAST_BASE_URL}/${BEAST_TAR_FILE}" -o "${DEST_DIR}/${BEAST_TAR_FILE}"

# Extract the BEAST tarball into the lib directory
tar -xvzf "${DEST_DIR}/${BEAST_TAR_FILE}" -C "$DEST_DIR"

# Set permissions for required BEAST files
chmod 750 "${DEST_DIR}/beast/bin/beast"
chmod 750 "${DEST_DIR}/beast/bin/beauti"
chmod 750 "${DEST_DIR}/beast/jre/bin/java"

echo "BEAST 2.7.3 setup completed."

# Download the Tracer tarball
curl -L "${TRACER_BASE_URL}/${TRACER_TAR_FILE}" -o "${DEST_DIR}/${TRACER_TAR_FILE}"

# Extract the Tracer tarball into the lib directory
mkdir "${DEST_DIR}/tracer"
tar -xvzf "${DEST_DIR}/${TRACER_TAR_FILE}" -C "${DEST_DIR}/tracer"

echo "Tracer 1.7.2 setup completed."
