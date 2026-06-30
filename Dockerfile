FROM python:3.11-slim

WORKDIR /app

# Install dependencies first — separate layer so rebuilds are fast when only code changes
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy only the files the server needs at runtime
COPY server.py snapshot.py .env.example ./

# IBM_QUANTUM_TOKEN must be passed at runtime via -e or --env-file
# Never bake secrets into the image
EXPOSE 3020

CMD ["python", "server.py", "--transport", "http", "--host", "0.0.0.0", "--port", "3020"]
