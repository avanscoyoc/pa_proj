FROM ghcr.io/prefix-dev/pixi:0.41.4

# Set working directory
WORKDIR /workspace

# Copy all project files including pixi.toml, pixi.lock, and source code
COPY . .

RUN apt-get update && apt-get install -y git

# Install all dependencies into the default Pixi environment (.pixi/env)
RUN pixi install --locked

# Generate an entrypoint that activates the Pixi environment
RUN pixi shell-hook -s bash > /shell-hook \
 && echo "#!/bin/bash" > /workspace/.devcontainer/entrypoint.sh \
 && cat /shell-hook >> /workspace/.devcontainer/entrypoint.sh \
 && echo 'exec "$@"' >> /workspace/.devcontainer/entrypoint.sh \
 && chmod +x /workspace/.devcontainer/entrypoint.sh

# Expose your app's port
EXPOSE 8000

# Activate Pixi env and run your app
ENTRYPOINT ["/workspace/.devcontainer/entrypoint.sh"]
CMD ["uvicorn", "my_project:app", "--host", "0.0.0.0"]

