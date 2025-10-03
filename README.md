# ViralSeq Utils

## Docker Images

Docker images are automatically built and published to GitHub Container Registry when new tags are created.

### Latest Release
```bash
docker pull ghcr.io/jonbra/viralseq_utils:latest
```

### Specific Version
```bash
docker pull ghcr.io/jonbra/viralseq_utils:v1.0.0
```

### Running the Container
```bash
docker run ghcr.io/jonbra/viralseq_utils:latest
```

## Development

### Building Locally
```bash
docker build -t viralseq_utils .
```

### Creating a Release
To trigger a new Docker image build:
1. Create and push a new tag:
   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```
2. The GitHub Actions workflow will automatically build and publish the image

### Legacy Docker Hub Image
~~Download docker image here: https://hub.docker.com/repository/docker/jonbra/create_samplesheet/general~~  
~~`docker pull jonbra/create_samplesheet:1.0`~~

*Note: This repository now uses GitHub Container Registry instead of Docker Hub.*
