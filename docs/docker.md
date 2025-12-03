# Docker Usage Guide for VCFX

This document explains how to use VCFX with Docker.

## Using the Pre-built Image (Recommended)

VCFX is available as a pre-built Docker image on GitHub Container Registry:

```bash
# Pull the image (only needed once)
docker pull ghcr.io/jorgemfs/vcfx:latest

# Run a VCFX tool
docker run --rm ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]

# Mount a directory with your data
docker run --rm -v /path/to/your/data:/data ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]

# Example: Process a VCF file (using tests/data/valid.vcf as an example)
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest 'cat /data/valid.vcf | VCFX_allele_freq_calc > /data/output.tsv'
```

Using the pre-built image is recommended for most users as it:

- Requires no build time
- Is automatically updated with each release
- Has been tested and verified to work correctly

## Testing the Docker Image

A test script is provided to verify the Docker image works correctly with VCFX:

```bash
# Run the Docker test script
./tests/test_docker.sh
```

This script will:

1. Pull the latest Docker image
2. Test several VCFX tools using test data
3. Verify the tools work correctly in the Docker environment

## Building Your Own Image

If you need to customize the Docker image, you can build it yourself:

```bash
# Clone the repository
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX

# Using Docker directly
docker build -t vcfx:local .

# Using Docker Compose
docker-compose build
```

## Entrypoint Script and Passing Commands

The image uses `/usr/local/bin/docker_entrypoint.sh` as its entrypoint. This script
adds all VCFX tools to the `PATH` and then executes whatever command you pass to
`docker run`.

```bash
docker run --rm ghcr.io/jorgemfs/vcfx:latest VCFX_variant_counter < input.vcf
```

You can substitute `VCFX_variant_counter` with any other tool or quote a more
complex shell command.

## Running VCFX Tools

There are several ways to run VCFX tools with Docker:

### Using Docker Directly

```bash
# With the pre-built image
docker run --rm ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]

# With a locally built image
docker run --rm vcfx:local VCFX_tool_name [options]

# Mount the tests/data directory to access test files
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]

# Process files in the tests/data directory
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest 'cat /data/valid.vcf | VCFX_validator'

# Example: Calculate allele frequencies for a VCF file
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest 'cat /data/valid.vcf | VCFX_allele_freq_calc > /data/output.tsv'
```

### Using Docker Compose

```bash
# Basic usage (with locally built image)
docker-compose run --rm vcfx VCFX_tool_name [options]

# Example: List all available tools
docker-compose run --rm vcfx 'ls -1 /usr/local/bin/VCFX_*'

# Example: Process a VCF file from tests/data
docker-compose run --rm vcfx 'cat /data/valid.vcf | VCFX_allele_freq_calc > /data/output.tsv'
```

## Data Management

When using Docker directly, you need to mount a directory to access your files:

```bash
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]
```

When using Docker Compose, the `tests/data` directory is mounted by default:

1. VCF files in the tests/data directory are accessible in the container at `/data`
2. Output files will be saved back to the tests/data directory

You can modify the docker-compose.yml file to mount a different directory if needed.

## Advanced Usage

### Creating Pipelines

You can create complex pipelines by chaining VCFX tools:

```bash
docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest 'cat /data/classifier_mixed.vcf | VCFX_variant_classifier --append-info | grep "VCF_CLASS=SNP" | VCFX_allele_freq_calc > /data/snp_frequencies.tsv'
```

### Creating Shell Scripts

For complex workflows, consider creating a shell script:

```bash
#!/bin/bash
# save as vcfx_workflow.sh

docker run --rm -v $(pwd)/tests/data:/data ghcr.io/jorgemfs/vcfx:latest 'cat /data/valid.vcf | \
  VCFX_validator | \
  VCFX_variant_classifier --append-info | \
  VCFX_allele_freq_calc > /data/pipeline_output.tsv'
```

Then make it executable and run it:

```bash
chmod +x vcfx_workflow.sh
./vcfx_workflow.sh
```

## Troubleshooting

### Permission Issues

If you encounter permission issues with files created in the container:

```bash
# Run the container with your user ID
docker run --rm -v $(pwd)/tests/data:/data -u $(id -u):$(id -g) ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name [options]
```

### Container Not Finding Commands

If the container can't find VCFX commands, ensure they were properly built in the image:

```bash
# List available VCFX tools in the container
docker run --rm ghcr.io/jorgemfs/vcfx:latest 'ls -1 /usr/local/bin/VCFX_*'
```

## Citation

If you use VCFX with Docker in your research, please cite:

```bibtex
@inproceedings{silva2025vcfx,
  title={VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis},
  author={Silva, Jorge Miguel and Oliveira, Jos{\'e} Luis},
  booktitle={12th International Work-Conference on Bioinformatics and Biomedical Engineering (IWBBIO 2025)},
  year={2025},
  organization={Springer}
}
```

For more citation formats and information, see the [Citation page](citation.md). 
