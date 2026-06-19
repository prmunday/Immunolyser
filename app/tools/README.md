# Bioinformatics Tools

This directory holds the bioinformatics prediction tools used by Immunolyser.

## Tools included in the repository (open-source)

The following tool archives are committed to this repo and extracted automatically during `docker build`:

| Tool | Archive |
|------|---------|
| seq2logo 2.1 | `seq2logo-2.1.all.tar.gz` |
| GibbsCluster 2.0f | `gibbscluster-2.0f.Linux.tar.gz` |

MixMHCpred, MixMHC2pred, and HLA-PepClust (MHC-TP) are downloaded from their public GitHub releases at build time.

## Licensed tools — you must download these before building

**netMHCpan 4.2** and **netMHCIIpan 4.3** require a free academic licence from DTU Health Tech.

1. Go to https://services.healthtech.dtu.dk/software.php
2. Register and accept the licence for each tool
3. Download the Linux tarballs:
   - `netMHCpan-4.2c.Linux.tar.gz`
   - `netMHCIIpan-4.3j.Linux.tar.gz`
4. Place them in this directory (`app/tools/`) before running `docker build`

The Dockerfile detects whether the tarballs are present and prints a warning (not an error) if they are missing. Without them, NetMHCpan-based predictions will not be available.

## Ref data for HLA-PepClust

The large reference databases (~1 GB each) are NOT in the repository and must be downloaded separately after `docker compose up`:

```bash
# Inside the running flask_app container:
docker compose exec flask_app bash
cd app/tools/HLA-PepClust/data/ref_data

# Class I human ref data
python3 -m gdown 'https://drive.google.com/uc?id=1iAAvir1woMOnURkP46zr_ETqpW2oUgGD'
unzip Gibbs_motifs_human.zip && rm Gibbs_motifs_human.zip

# Class II human ref data
python3 -m gdown 'https://drive.google.com/uc?id=1Dd9_63vCHezgculDycqfbCDzXsARsvbi'
unzip Gibbs_motifs_human_classII.zip && rm Gibbs_motifs_human_classII.zip

# Build the Class II database
source ../hlapepclust-env/bin/activate
python -m cli.main --database scripts/classII_config.json
deactivate
```

After the first Class II job completes, warm the Numba cache once (or MHC-TP output will be empty):
```bash
# Use the gibbscluster output dir from a completed Class II job:
clust-search /pvol/<taskId>/<sample>/gibbscluster/<replicate>/<subdir>/ \
  data/ref_data -im --output /tmp/cache_warmup \
  --NumbaDB data/ref_data -s human_classii -t 0.1 --processes 4
rm -rf /tmp/cache_warmup
```
