FROM python:3.12-slim-bookworm

WORKDIR /app

# ── Layer 1: OS dependencies ────────────────────────────────────────────────
# apt-get upgrade -y resolves the bulk of SCA/SCD CVE findings from Checkmarx
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    bash git wget curl unzip tar \
    build-essential libssl-dev libbz2-dev libreadline-dev libsqlite3-dev zlib1g-dev \
    sqlite3 r-base ncompress tcsh perl \
    fonts-urw-base35 gsfonts \
  && rm -rf /var/lib/apt/lists/*

# ── Layer 2: Python 2.7 (required by seq2logo and some legacy tool helpers) ─
# Debian bookworm dropped python2 from main repos — compile from source
RUN wget -q https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz \
  && tar xf Python-2.7.18.tgz \
  && cd Python-2.7.18 \
  && ./configure --prefix=/usr/local/python2 > /dev/null \
  && make -j"$(nproc)" > /dev/null \
  && make install > /dev/null \
  && cd .. && rm -rf Python-2.7.18 Python-2.7.18.tgz \
  && ln -sf /usr/local/python2/bin/python2 /usr/local/bin/python2

RUN wget -q https://bootstrap.pypa.io/pip/2.7/get-pip.py \
  && python2 get-pip.py > /dev/null && rm get-pip.py \
  && python2 -m pip install --quiet numpy matplotlib

# ── Layer 3: Ghostscript 9.53.3 pre-built binary ──────────────────────────
# System gs (9.06) has a font rendering bug — must use 9.53.3
RUN wget -q https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs9533/ghostscript-9.53.3-linux-x86_64.tgz \
  && tar -xzf ghostscript-9.53.3-linux-x86_64.tgz -C /app/app/tools/ \
  && rm ghostscript-9.53.3-linux-x86_64.tgz \
  && chmod +x /app/app/tools/ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64

# ── Layer 4: App code ───────────────────────────────────────────────────────
COPY . .

# ── Layer 5: seq2logo and GibbsCluster (tarballs in repo) ──────────────────
RUN tar -xzf app/tools/seq2logo-2.1.all.tar.gz -C app/tools/ \
  && rm app/tools/seq2logo-2.1.all.tar.gz

RUN tar -xvf app/tools/gibbscluster-2.0f.Linux.tar.gz -C app/tools/ \
  && rm app/tools/gibbscluster-2.0f.Linux.tar.gz

# Patch GibbsCluster GIBBS path
RUN sed -i 's|setenv GIBBS .*|setenv GIBBS /app/app/tools/gibbscluster-2.0|' \
    app/tools/gibbscluster-2.0/gibbscluster

# Patch seq2logo to use bundled Ghostscript 9.53.3
RUN sed -i "s|gsPath='gs'|gsPath='/app/app/tools/ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64'|" \
    app/tools/seq2logo-2.1/Seq2Logo.py

# Apply GibbsCluster SA script patches (output path + seqlogo command)
RUN sed -i \
    -e 's|^\(\s*\$resdir .= "/\$prefix";\)|# \1  # removed prefix from output path|' \
    -e 's|^\(my \$barplot = "\$resdir/images/\$prefix.gibbs.KLDvsCluster.barplot.png";\)|my \$barplot = "\$resdir/images/gibbsKLDvsCluster.barplot.JPG";|' \
    -e '530s@.*@$cmd .= "$seq2logo -f $corefile -o $logofile -I 2 --format [JPEG] -b $wlc -C 2 -S 2 -t $title \&>/dev/null; ";@' \
    app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl

# ── Layer 6: Licensed tools (netMHCpan / netMHCIIpan) ──────────────────────
# IMPORTANT: These tools require a DTU academic license.
# Download from https://services.healthtech.dtu.dk/software.php and place the
# tarballs in app/tools/ BEFORE running docker build.
# If the tarballs are absent the build continues and those prediction methods
# will be unavailable at runtime.
RUN if [ -f app/tools/netMHCpan-4.2c.Linux.tar.gz ]; then \
      gunzip -c app/tools/netMHCpan-4.2c.Linux.tar.gz | tar xf - -C app/tools/ \
      && mkdir -p app/tools/netMHCpan-4.2/tmp \
      && sed -i 's|setenv  *NMHOME  *[^ ]*|setenv NMHOME /app/app/tools/netMHCpan-4.2|g' \
             app/tools/netMHCpan-4.2/netMHCpan \
      && sed -i 's|setenv  *TMPDIR  *[^ ]*|setenv TMPDIR /app/app/tools/netMHCpan-4.2/tmp|g' \
             app/tools/netMHCpan-4.2/netMHCpan \
      && rm app/tools/netMHCpan-4.2c.Linux.tar.gz; \
    else echo "WARNING: netMHCpan-4.2c.Linux.tar.gz not found — netMHCpan will be unavailable"; fi

RUN if [ -f app/tools/netMHCIIpan-4.3j.Linux.tar.gz ]; then \
      tar -xvf app/tools/netMHCIIpan-4.3j.Linux.tar.gz -C app/tools/ \
      && mkdir -p app/tools/netMHCIIpan-4.3/tmp \
      && sed -i 's|setenv  *NMHOME  *[^ ]*|setenv NMHOME /app/app/tools/netMHCIIpan-4.3|g' \
             app/tools/netMHCIIpan-4.3/netMHCIIpan \
      && rm app/tools/netMHCIIpan-4.3j.Linux.tar.gz; \
    else echo "WARNING: netMHCIIpan-4.3j.Linux.tar.gz not found — netMHCIIpan will be unavailable"; fi

# ── Layer 7: MixMHCpred and MixMHC2pred ────────────────────────────────────
RUN wget -q https://github.com/GfellerLab/MixMHCpred/archive/refs/tags/v3.0.tar.gz -O /tmp/mixmhcpred.tar.gz \
  && tar -xzf /tmp/mixmhcpred.tar.gz -C app/tools/ \
  && mv app/tools/MixMHCpred-3.0 app/tools/MixMHCpred \
  && chmod +x app/tools/MixMHCpred/MixMHCpred \
  && rm /tmp/mixmhcpred.tar.gz

RUN wget -q https://github.com/GfellerLab/MixMHC2pred/releases/download/v2.0.2.2/MixMHC2pred-2.0.zip -O /tmp/mixmhc2pred.zip \
  && unzip -q /tmp/mixmhc2pred.zip -d app/tools/ \
  && chmod +x app/tools/MixMHC2pred-2.0/MixMHC2pred_unix \
  && rm /tmp/mixmhc2pred.zip

# ── Layer 8: HLA-PepClust (MHC-TP) ─────────────────────────────────────────
# Clone into HLA-PepClust/ — this directory name is hardcoded in app/utils.py
RUN git clone --depth 1 --branch immunolyser/class2-mhctp \
      https://github.com/PurcellLab/MHC-TP.git app/tools/HLA-PepClust \
  && cd app/tools/HLA-PepClust \
  && python3 -m venv hlapepclust-env \
  && hlapepclust-env/bin/pip install --quiet -e .

# ── Layer 9: Python 3 virtualenv + app dependencies ────────────────────────
ENV MHCFLURRY_DATA_PATH=/app/.mhcflurry

RUN python3 -m venv lenv \
  && lenv/bin/pip install --quiet --upgrade pip \
  && lenv/bin/pip install --quiet -r requirements_python3.txt \
  && lenv/bin/mhcflurry-downloads fetch

# Run any hotfix patching needed after package install
RUN lenv/bin/python hotfix_package_files.py 2>/dev/null || true

# ── Layer 10: Non-root user (fixes IaC/CON findings) ───────────────────────
RUN useradd -m --uid 1000 appuser \
  && chown -R appuser:appuser /app

USER appuser

# ── Healthcheck (fixes IaC finding) ─────────────────────────────────────────
HEALTHCHECK --interval=30s --timeout=10s --start-period=90s --retries=3 \
  CMD curl -f http://localhost:5000/healthz || exit 1

EXPOSE 5000

ENTRYPOINT ["./entrypoint.sh"]
