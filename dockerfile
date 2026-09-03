FROM ubuntu:22.04
ENV PATH="/root/.pixi/bin:${PATH}"
RUN mkdir -p scripts
COPY install.sh /scripts/
COPY joint_genotyping.sh ./
COPY merge_tables.py ./
COPY clinvar.toml ./
COPY clinvar.lua ./
RUN bash /scripts/install.sh
ENTRYPOINT ["bash"]
CMD ["bash", "joint_genotyping.sh", "--help"]