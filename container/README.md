Build container:
```
sudo -E singularity build sc-adt-rna-10x.v6.sif sc-adt-rna-10x.v6-build
```

#### Run command with container
```
singularity exec --bind /fs1 sc-adt-rna-10x.v6.sif cellranger count (..arguments..)
```
