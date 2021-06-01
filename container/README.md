Build container:
```
sudo -E singularity build sc-cite-seq-10x.sif sc-cite-seq-10x-builder
```

#### Run command with container
```
singularity exec --bind /fs1 sc-cite-seq-10x.sif cellranger count (..arguments..)
```
