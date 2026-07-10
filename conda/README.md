# Conda Packaging

This directory contains the conda recipe for bioconda submission.

## Files

- `basevar/meta.yaml` - Package metadata and dependencies
- `basevar/build.sh` - Build script for conda

## Updating to a new version

1. Update `basevar/meta.yaml`:
   - Change `version` to new version
   - Update `sha256` with new release tarball hash

2. Reset `build.number` to 0

3. Copy to bioconda-recipes fork and submit PR:
   ```bash
   cd ~/Projects/bioconda-recipes
   git checkout master && git pull origin master
   git checkout -b update-basevar-<version>
   
   cp <this-repo>/conda/basevar/meta.yaml recipes/basevar/
   cp <this-repo>/conda/basevar/build.sh recipes/basevar/
   chmod +x recipes/basevar/build.sh
   
   git add recipes/basevar/
   git commit -m "Update basevar to v<version>"
   git push origin update-basevar-<version>
   ```

4. Create PR at https://github.com/bioconda/bioconda-recipes/compare

## References

- [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
- [bioconda contribution guide](https://bioconda.github.io/contributor/)
- [Conda build documentation](https://docs.conda.io/projects/conda-build/)
