# Conda Packaging

This directory contains a **local copy** of the conda recipe for reference.

> **Note**: The authoritative recipe is maintained at [bioconda-recipes/recipes/basevar](https://github.com/bioconda/bioconda-recipes/tree/main/recipes/basevar).
> Version updates are handled automatically by the **bioconda autobump bot** when a new GitHub Release is published.

## Files

- `basevar/meta.yaml` - Package metadata and dependencies
- `basevar/build.sh` - Build script for conda

## How updates work

**For standard version updates** (no recipe changes needed):
- Simply publish a new GitHub Release
- The autobump bot will automatically detect it, update the recipe, and create a PR
- No action needed in this directory

**For recipe changes** (adding dependencies, modifying build options, etc.):
1. Edit files in this directory
2. Copy to bioconda-recipes fork and submit PR:
   ```bash
   cd ~/Projects/bioconda-recipes
   git checkout master && git pull origin master
   git checkout -b update-basevar-<description>
   
   cp <this-repo>/conda/basevar/meta.yaml recipes/basevar/
   cp <this-repo>/conda/basevar/build.sh recipes/basevar/
   chmod +x recipes/basevar/build.sh
   
   git add recipes/basevar/
   git commit -m "<description>"
   git push origin update-basevar-<description>
   ```
3. Create PR at https://github.com/bioconda/bioconda-recipes/compare

## References

- [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
- [bioconda contribution guide](https://bioconda.github.io/contributor/)
- [Conda build documentation](https://docs.conda.io/projects/conda-build/)
