# Releasing Pyorbital

1. checkout master
2. pull from repo
3. run the unittests
4. run `loghub` and update the `CHANGELOG.md` file:

```
loghub pytroll/pysorbital --token <personal access token (see https://github.com/settings/tokens)>  -st v<previous version> -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes"

```

Don't forget to commit!

5. Create a tag with the new version number, starting with a 'v', eg:

```
git tag -a v0.22.45 -m "Version 0.22.45"
```

See [semver.org](http://semver.org/) on how to write a version number.

6. push changes to github `git push --follow-tags`

7. Verify github actions passed and code deployed to PyPI
